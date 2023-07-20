import pandas as pd
import numpy as np
import xarray as xr


def preprocess_csv(input_csv):

    # The original decimal points are ',', which need to be replaced by '.'
    input_df = pd.read_csv(input_csv, sep=';')
    cols = input_df.columns
    for i in cols:
        try:
            input_df[i] = input_df[i].apply(lambda x: float(x.replace(',', '.')))
        except:
            continue

    return input_df


def get_hourly_factors(input_df, daily_temperature_array, daily_demand_array, building_type, weekday=None):

    # This function expands the daily demand factors to hourly resolution using the hourly factors
    # daily_temperature_array and daily_demand_array are only for one day

    def get_column_name(x):
        if x <= -10:
            return '-15'
        elif -10 <= x < -5:
            return '-10'
        elif -5 <= x < 0:
            return '-5'
        elif 0 <= x < 5:
            return '0'
        elif 5 <= x < 10:
            return '5'
        elif 10 <= x < 15:
            return '10'
        elif 15 <= x < 20:
            return '15'
        elif 20 <= x < 25:
            return '20'
        elif 25 <= x < 30:
            return '25'
        elif x >= 30:
            return '30'
    
    # Here the array is reshaped to one row to facilitate analysis
    x, y = daily_temperature_array.shape
    temperature_one_row = daily_temperature_array.reshape(np.prod(daily_temperature_array.shape))
    col_name_list = [get_column_name(i) for i in temperature_one_row]
    daily_demand_one_row = daily_demand_array.reshape(np.prod(daily_demand_array.shape))

    hourly_factors_list = []
    for n, i in enumerate(col_name_list):   
        if building_type in ['SFH', 'MFH']:
            hourly_factors_list.append(
                daily_demand_one_row[n] * input_df.loc[:, i]) 
        elif building_type == 'COM':
            hourly_factors_list.append(
                daily_demand_one_row[n] * input_df.loc[lambda x: x.weekday == weekday, i])
        else:
            print('input error')

    # After processing, the array is reshaped back 
    hourly_factors_list = np.transpose(hourly_factors_list, (1, 0)).reshape(24, x, y)

    return hourly_factors_list


def get_hourly_heat_demand(input_csv, input_demand_nc, building_type, input_temperature_nc):

    # This function aggregates the hourly demand factors of each day to one year

    assert building_type in ['SFH', 'MFH', 'COM']

    input_df = preprocess_csv(input_csv)

    input_demand_xarray = xr.open_dataset(input_demand_nc)
    var_name1 = list(input_demand_xarray.data_vars)[0]
    total_daily_demand_array = input_demand_xarray.data_vars[var_name1].values

    input_temperature_xarray = xr.open_dataset(input_temperature_nc)
    var_name2 = list(input_temperature_xarray.data_vars)[0]
    total_daily_temp_array = input_temperature_xarray.data_vars[var_name2].values

    length = len(total_daily_demand_array) 
    assert length == 365

    # In the iteration, each daily factor number is multiplied with an array with all the hourly factors
    # To accelerate processing, 2darrays are iteratively concatenated
    # First 2darray
    if building_type in ['SFH', 'MFH']:
        result_ndarray = (
            get_hourly_factors(input_df, total_daily_temp_array[0], total_daily_demand_array[0], building_type))
    else:
        result_ndarray = (
            get_hourly_factors(input_df, total_daily_temp_array[0], total_daily_demand_array[0], building_type, 0))

    # Begin to iterate
    for n in range(1, length):
        if building_type in ['SFH', 'MFH']:
            next_ndarray = (
                get_hourly_factors(input_df, total_daily_temp_array[n], total_daily_demand_array[n], building_type)
                )
        else:
            next_ndarray = (
                get_hourly_factors(input_df, total_daily_temp_array[n], total_daily_demand_array[n], building_type, n % 7)
                )

        result_ndarray = np.concatenate((result_ndarray, next_ndarray), axis=0)

    lons = input_demand_xarray.coords['longitude']
    lats = input_demand_xarray.coords['latitude']
    time = [f'h_{x}' for x in range(8760)]

    result_xarray = xr.DataArray(result_ndarray, coords=[time, lats, lons], name='hourly_factor')

    return result_xarray


def synthesize(input_temperature_nc, input_water_nc_list, input_total_nc_list, input_hourly_factor_csv_list,
        output_hourly_water_nc_list, output_hourly_total_nc_list, building_list):

    # Batch calculate residential and commercial water and total heat demand
    for n in range(len(building_list)):

        total_xarray = get_hourly_heat_demand(input_hourly_factor_csv_list[n], input_total_nc_list[n], building_list[n], input_temperature_nc)
        water_xarray = get_hourly_heat_demand(input_hourly_factor_csv_list[n], input_water_nc_list[n], building_list[n], input_temperature_nc)

        total_xarray.to_netcdf(output_hourly_total_nc_list[n])
        water_xarray.to_netcdf(output_hourly_water_nc_list[n])


if __name__ == "__main__":
    synthesize(
        input_temperature_nc = snakemake.input.temperature, 
        input_water_nc_list = snakemake.input.water, 
        input_total_nc_list = snakemake.input.total, 
        input_hourly_factor_csv_list = snakemake.input.hourly_factor,
        output_hourly_water_nc_list = snakemake.output.water, 
        output_hourly_total_nc_list = snakemake.output.total,
        building_list = snakemake.params.building_types
        )