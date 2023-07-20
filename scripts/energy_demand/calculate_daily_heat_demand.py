import pandas as pd
import numpy as np
import xarray as xr
import geopandas as gpd
from pathlib import Path


def calculate_daily_heat_demand(input_nc, daily_demand_csv, temp_correction_csv, output_total_list, output_water_list):

    # This function calculates the daily heat total and water demand for single-family, multi-family and commercial houses
    # in windy and normal weather
    def get_bdew_params(input_df, type):

        # Extract parameters of BDEW heat demand function
        wanted_parameters = [float(str(i).replace(',', '.')) for i in list(input_df.loc[1: 8, type])]

        A = wanted_parameters[0]
        B = wanted_parameters[1]
        C = wanted_parameters[2]
        D = wanted_parameters[3]
        m_s = wanted_parameters[4]
        b_s = wanted_parameters[5]
        m_w = wanted_parameters[6]
        b_w = wanted_parameters[7]

        return A, B, C, D, m_s, b_s, m_w, b_w

    def total_function(A, B, C, D, m_s, b_s, m_w, b_w, delta_T, array):

        # Calculate total heat demand using BDEW function
        # The function needs to be vectorized, since it needs to process array rather than a number
        sigmoid = A / (1 + (B / (array - 40 - delta_T)) ** C) + D
        linear = max(m_s * array + b_s, m_w * array + b_w)

        return sigmoid + linear
    vect_total_function = np.vectorize(total_function)

    def water_function(D, m_w, b_w, array):

        # Calculate total heat demand using BDEW function
        # The function needs to be vectorized, since it needs to process array rather than a number
        if array > 15:
            return D + m_w * array + b_w
        else:
            return D + m_w * 15 + b_w
    vect_water_function = np.vectorize(water_function)


    code3 = Path(input_nc).stem[: 3]
    daily_demand_df = pd.read_csv(daily_demand_csv, sep=';')

    temp_correction_df = pd.read_csv(temp_correction_csv)
    delta_T = float(temp_correction_df.loc[lambda x: x.country == code3, 'heating_threshold_temperature']) - \
        float(temp_correction_df.loc[lambda x: x.country == 'DEU', 'heating_threshold_temperature'])

    input_xarray = xr.open_dataset(input_nc)
    lons = input_xarray.coords['longitude']
    lats = input_xarray.coords['latitude']
    time = input_xarray.coords['dim_0']
    array = input_xarray.data_vars['ref_temp']

    type_list = ['SFH', 'MFH', 'COM', 'SFH.1', 'MFH.1', 'COM.1']
    for n in range(len(type_list)):
        A, B, C, D, m_s, b_s, m_w, b_w = get_bdew_params(daily_demand_df, type_list[n])
        total_array = vect_total_function(A, B, C, D, m_s, b_s, m_w, b_w, delta_T, array)
        xr.DataArray(total_array, coords=[time, lats, lons], name='total_daily_heat_{}'.format(type_list[n])).to_netcdf(output_total_list[n])

        water_array = vect_water_function(D, m_w, b_w, array)
        xr.DataArray(water_array, coords=[time, lats, lons], name='water_daily_heat_{}'.format(type_list[n])).to_netcdf(output_water_list[n])



if __name__ == "__main__":
    calculate_daily_heat_demand(
        input_nc = snakemake.input[0], 
        daily_demand_csv = snakemake.input[1], 
        temp_correction_csv = snakemake.input[2], 
        output_total_list = snakemake.output.total, 
        output_water_list = snakemake.output.water
        )
