import pandas as pd
import numpy as np
import xarray as xr


def calculate_hourly_space_heat_demand(total_nc, water_nc, space_nc):

    # This function caculates the hourly space heat demand
    
    # Space heat demand can't be lower than 0
    def get_space_heat(ndarray_1, ndarray_2):
        return max(0, ndarray_1 - ndarray_2)
    vect_get_space_heat = np.vectorize(get_space_heat)

    total_xarray = xr.open_dataset(total_nc)
    var_name1 = list(total_xarray.data_vars)[0]
    total_ndarray = total_xarray.data_vars[var_name1].values

    water_xarray = xr.open_dataset(water_nc)
    var_name2 = list(total_xarray.data_vars)[0]
    water_ndarray = water_xarray.data_vars[var_name2].values

    space_ndarray = []
    length = len(total_ndarray)
    for n in range(length):
        space_ndarray.append(vect_get_space_heat(total_ndarray[n], water_ndarray[n]))
    space_ndarray = np.array(space_ndarray)

    lons = total_xarray.coords['longitude']
    lats = total_xarray.coords['latitude']
    time = [f'h_{x}' for x in range(8760)]

    space_xarray = xr.DataArray(space_ndarray, coords=[time, lats, lons], name='hourly_factor')

    space_xarray.to_netcdf(space_nc)


def synthesize(total_nc_list, water_nc_list, space_nc_list):

    # Batch process
    length = len(total_nc_list)
    for n in range(length):
        calculate_hourly_space_heat_demand(total_nc_list[n], water_nc_list[n], space_nc_list[n])


if __name__ == "__main__":
    synthesize(
        total_nc_list = snakemake.input.total, 
        water_nc_list = snakemake.input.water, 
        space_nc_list = snakemake.output.space
        )