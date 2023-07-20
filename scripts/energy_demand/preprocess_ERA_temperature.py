import xarray as xr
import numpy as np


def get_average_nc(input_nc_list, output_nc):

    # This function calculates the average temperature from existing nc files
    wanted_array = np.array(xr.open_dataset(input_nc_list[0]).data_vars['t2m'])[: 8760]

    length = len(input_nc_list)
    for n in range(1, length):
        next_array = np.array(xr.open_dataset(input_nc_list[n]).data_vars['t2m'])[: 8760]
        wanted_array = (n * next_array + wanted_array) / (n + 1)

    output_array = xr.open_dataset(input_nc_list[0])
    lons = output_array.coords['longitude']
    lats = output_array.coords['latitude']
    time = output_array.coords['time'][: 8760]
    new_xarray = xr.DataArray(wanted_array[: 8760], coords=[time, lats, lons], name='ave_temp')

    new_xarray.to_netcdf(output_nc)


if __name__ == "__main__":
    get_average_nc(
        input_nc_list = snakemake.input, 
        output_nc = snakemake.output[0]
        )