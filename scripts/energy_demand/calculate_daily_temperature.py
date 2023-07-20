import xarray as xr
import numpy as np
import pandas as pd


def calculate_raw_average_daily_temperature(input_total_ndarray):

    # This function calculates the average daily temperature from 2000 to 2018
    name = [i for i in input_total_ndarray.data_vars][0]
    input_total_ndarray = input_total_ndarray.data_vars[name]
    shape = input_total_ndarray.shape

    selected_ndarray = np.mean(input_total_ndarray[0: 24], axis=0)
    for i in range(1, 365):
        selected_ndarray = np.concatenate((selected_ndarray, np.mean(input_total_ndarray[24 * i: 24 * (i + 1)], axis=0)), axis=0)
    selected_ndarray = selected_ndarray.reshape(365, shape[1], shape[2]) - 273.15

    return selected_ndarray


def get_reference_daily_temperature(input_total_ndarray):

    # This function calculates the reference average daily temperature
    shape = input_total_ndarray.shape

    ndarray_0 = (input_total_ndarray[0] + 0.5 * input_total_ndarray[-1] + 0.25 * input_total_ndarray[-2] + \
        0.125 * input_total_ndarray[-3]) / (1 + 0.5 + 0.25 + 0.125)
    ndarray_1 = (input_total_ndarray[1] + 0.5 * input_total_ndarray[0] + 0.25 * input_total_ndarray[-1] + \
        0.125 * input_total_ndarray[-2]) / (1 + 0.5 + 0.25 + 0.125)
    ndarray_2 = (input_total_ndarray[2] + 0.5 * input_total_ndarray[1] + 0.25 * input_total_ndarray[0] + \
        0.125 * input_total_ndarray[-1]) / (1 + 0.5 + 0.25 + 0.125)

    wanted_ndarray = np.concatenate((ndarray_0, ndarray_1, ndarray_2), axis=0)
    for i in range(3, 365):
        wanted_ndarray = np.concatenate((wanted_ndarray, (input_total_ndarray[i] + 0.5 * input_total_ndarray[i-1] + 0.25 * input_total_ndarray[i-2] + \
            0.125 * input_total_ndarray[-3]) / (1 + 0.5 + 0.25 + 0.125)), axis=0)
    
    wanted_ndarray = wanted_ndarray.reshape(shape) 

    return wanted_ndarray


def calculate_daily_temperature(input_nc, output_nc):

    input_xarray = xr.open_dataset(input_nc)

    wanted_ndarray = get_reference_daily_temperature(
        calculate_raw_average_daily_temperature(input_xarray))

    lons = input_xarray.coords['longitude']
    lats = input_xarray.coords['latitude']
    time = [f'd_{x}' for x in range(365)]

    output_ndarray = xr.DataArray(wanted_ndarray, coords=[time, lats, lons], name = 'ref_temp')

    output_ndarray.to_netcdf(output_nc)


if __name__ == "__main__":
    calculate_daily_temperature(
        input_nc = snakemake.input[0], 
        output_nc = snakemake.output[0]
        )