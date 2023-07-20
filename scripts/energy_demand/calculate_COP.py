import xarray as xr
import numpy as np


def calculate_air_source_heat_pump_cop(x, type):

    assert type in ['radiator', 'floor', 'water']

    celsius = x - 273.15
    
    if type == 'radiator':
        t_sink = 40 - celsius
    elif type == 'floor':
        t_sink = 30 - 0.5 * celsius
    elif type == 'water':
        t_sink = 50

    delta_t = max(t_sink - celsius, 15)

    # This function comes from When2heat paper
    cop = 6.08 - 0.09 * delta_t + 0.0005 * delta_t ** 2

    return cop


def calculate_cop_for_each_point(input_nc, output_nc, type):

    # Apply the above-mentioned function to each existing point
    input_xarray = xr.open_dataset(input_nc)

    lons = input_xarray.coords['longitude']
    lats = input_xarray.coords['latitude']
    time = input_xarray.coords['time']
    array = input_xarray.data_vars['ave_temp']

    vectorized_function = np.vectorize(calculate_air_source_heat_pump_cop)
    wanted_array = vectorized_function(array[0], type)

    for n in range(1, len(array)):
        wanted_array = np.concatenate(
            (wanted_array, vectorized_function(array[n], type)), axis=0
            )
    wanted_array = wanted_array.reshape(array.shape)

    xr.DataArray(wanted_array, coords=[time, lats, lons], name=f'{type}_cop').to_netcdf(output_nc)


if __name__ == "__main__":
    calculate_cop_for_each_point(
        input_nc = snakemake.input[0], 
        output_nc = snakemake.output[0],
        type = 'radiator'
        )
    calculate_cop_for_each_point(
        input_nc = snakemake.input[0], 
        output_nc = snakemake.output[1],
        type = 'water'
        )
    