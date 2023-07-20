import xarray as xr
import numpy as np


def get_monthly_average_wind(input_nc, output_nc):

    # This function calculates the average wind speed from existing nc files
    input_xarray = xr.open_dataset(input_nc)
    lons = input_xarray.coords['longitude']
    lats = input_xarray.coords['latitude']
    time = [f'm_{x}' for x in range(12)]
    array = input_xarray['si10'].values

    sum_array = []
    for k in range(12):
        relevent_years = [i * 12 + k for i in range(int(len(array) / 12))]
        sum_array.append(np.mean(np.array([array[n] for n in relevent_years]), axis=0))
    sum_array = np.array(sum_array)

    xr.DataArray(sum_array, coords=[time, lats, lons], name='monthly_average_wind').to_netcdf(output_nc)


if __name__ == "__main__":
    get_monthly_average_wind(
        input_nc = snakemake.input[0], 
        output_nc = snakemake.output[0]
        )

