import pandas as pd
import numpy as np
import xarray as xr


def merge_windy_and_normal(wind_nc, windy_nc, normal_nc, output_nc):

    # This function merges the heat demand series in windy and normal conditions
    # If wind speed > 4.4, then adopt the value in windy condition; else, adopt the value in normal condition
    wind_xarray = xr.open_dataset(wind_nc)
    windy_xarray = xr.open_dataset(windy_nc)
    normal_xarray = xr.open_dataset(normal_nc)

    wind_ndarray = wind_xarray[[i for i in wind_xarray.data_vars][0]].values
    windy_ndarray = windy_xarray[[i for i in windy_xarray.data_vars][0]].values
    normal_ndarray = normal_xarray[[i for i in normal_xarray.data_vars][0]].values

    days_in_months = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    assert len(wind_ndarray) == 12

    wanted_array = []
    for n in range(len(wind_ndarray)):
        monthly_wind_ndarray = wind_ndarray[n]
        for m in range(sum(days_in_months[: n]), sum(days_in_months[: (n + 1)])):
            monthly_windy_ndarray = windy_ndarray[m]
            monthly_normal_ndarray = normal_ndarray[m]
            rows, cols = monthly_windy_ndarray.shape
            model = monthly_windy_ndarray.copy()
            for i in range(rows):
                for j in range(cols):
                    if monthly_wind_ndarray[i, j] <= 4.4:
                        model[i, j] = monthly_normal_ndarray[i, j]
            wanted_array.append(model)
    wanted_array = np.array(wanted_array)

    lons = normal_xarray.coords['longitude']
    lats = normal_xarray.coords['latitude']
    time = normal_xarray.coords['dim_0']

    xr.DataArray(wanted_array, coords=[time, lats, lons], name='merged_daily_factor').to_netcdf(output_nc)


if __name__ == "__main__":
    merge_windy_and_normal(
        wind_nc = snakemake.input.wind, 
        windy_nc = snakemake.input.water[0], 
        normal_nc = snakemake.input.water[3], 
        output_nc = snakemake.output.water[0]
        )
    merge_windy_and_normal(
        wind_nc = snakemake.input.wind, 
        windy_nc = snakemake.input.water[1], 
        normal_nc = snakemake.input.water[4], 
        output_nc = snakemake.output.water[1]
        )
    merge_windy_and_normal(
        wind_nc = snakemake.input.wind, 
        windy_nc = snakemake.input.water[2], 
        normal_nc = snakemake.input.water[5], 
        output_nc = snakemake.output.water[2]
        )
    merge_windy_and_normal(
        wind_nc = snakemake.input.wind, 
        windy_nc = snakemake.input.total[0], 
        normal_nc = snakemake.input.total[3], 
        output_nc = snakemake.output.total[0]
        )
    merge_windy_and_normal(
        wind_nc = snakemake.input.wind, 
        windy_nc = snakemake.input.total[1], 
        normal_nc = snakemake.input.total[4], 
        output_nc = snakemake.output.total[1]
        )
    merge_windy_and_normal(
        wind_nc = snakemake.input.wind, 
        windy_nc = snakemake.input.total[2], 
        normal_nc = snakemake.input.total[5], 
        output_nc = snakemake.output.total[2]
        )