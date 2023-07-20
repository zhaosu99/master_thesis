import numpy as np
import xarray as xr
import pandas as pd
import geopandas as gpd


def get_sum_2darray(input_xarray):

    # Summarize daily data to a single average one
    var_name = list(input_xarray.data_vars)[0]
    total_ndarray = input_xarray.data_vars[var_name].values

    # In the last step, the values are multiplied with 1000000 to increase accuracy, now it's deducted again
    ratio_2darray = np.sum(total_ndarray, axis=0) / 1000000
    return ratio_2darray


def get_annual_heat_electricity_ratio(input_nc_list, output_nc, output_geojson):

    # Batch-process nc files, save all ratios into one netcdf dataset ile
    name_list = [f'{x}_{y}' for x in ['residental', 'commercial'] for y in ['radiator', 'water']]

    result_dataset = xr.Dataset({})

    for n in range(len(input_nc_list)):
        selected_xarray = xr.open_dataset(input_nc_list[n])  
        lons = selected_xarray.coords['longitude']
        lats = selected_xarray.coords['latitude']
        name = name_list[n]
        ratio_2darray = get_sum_2darray(selected_xarray)
        result_dataset[name] = xr.DataArray(ratio_2darray, coords=[lats, lons], dims=['latitude', 'longitude'], name=f'{name}')

    result_dataset.to_netcdf(output_nc)
    
    # Convert nc file to geojson file
    df = result_dataset.to_dataframe().reset_index()
    gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.longitude, df.latitude), crs=4326).to_crs(3035)
    cols = ['geometry'] + [i for i in list(gdf.columns) if i not in ['latitude', 'longitude', 'geometry']]
    gdf = gdf.loc[:, cols]
    gdf.to_file(output_geojson)


if __name__ == "__main__":
    get_annual_heat_electricity_ratio(
        input_nc_list = snakemake.input, 
        output_nc = snakemake.output[0],
        output_geojson = snakemake.output[1]
        )