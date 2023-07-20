import geopandas as gpd
import pandas as pd
from shapely.geometry import Point
import numpy as np


def convert_csv_to_geodataframe(input_csv):

    # This function converts the csv with coordinates to a geodataframe
    # The prerequisite of this function is to make sure the coordinates are in EPSG: 3035
    gdf = pd.read_csv(input_csv).iloc[:, 1:]
    gdf = gdf.apply(lambda n: Point(n.x_m, n.y_m), axis=1)
    gdf = pd.DataFrame(gdf).set_axis(['geometry'], axis=1)
    gdf = gpd.GeoDataFrame(gdf, geometry='geometry', crs=3035)

    return gdf


def add_info_to_points(input_csv, site_id_geojson, municipal_geojson, capacity_factor_csv, output_geojson):

    # This function calculates the average capacity factor and offshore turbine amount for each municipality
    input_points = convert_csv_to_geodataframe(input_csv)
    
    site_id_geojson = gpd.read_file(site_id_geojson).to_crs(3035)
    municipal_geojson = gpd.read_file(municipal_geojson).to_crs(3035).loc[:, ['municipal_id', 'geometry']]

    # Add the site ID of nearest point from already-made renewable.ninja simulations
    output = input_points.sjoin(site_id_geojson, how="inner", predicate='intersects')
    output = output.loc[:, ['geometry', 'site_id']]
    
    # Add the municipal_id of nearest municipal
    output = output.sjoin_nearest(municipal_geojson, how="inner")
    output = output.loc[:, ['municipal_id', 'site_id']]

    # Calculate the turbine amount in each (municipal_id, site_id) pair
    output['turbine_amount'] = 1
    output = output.groupby(['municipal_id', 'site_id'])['turbine_amount'].sum()

    # Reorganize the data
    output = output.reset_index()
    output = output.loc[:, ['municipal_id', 'site_id', 'turbine_amount']]

    # Add the capacity factors in 12 months to each (municipal_id, site_id) pair
    capacity_factor = pd.read_csv(capacity_factor_csv)
    output = output.merge(capacity_factor, how="inner", on='site_id')

    # Calculate the average capacity factor and total turbine amount for each municipality
    month_list = [str(i).zfill(2) for i in range(1, 13)]
    for i in month_list:
        output[i] *= output['turbine_amount']
        
    output = output.drop(['site_id'], axis=1).groupby(['municipal_id']).sum().reset_index()

    for i in month_list:
        output[i] /= output['turbine_amount']

    # Add the geometry of each municipal_id
    output = municipal_geojson.merge(output, how="inner", on="municipal_id")   

    output.to_file(output_geojson)


def create_empty_geodataframe(output_geojson):

    # Inland countries don't have offshore wind turbine sites, to avoid error, this functiom creates empty dataframe for them
    month_list = [str(i).zfill(2) for i in range(1, 13)]
    other_columns = ['geometry', 'municipal_id', 'turbine_amount']
    output = gpd.GeoDataFrame(columns=other_columns + month_list, geometry='geometry')

    output.to_file(output_geojson)


def calculate_offshore_wind_potential(input_csv, site_id_geojson, municipal_geojson, capacity_factor_csv, output_geojson):
    try:
        add_info_to_points(input_csv, site_id_geojson, municipal_geojson, capacity_factor_csv, output_geojson)
    except:
        create_empty_geodataframe(output_geojson)


if __name__ == "__main__":
    calculate_offshore_wind_potential(
        input_csv = snakemake.input.offshore_turbines, 
        site_id_geojson = snakemake.input.site_id_geojson, 
        municipal_geojson = snakemake.input.municipal_geojson, 
        capacity_factor_csv = snakemake.input.capacity_factor_csv, 
        output_geojson = snakemake.output[0]
        )