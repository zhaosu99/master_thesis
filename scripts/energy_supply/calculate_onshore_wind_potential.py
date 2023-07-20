import geopandas as gpd
import pandas as pd


def calculate_onshore_wind_potential(turbine_points_geojson, municipal_geojson_with_site_id, municipal_geojson, capacity_factor_csv, output_geojson):

    # This function calculates the average capacity factor and onshore turbine amount for each municipality
    input_points = gpd.read_file(turbine_points_geojson).to_crs(3035)
    municipal_geojson_with_site_id = gpd.read_file(municipal_geojson_with_site_id).to_crs(3035)
    municipal_geojson = gpd.read_file(municipal_geojson).to_crs(3035).loc[:, ['geometry', 'municipal_id']]

    # Add the municipal_id and site_id of nearest point from already-made renewable.ninja simulations
    output = input_points.sjoin(municipal_geojson_with_site_id, how="inner", predicate='intersects')
    output = output.loc[:, ['municipal_id', 'site_id']]

    # Calculate the turbine amount in each (municipal_id, site_id) pair
    output['turbine_amount'] = 1
    output = output.groupby(['municipal_id', 'site_id'])['turbine_amount'].sum().reset_index()
    output = output.loc[:, ['municipal_id', 'site_id', 'turbine_amount']]

    # Add the capacity factors in 12 months to each (municipal_id, site_id) pair
    capacity_factor = pd.read_csv(capacity_factor_csv)
    output = output.merge(capacity_factor, how="inner", on='site_id') 

    # Calculate the average capacity factor and total turbine amount for each municipal_id (municipal)
    month_list = [str(i).zfill(2) for i in range(1, 13)]
    for i in month_list:
        output[i] *= output['turbine_amount']
    output = output.drop(['site_id'], axis=1).groupby(['municipal_id']).sum().reset_index()
    for i in month_list:
        output[i] /= output['turbine_amount']

    # Add the geometry of each municipal_id (municipal)
    output = municipal_geojson.merge(output, how="inner", on="municipal_id").dropna()

    output.to_file(output_geojson)


if __name__ == "__main__":
    calculate_onshore_wind_potential(
        turbine_points_geojson = snakemake.input.onshore_turbines, 
        municipal_geojson_with_site_id = snakemake.input.municipal_geojson_with_site_id, 
        municipal_geojson = snakemake.input.national_municipalities,
        capacity_factor_csv = snakemake.input.capacity_factor_csv, 
        output_geojson = snakemake.output[0]
        )