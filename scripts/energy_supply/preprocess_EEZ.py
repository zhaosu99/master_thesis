import geopandas as gpd
import pandas as pd


def select_eez(input_geojson, coastline_geojson, country_list, max_distance_to_land, outoput_geojson_list):

    # This function preprocess raw EEZ data to the study scope

    # Simplify the data, exclude useless columns
    studied_eezs = gpd.read_file(input_geojson).to_crs(3035)
    studied_eezs = studied_eezs[studied_eezs['ISO_Ter1'].isin(country_list)].loc[:, ['geometry', 'fid', 'ISO_Ter1']]

    # Exlcude overly far EEZ areas
    coastline_geojson = gpd.read_file(coastline_geojson).to_crs(3035)
    coastline_geojson = coastline_geojson.buffer(distance=max_distance_to_land, resolution=1).unary_union
    studied_eezs = studied_eezs.clip(coastline_geojson)

    # Generate qualified EEZ areas for each country
    length = len(country_list)
    for n in range(length):
        selected_eez = studied_eezs[studied_eezs['ISO_Ter1']==country_list[n]]
        selected_eez.to_file(outoput_geojson_list[n])


if __name__ == "__main__":
    select_eez(
        input_geojson = snakemake.input[0], 
        coastline_geojson = snakemake.input[1],
        country_list = snakemake.params.country_codes, 
        max_distance_to_land = snakemake.params.max_distance_to_land,
        outoput_geojson_list = snakemake.output.eez_by_country
        )