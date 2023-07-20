import geopandas as gpd
import pandas as pd


def concat_geojson_list(input_geojson_list, output_geojson):

    input_gdf_list = [gpd.read_file(i).to_crs(3035) for i in input_geojson_list]
    pd.concat(input_gdf_list).reset_index(drop=True).to_file(output_geojson)


if __name__ == "__main__":
    concat_geojson_list(
        input_geojson_list = snakemake.input, 
        output_geojson = snakemake.output[0]
        )