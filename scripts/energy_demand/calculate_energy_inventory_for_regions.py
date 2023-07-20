import geopandas as gpd
import pandas as pd


def groupby_region(region_geojson, municipality_geojson, output_geojson, output_csv):

    # This function summarizes the energy demand or supply inventory of each region

    # Convert the ploygons in geometry columns into centroid points
    def polygon_to_centroid(input_geojson):
        input_gdf = gpd.read_file(input_geojson)
        input_gdf['geometry'] = input_gdf['geometry'].centroid
        return input_gdf

    municipality_gdf = polygon_to_centroid(municipality_geojson)

    region_gdf = gpd.read_file(region_geojson)
    result_df = municipality_gdf.sjoin(region_gdf, how='inner', predicate='intersects')
    result_df = result_df.drop(['geometry', 'index_right'], axis=1)
    result_df = result_df.groupby(['regional_id']).sum()

    result_gdf = pd.merge(region_gdf, result_df, on='regional_id')
    result_gdf.to_file(output_geojson)

    result_gdf.drop('geometry', axis=1).to_csv(output_csv)



if __name__ == "__main__":
    groupby_region(
        region_geojson = snakemake.input[0], 
        municipality_geojson = snakemake.input[1], 
        output_geojson = snakemake.output[0], 
        output_csv = snakemake.output[2]
        )
    groupby_region(
        region_geojson = snakemake.input[0], 
        municipality_geojson = snakemake.input[2], 
        output_geojson = snakemake.output[1], 
        output_csv = snakemake.output[3]
        )