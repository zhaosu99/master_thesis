import geopandas as gpd
import pandas as pd


def concat_national_data(input_geojson_list, output_geojson):

    input_gdf_list = [gpd.read_file(i).to_crs(3035) for i in input_geojson_list]
    raw_output_gdf = pd.concat(input_gdf_list)
    if 'suitable_area_ha' in raw_output_gdf.columns:
        # Because in the last steps, buffer(0) was used to avoid self-intersection error, but this process lowers the accuracy of shapes
        # Thus, some small edges resulted from buffer(0) should be excluded now
        raw_output_gdf = raw_output_gdf.loc[lambda x: x.suitable_area_ha >= 0.5, :]
    raw_output_gdf = raw_output_gdf.reset_index(drop=True)

    raw_output_gdf.to_file(output_geojson)
    

if __name__ == "__main__":
    concat_national_data(
        input_geojson_list = snakemake.input, 
        output_geojson = snakemake.output[0]
        )