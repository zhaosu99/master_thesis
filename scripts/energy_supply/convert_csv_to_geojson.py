import geopandas as gpd
import pandas as pd


def from_csv_to_geojson(input_csv, output_geojson):

    # This function converts the csv with coordinates to geodataframe
    gdf = pd.read_csv(input_csv)
    gdf = gpd.GeoDataFrame(gdf, geometry=gpd.points_from_xy(gdf.x_m, gdf.y_m, 3035))['geometry']
    gdf = gpd.GeoDataFrame({'geometry': gdf}, crs=3035)
    gdf = gdf.rename(columns={'geometry': 'coordinates'}).set_geometry('coordinates')
    gdf.to_file(output_geojson)


if __name__ == "__main__":
    from_csv_to_geojson(
        input_csv = snakemake.input[0], 
        output_geojson = snakemake.output[0]
        )