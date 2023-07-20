from rasterstats import zonal_stats
import rasterio as rio
import geopandas as gpd


def zonal_stats_tif_on_municipality_gdf(railway_tif, municipality_geojson, output_geojson):

    # This function calculates the number of wanted raster pixel values on vectors, and convert the number to ratios
    with rio.open(railway_tif) as suitable_area:
        assert suitable_area.crs == 3035
    municipality_gdf = gpd.read_file(municipality_geojson)

    stat = zonal_stats(municipality_gdf, railway_tif, categorical=True)

    suitable_area_list = []
    for m in range(len(stat)):
        existing_value_list = []
        try:
            existing_value_list.append(stat[m][0])
        except:
            existing_value_list.append(0)
        suitable_area_list.append(sum(existing_value_list))

    municipality_gdf['railway_length_ratio'] = suitable_area_list
    column_sum = sum(suitable_area_list)

    if column_sum != 0:
        municipality_gdf['railway_length_ratio'] = municipality_gdf['railway_length_ratio'] / column_sum
    else:
        municipality_gdf['railway_length_ratio'] = 0

    municipality_gdf.to_file(output_geojson)


if __name__ == "__main__":
    zonal_stats_tif_on_municipality_gdf(
        railway_tif = snakemake.input[1], 
        municipality_geojson = snakemake.input[0], 
        output_geojson = snakemake.output[0]
        )