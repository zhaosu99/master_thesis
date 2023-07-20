from rasterstats import zonal_stats
import rasterio as rio
import geopandas as gpd


def zonal_stats_tif_on_municipality_gdf(clc_tif, municipality_geojson, value_list_of_wanted_cells, area_name, output_geojson):

    # This function adds the ratio of certain type(s) of land cover to each municipality
    with rio.open(clc_tif) as suitable_area:
        assert suitable_area.crs == 3035
        
    municipality_gdf = gpd.read_file(municipality_geojson)

    stat = zonal_stats(municipality_gdf, clc_tif, categorical=True)

    suitable_area_list = []
    for m in range(len(stat)):
        existing_value_list = []
        for n in value_list_of_wanted_cells:
            try:
                existing_value_list.append(stat[m][n])
            except:
                continue
        suitable_area_list.append(sum(existing_value_list))

    municipality_gdf[area_name] = suitable_area_list
    column_sum = sum(suitable_area_list)

    if column_sum != 0:
        municipality_gdf[area_name] = municipality_gdf[area_name] / column_sum
    else:
        municipality_gdf[area_name] = 0

    municipality_gdf.to_file(output_geojson)


if __name__ == "__main__":
    zonal_stats_tif_on_municipality_gdf(
        clc_tif = snakemake.input[1], 
        municipality_geojson = snakemake.input[0], 
        value_list_of_wanted_cells = snakemake.params.target_cell_values, 
        area_name = snakemake.params.area_name, 
        output_geojson = snakemake.output[0]
        )