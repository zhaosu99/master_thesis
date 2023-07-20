from rasterstats import zonal_stats
import rasterio as rio
import geopandas as gpd


def zonal_stats_tif_on_municipality_gdf(demograph_tif, municipality_geojson, output_geojson):

    # This function adds the population ratio to each municipality
    with rio.open(demograph_tif) as suitable_area:
        assert suitable_area.crs == 3035
    municipality_gdf = gpd.read_file(municipality_geojson)

    stat = zonal_stats(municipality_gdf, demograph_tif, stats=['sum'])

    suitable_area_list = []
    for m in range(len(stat)):
        try:
            s = stat[m]['sum']
            suitable_area_list.append(s)
        except:
            suitable_area_list.append(0)

    municipality_gdf['population_ratio'] = suitable_area_list

    municipality_gdf['population_ratio'] = \
        municipality_gdf['population_ratio'] / municipality_gdf['population_ratio'].sum()

    municipality_gdf.to_file(output_geojson)


if __name__ == "__main__":
    zonal_stats_tif_on_municipality_gdf(
        demograph_tif = snakemake.input[1], 
        municipality_geojson = snakemake.input[0], 
        output_geojson = snakemake.output[0]
        )