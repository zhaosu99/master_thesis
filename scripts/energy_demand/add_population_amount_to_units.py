from rasterstats import zonal_stats
import rasterio as rio
import geopandas as gpd


def zonal_stats_tif_on_unit_gdf(demograph_tif, unit_geojson, output_geojson):

    # This functuin zonal_stats the sum of raster pixel values on vectors
    with rio.open(demograph_tif) as suitable_area:
        assert suitable_area.crs == 3035
    unit_gdf = gpd.read_file(unit_geojson)

    stat = zonal_stats(unit_gdf, demograph_tif, stats=['sum'])

    suitable_area_list = []
    for m in range(len(stat)):
        try:
            s = stat[m]['sum']
            suitable_area_list.append(s)
        except:
            suitable_area_list.append(0)

    unit_gdf['population_amount'] = suitable_area_list

    unit_gdf.to_file(output_geojson)


if __name__ == "__main__":
    zonal_stats_tif_on_unit_gdf(
        demograph_tif = snakemake.input[1], 
        unit_geojson = snakemake.input[0], 
        output_geojson = snakemake.output[0]
        )