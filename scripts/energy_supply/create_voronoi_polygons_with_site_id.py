import xarray as xr
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point
from shapely.ops import voronoi_diagram
from shapely.geometry import MultiPoint


def create_voronoi_polygons_with_site_id(input_nc, output_geojson):

    '''
    This function extracts the site_id and coordinates of the points in existing renewable.ninja simulations, 
    then create voronoi polygons for existing points.

    The created geodataframe has two columns: geometry and site_id, 
    where the geometry column is made of the voronoi polygons from existing site id coordinates.

    The created geodataframe is also cliped to the study scope
    '''

    raw_data = xr.open_dataset(input_nc)
    site_id_list, coord_list = [int(i) for i in raw_data.coords['site_id']], []
    for id in site_id_list:
        coord = Point(float(raw_data.sel(site_id=id).lon), float(raw_data.sel(site_id=id).lat))
        coord_list.append(coord)
    coord_list = voronoi_diagram(MultiPoint(coord_list))

    created_geodataframe = gpd.GeoDataFrame(
        {'geometry': coord_list, 'site_id': site_id_list}, geometry='geometry', crs=4326
        ).to_crs(3035)
    
    created_geodataframe.to_file(output_geojson)


if __name__ == "__main__":
    create_voronoi_polygons_with_site_id(
        input_nc = snakemake.input[0],
        output_geojson = snakemake.output[0]
        )