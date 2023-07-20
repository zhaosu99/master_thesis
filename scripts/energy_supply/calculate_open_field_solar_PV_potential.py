from rasterstats import zonal_stats
import rasterio as rio
import geopandas as gpd
import pandas as pd


def get_regional_municipalities(municipal_geojson_with_site_id, region_geojson):

    # This function clips the national municipalitis to the boundary of selected region,
    # The results are regional municipalities

    municipal_geojson_with_site_id = gpd.read_file(municipal_geojson_with_site_id).to_crs(3035)
    region_geojson = gpd.read_file(region_geojson).to_crs(3035)

    # Here buffer(0) is to avoid self-intersection error
    try:
        regional_municipalities_gdf = gpd.clip(municipal_geojson_with_site_id, region_geojson)
    except:
        regional_municipalities_gdf = gpd.clip(municipal_geojson_with_site_id, region_geojson.buffer(0))

    return regional_municipalities_gdf


def zonal_stats_tif_on_atom_gdf(glaes_tif, atom_gdf, value_of_wanted_cells):

    # This function zonal_stats the raster on a vector
    # Suitable areas are added as a new column to the geodataframe
    suitable_area = rio.open(glaes_tif)
    cell_size_x = suitable_area.transform[0]
    cell_size_y = -suitable_area.transform[4]
    pixel_area = float(cell_size_x * cell_size_y / 10_000)

    stat = zonal_stats(atom_gdf, glaes_tif, categorical=True)

    suitable_area_list = []
    for n in range(len(stat)):
        try:
            suitable_area_list.append(
                sum([stat[n][i] * pixel_area for i in value_of_wanted_cells])
                )
        except:
            suitable_area_list.append(0)
            
    atom_gdf['suitable_area_ha'] = suitable_area_list

    atom_gdf_with_suitable_area = atom_gdf[atom_gdf['suitable_area_ha'] != 0]

    return atom_gdf_with_suitable_area


def add_capacity_factors(atom_gdf_with_suitable_area, capacity_factor_csv):

    # This function adds capacity factors to municipalities
    capacity_factor = pd.read_csv(capacity_factor_csv)
    atom_gdf_with_capacity_factors = atom_gdf_with_suitable_area.merge(capacity_factor, how="inner", on='site_id')

    return atom_gdf_with_capacity_factors


def summarize_suitable_area_on_atom_gdf(atom_gdf_with_capacity_factors):

    # This function summarizes the suitable areas and weighted average capacity factor of each municipality
    month_list = [str(i).zfill(2) for i in range(1, 13)]
    for i in month_list:
        atom_gdf_with_capacity_factors[i] *= atom_gdf_with_capacity_factors['suitable_area_ha']

    atom_gdf_with_capacity_factors = atom_gdf_with_capacity_factors.drop('site_id', axis=1)
    result = atom_gdf_with_capacity_factors.dissolve(by='municipal_id', aggfunc='sum')

    for i in month_list:
        result[i] /= result['suitable_area_ha']

    return result


def synthesize(national_municipalities_geojson, municipal_geojson_with_site_id, regional_geojson, glaes_tif, value_of_wanted_cells, capacity_factor_csv, output_geojson):
    
    # This function aggregates the above-mentioned functions
    # This purpose of creating a temporal index column is to guarantee the selected parts still remain to be a geodataframe rather than a geoseries
    regional_municipalities_gdf = get_regional_municipalities(municipal_geojson_with_site_id, regional_geojson)
    length = regional_municipalities_gdf.shape[0]
    regional_municipalities_gdf['index'] = [str(i) for i in range(length)]

    temp_list = []
    for n in range(length):
        try:
            atom_gdf = regional_municipalities_gdf[regional_municipalities_gdf['index']==str(n)].copy()     
            temp_list.append(
                add_capacity_factors(
                    zonal_stats_tif_on_atom_gdf(
                        glaes_tif, atom_gdf, value_of_wanted_cells
                        ), 
                    capacity_factor_csv) 
                )
        except:
            print(f"Error occurs when adding areas and capacity factors for the {n}th municipality")
    try:
        aggregated_result = summarize_suitable_area_on_atom_gdf(
            pd.concat(temp_list)
            )
    except:
        columns = ['geometry', 'municipal_id', 'suitable_area_ha'] + [str(i).zfill(2) for i in range(1, 13)]
        aggregated_result = gpd.GeoDataFrame(columns=columns, geometry='geometry')
    
    # Reset the geometry column
    aggregated_result = aggregated_result.drop('geometry', axis=1)
    geometry_gdf = gpd.read_file(national_municipalities_geojson).loc[:, ['geometry', 'municipal_id']]
    aggregated_result = geometry_gdf.merge(aggregated_result, on='municipal_id')
    aggregated_result.to_file(output_geojson)


if __name__ == "__main__":
    synthesize(
        national_municipalities_geojson = snakemake.input.national_municipalities_geojson,
        municipal_geojson_with_site_id = snakemake.input.municipal_geojson_with_site_id, 
        regional_geojson = snakemake.input.region,
        glaes_tif = snakemake.input.open_field_solar_PV_areas,
        value_of_wanted_cells = snakemake.params.value_of_wanted_cells, 
        capacity_factor_csv = snakemake.input.capacity_factor_csv, 
        output_geojson = snakemake.output[0]
        )

