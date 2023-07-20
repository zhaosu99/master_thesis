import geopandas as gpd
import rasterio as rio
import pandas as pd
from rasterstats import zonal_stats


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


def return_tif_bounds_of_point(x, y):

    # This function returns the relevant ESM bounds of one point
    x, y = int(x), int(y)
    min_x = int((x - 44000)/10 ** 5) * 10 ** 5 + 44000
    min_y = int((y - 42000) / 10 ** 5) * 10 ** 5 + 42000
    max_x = (int((x - 44000) / 10 ** 5) + 1) * 10 ** 5 + 44000
    max_y = (int((y - 42000) / 10 ** 5) + 1) * 10 ** 5 + 42000
    output_tuple = (min_x, min_y, max_x, max_y)

    return output_tuple


def calculate_bounds_of_relevant_tifs_for_gdf(examined_gdf):

    # This function returns the relevant ESM bounds list of the geojson bounds points
    bounds = list(examined_gdf.bounds.iloc[0, :])

    bottom_left = return_tif_bounds_of_point(bounds[0], bounds[1])
    top_left = return_tif_bounds_of_point(bounds[0], bounds[3])
    bottom_right = return_tif_bounds_of_point(bounds[2], bounds[1])
    top_right = return_tif_bounds_of_point(bounds[2], bounds[3])

    x_min = min(bottom_left[0], top_left[0])
    y_min = min(bottom_left[1], bottom_right[1])
    x_max = max(bottom_right[2], top_right[2])
    y_max = max(top_left[3], top_right[3])

    # Get the coordinates of bottom left points of related tifs
    x_list = [i for i in range(x_min, x_max, 10 ** 5)]
    y_list = [i for i in range(y_min, y_max, 10 ** 5)]

    relevant_tifs_bounds = [(x, y, x + 10 ** 5, y + 10 ** 5) for x in x_list for y in y_list]

    return relevant_tifs_bounds


def identify_relevant_ESM_for_gdf(examined_gdf, ESM_bounds_csv, total_ESM_tif_list):

    # This function returns the list of the relevant ESM tif(s) for an examined gdf
    bounds_list = pd.read_csv(ESM_bounds_csv).set_index(['min_x', 'min_y', 'max_x', 'max_y'])

    index_list = []
    for i in calculate_bounds_of_relevant_tifs_for_gdf(examined_gdf):
        try:
            # Some related ESM data may not exist
            m = bounds_list.loc[i, 'index']
            index_list.append(m)
        except:
            continue
    index_list = list(set(index_list))
    ESM_tif_list = [total_ESM_tif_list[i] for i in index_list]

    return ESM_tif_list


def zonal_stats_tif_on_atom_gdf(atom_gdf, ESM_bounds_csv, total_ESM_tif_list, value_of_wanted_cells):

    # This function zonal_stats the raster on a vector
    # Suitable areas are added as a new column to the geodataframe

    ESM_tif_list = identify_relevant_ESM_for_gdf(atom_gdf, ESM_bounds_csv, total_ESM_tif_list)
    temp_list = []
    for ESM_tif in ESM_tif_list:

        suitable_area = rio.open(ESM_tif)
        cell_size_x = suitable_area.transform[0]
        cell_size_y = -suitable_area.transform[4]
        pixel_area = float(cell_size_x * cell_size_y / 10000)

        stat = zonal_stats(atom_gdf, ESM_tif, categorical=True)

        try:
            atom_gdf.loc[:, 'suitable_area_ha'] = sum([stat[0][i] * pixel_area for i in value_of_wanted_cells])
        except:
            atom_gdf.loc[:, 'suitable_area_ha'] = 0
        
        temp_list.append(atom_gdf[atom_gdf['suitable_area_ha'] != 0])

    atom_gdf_with_suitable_area = pd.concat(temp_list)

    return atom_gdf_with_suitable_area
    

def add_capacity_factors(atom_gdf_with_suitable_area, capacity_factor_csv):

    # This function adds capacity factors to municipalities
    capacity_factor = pd.read_csv(capacity_factor_csv)
    atom_gdf_with_capacity_factors = atom_gdf_with_suitable_area.merge(capacity_factor, how="inner", on='site_id')
    atom_gdf_with_capacity_factors = atom_gdf_with_capacity_factors.drop(['Unnamed: 0'], axis=1)

    return atom_gdf_with_capacity_factors


def summarize_suitable_area_on_atom_gdf(aggregated_gdf_with_capacity_factors):

    # This function summarizes the suitable areas and weighted average capacity factor of each municipality
    month_list = [str(i).zfill(2) for i in range(1, 13)]
    for i in month_list:
        aggregated_gdf_with_capacity_factors[i] = aggregated_gdf_with_capacity_factors[i] * aggregated_gdf_with_capacity_factors['suitable_area_ha']  

    aggregated_result = aggregated_gdf_with_capacity_factors.dissolve(by='municipal_id', aggfunc='sum')

    for i in month_list:
        aggregated_result[i] /= aggregated_result['suitable_area_ha']

    return aggregated_result


def synthesize(national_municipalities_geojson, municipal_geojson_with_site_id, regional_geojson, ESM_bounds_csv, total_ESM_tif_list, value_of_wanted_cells, capacity_factor_csv, output_geojson):

    # This function aggregates the above-mentioned functions
    # This purpose of creating a temporal index column is to guarantee the selected parts still remain to be geodataframe rather than geoseries
    
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
                        atom_gdf, ESM_bounds_csv, total_ESM_tif_list, value_of_wanted_cells
                        ), 
                    capacity_factor_csv
                    )
                )
        except:
            continue

    try:
        aggregated_gdf_with_capacity_factors = pd.concat(temp_list)
        aggregated_result = summarize_suitable_area_on_atom_gdf(aggregated_gdf_with_capacity_factors)
        aggregated_result = aggregated_result.drop('site_id', axis=1)
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
        ESM_bounds_csv = snakemake.input.ESM_bounds_csv, 
        total_ESM_tif_list = snakemake.input.ESM_list, 
        value_of_wanted_cells = snakemake.params.value_of_wanted_cells, 
        capacity_factor_csv = snakemake.input.synthesized_capacity_factor_csv, 
        output_geojson = snakemake.output[0]
        )