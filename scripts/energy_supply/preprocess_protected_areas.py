import geopandas as gpd
import pandas as pd


def select_protected_areas(input_geojson_list, input_eez_list, country_list, outoput_geojson_list):

    # This function selects qualified marine protected areas

    # The filter is in accordance to the way UNEP-WCMC calculates statistics:
    # https://www.protectedplanet.net/c/calculating-protected-area-coverage
    proteced_areas_list= []
    for i in input_geojson_list:
        proteced_areas = gpd.read_file(i)
        proteced_areas = proteced_areas[proteced_areas['MARINE'].isin(['1', '2'])]
        proteced_areas = proteced_areas[proteced_areas['ISO3'].isin(country_list)]
        proteced_areas = proteced_areas[proteced_areas['STATUS'].isin(['Designated', 'Inscribed', 'Established'])]
        proteced_areas = proteced_areas[proteced_areas['DESIG_ENG'] != 'UNESCO-MAB Biosphere Reserve']
        proteced_areas = proteced_areas.loc[:, ['geometry', 'ISO3']]
        proteced_areas_list.append(proteced_areas)
    merged_geojson = pd.concat(proteced_areas_list)

    # Pick out marine proteced areas for each country
    length = len(input_eez_list)
    for n in range(length):
        selected_geojson = merged_geojson[merged_geojson['ISO3']==country_list[n]]
        selected_geojson = selected_geojson.to_crs(3035)
        selected_geojson.to_file(outoput_geojson_list[n])


if __name__ == "__main__":
    select_protected_areas(
        input_geojson_list = snakemake.input.protected_areas, 
        input_eez_list = snakemake.input.eez_by_country,
        country_list = snakemake.params[0], 
        outoput_geojson_list = snakemake.output.offshore_protected_areas_by_country
        )