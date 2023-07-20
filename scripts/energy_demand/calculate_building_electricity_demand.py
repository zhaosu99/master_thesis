import geopandas as gpd
import pandas as pd
from pathlib import Path
import pycountry


def sjoin_municipalities_to_when2heat(municipalities_geojson, when2heat_geojson):

    # This function appends the COP values to nearest municipalities 
    municipalities_gdf = gpd.read_file(municipalities_geojson)
    when2heat_gdf = gpd.read_file(when2heat_geojson)

    centroid_gdf = municipalities_gdf.copy().loc[:, ['municipal_id', 'population_ratio']]
    centroid_gdf['geometry'] = municipalities_gdf.geometry.centroid
    centroid_gdf = centroid_gdf.set_geometry('geometry').drop('population_ratio', axis=1)

    preprocessed_df = centroid_gdf.sjoin_nearest(when2heat_gdf, how='inner')
    preprocessed_df = preprocessed_df.drop(['index_right', 'geometry'], axis=1)

    result_gdf = municipalities_gdf.merge(preprocessed_df, on='municipal_id')

    return result_gdf


def get_code2(code_3):
    if code_3 == 'GRC':
        code_2 = 'EL'
    elif code_3 == "GBR":
        code_2 = 'UK'
    else:
        code_2 = pycountry.countries.get(alpha_3=code_3).alpha_2
    return code_2


def get_total_pop(country_codes, pop_csv):
    pop_df = pd.read_csv(pop_csv)
    total_pop = 0
    for i in country_codes:
        code_2 = get_code2(i)
        try:
            total_pop += float(pop_df.loc[lambda x: x.geo == code_2 and x.TIME_PERIOD == 2015, 'OBS_VALUE'])
        except:
            try:
                total_pop += pop_df.loc[lambda x: x.geo == code_2, 'OBS_VALUE'].mean()
            except:
                continue
    return total_pop


def get_energy_demand_of_each_country(total_pop, code_3, res_jrc_xlsx_list, com_jrc_xlsx_list, regression_csv, pop_csv):

    '''
    This functions calculates the national heat demands.

    Six results in sequence are: residential space heat demand, residential water heat demand, residential other energy demand, 
    commercial space heat demand, commercial water heat demand, commercial other energy demand
    '''

    code_2 = get_code2(code_3)

    regression_df = pd.read_csv(regression_csv)
    pop_df = pd.read_csv(pop_csv)
    jrc_codes2_list = [Path(i).stem[-2: ] for i in res_jrc_xlsx_list]
    m = 11630 # 1 ktoe = 11630 MWh

    # EU countries have direct values
    def get_jrc_valules(wanted_code2, jrc_codes2_list, res_jrc_xlsx_list, com_jrc_xlsx_list, m):

        index = jrc_codes2_list.index(wanted_code2)
        res_sum_df = pd.read_excel(open(res_jrc_xlsx_list[index], 'rb'), sheet_name='RES_summary')
        res_use_df = pd.read_excel(open(res_jrc_xlsx_list[index], 'rb'), sheet_name='RES_hh_tes')
        com_sum_df = pd.read_excel(open(com_jrc_xlsx_list[index], 'rb'), sheet_name='SER_summary')
        com_use_df = pd.read_excel(open(com_jrc_xlsx_list[index], 'rb'), sheet_name='SER_hh_tes')

        # To avoid double-counting, existing heat pump electricity demands are excluded
        res_space_heat_demand = float(res_use_df.loc[2, 2015]) * m
        res_water_heat_demand = float(res_use_df.loc[15, 2015]) * m
        res_other_energy_demand = (float(res_use_df.loc[13, 2015]) + float(res_use_df.loc[25, 2015]) + float(res_sum_df.loc[55, 2015])) * m

        com_space_heat_demand = float(com_use_df.loc[2, 2015]) * m
        com_water_heat_demand = float(com_use_df.loc[17, 2015]) * m
        com_other_energy_demand = (float(com_use_df.loc[14, 2015]) + float(com_use_df.loc[27, 2015]) + float(com_sum_df.loc[58, 2015])) * m

        return res_space_heat_demand, res_water_heat_demand, res_other_energy_demand, com_space_heat_demand, com_water_heat_demand, com_other_energy_demand

    # Values of non-EU countries are infered from population and regression function
    def get_regression_values(total_pop, wanted_code2, regression_df, pop_df, m):

        name_list = ['res_space_heat', 'res_water_heat', 'res_other_energy', 'com_space_heat', 'com_water_heat', 'com_other_energy']
        try:
            pop = float(pop_df.loc[lambda x: x.geo == wanted_code2 and x.TIME_PERIOD == 2015, 'OBS_VALUE'])
        except:
            pop = pop_df.loc[lambda x: x.geo == wanted_code2, 'OBS_VALUE'].mean()

        value_list = []
        for i in name_list:
            slope = float(regression_df.loc[lambda x: x.demand_type == i, 'slope'])
            intercept = float(regression_df.loc[lambda x: x.demand_type == i, 'intercept'])
            value_list.append((slope * total_pop + intercept) * m * pop / total_pop)

        return tuple(value_list)

    # Synthesize two functions
    if code_2 in jrc_codes2_list:
        return get_jrc_valules(code_2, jrc_codes2_list, res_jrc_xlsx_list, com_jrc_xlsx_list, m)
    else:
        return get_regression_values(total_pop, code_2, regression_df, pop_df, m)

        
def calculate_building_energy_demand(country_codes, municipalities_geojson, when2heat_geojson, res_jrc_xlsx_list, com_jrc_xlsx_list, 
    regression_csv, pop_csv, output_geojson):

    # This function calculates six types of electricity demand and their sum in the studied country
    code_3 = Path(municipalities_geojson).stem[0: 3]
    total_pop = get_total_pop(country_codes, pop_csv)

    result_gdf = sjoin_municipalities_to_when2heat(municipalities_geojson, when2heat_geojson)

    res_space_heat_demand, res_water_heat_demand, res_other_energy_demand, com_space_heat_demand, com_water_heat_demand, com_other_energy_demand = \
        get_energy_demand_of_each_country(total_pop, code_3, res_jrc_xlsx_list, com_jrc_xlsx_list, regression_csv, pop_csv)
    
    # Residential space heat demand
    result_gdf['res_space_heat_electricity_demand'] = res_space_heat_demand * result_gdf['population_ratio'] * result_gdf['residental_radiator'] 

    # Commercial space heat demand
    result_gdf['com_space_heat_electricity_demand'] = com_space_heat_demand * result_gdf['population_ratio'] * result_gdf['commercial_radiator'] 

    # Residential water heat demand
    result_gdf['res_water_heat_electricity_demand'] = res_water_heat_demand * result_gdf['population_ratio'] * result_gdf['residental_water']
    
    # Commercial water heat demand
    result_gdf['com_water_heat_electricity_demand'] = com_water_heat_demand * result_gdf['population_ratio'] * result_gdf['commercial_water']
    
    # Residential other demands
    result_gdf['res_other_electricity_demand'] = result_gdf['population_ratio'] * res_other_energy_demand

    # Commercial other demands
    result_gdf['com_other_electricity_demand'] = result_gdf['population_ratio'] * com_other_energy_demand

    # Total electricity demand
    result_gdf['total_building_electricity_demand'] = \
        result_gdf['res_space_heat_electricity_demand'] + result_gdf['res_water_heat_electricity_demand'] + result_gdf['res_other_electricity_demand'] + \
        result_gdf['com_space_heat_electricity_demand'] + result_gdf['com_water_heat_electricity_demand'] + result_gdf['com_other_electricity_demand']

    kept_cols = ['geometry', 'municipal_id', 'res_space_heat_electricity_demand', 'res_water_heat_electricity_demand', 'res_other_electricity_demand',\
        'com_space_heat_electricity_demand', 'com_water_heat_electricity_demand', 'com_other_electricity_demand', 'total_building_electricity_demand']
    result_gdf = result_gdf.loc[:, kept_cols]

    result_gdf.to_file(output_geojson)


if __name__ == "__main__":
    calculate_building_energy_demand(
        country_codes = snakemake.params[0],
        municipalities_geojson = snakemake.input.pop_ratio_geojson, 
        when2heat_geojson = snakemake.input.heat_electricity_ratio_geojson, 
        res_jrc_xlsx_list = snakemake.input.res_jrc_xlsx_list, 
        com_jrc_xlsx_list = snakemake.input.com_jrc_xlsx_list, 
        regression_csv = snakemake.input.regression_csv, 
        pop_csv = snakemake.input.pop_csv, 
        output_geojson = snakemake.output[0]
        )