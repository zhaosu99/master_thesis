import geopandas as gpd
import pandas as pd
from pathlib import Path
import pycountry
import numpy as np


def calculate_other_transport_electricity_demand(pop_amount_geojson, regression_csv, jrc_xlsx_list, output_geojson):
    
    # Energy conversion efficiency of different types of engines
    efficiency_diesel_oil_engine = 0.5
    efficiency_electric_engine = 0.85
    efficiency_hydrogen_engine = 0.4

    car_ownership_rate = 0.7 # Unit number of vehicles/population
    annual_mileage = 11000 # Unit km/a
    electricity_consumption = 0.2 # Unit kWh/km

    # Get the two-letter country codes of studied countries
    code_3 = Path(pop_amount_geojson).stem[:3]
    if code_3 == 'GRC':
        code_2 = 'EL'
    elif code_3 == "GBR":
        code_2 = 'UK'
    else:
        code_2 = pycountry.countries.get(alpha_3=code_3).alpha_2

    # These countries have detailed data about their energy consumption in transport sector
    existing_codes2_list = [Path(i).stem[-2: ] for i in jrc_xlsx_list]
    
    m = 11630 # 1 ktoe = 11630 MWh

    if code_2 in existing_codes2_list:

        index = existing_codes2_list.index(code_2)
        road_df = pd.read_excel(open(jrc_xlsx_list[index], 'rb'), sheet_name='TrRoad_ene')
        railway_list_df = pd.read_excel(open(jrc_xlsx_list[index], 'rb'), sheet_name='TrRail_ene')

        # Electricity demand of public passenger transport
        pp_sum = float(road_df.loc[31, 2015]) / efficiency_diesel_oil_engine * efficiency_electric_engine 

        # Electricity demand of light freight transport
        lf_sum = float(road_df.loc[41, 2015]) / efficiency_diesel_oil_engine * efficiency_electric_engine 

        # Electricity demand of heavy freight transport
        hf_sum = float(road_df.loc[44, 2015]) / efficiency_diesel_oil_engine * efficiency_hydrogen_engine

        # Electricity demand of rail tranport
        r_diesel_oil_engine = (float(railway_list_df.loc[18, 2015]) + float(railway_list_df.loc[22, 2015])) / efficiency_diesel_oil_engine
        r_electric_engine = float(railway_list_df.loc[19, 2015]) + float(railway_list_df.loc[23, 2015]) 
        r_metro = float(railway_list_df.loc[16, 2015])
        r_sum = (r_diesel_oil_engine + r_electric_engine + r_metro) * efficiency_electric_engine

        # Regionalize national data to municipal level
        pop_amount_gdf = gpd.read_file(pop_amount_geojson)
        pop_amount_gdf['population_ratio'] = pop_amount_gdf['population_amount'] / pop_amount_gdf.population_amount.sum()
        pop_amount_gdf['private_passenger_transport_MWh'] = pop_amount_gdf['population_amount'] * car_ownership_rate * annual_mileage * electricity_consumption / 1000
        pop_amount_gdf['public_passenger_transport_MWh'] = pop_amount_gdf['population_ratio'] * pp_sum * m
        pop_amount_gdf['light_freight_transport_MWh'] = pop_amount_gdf['population_ratio'] * lf_sum * m
        pop_amount_gdf['heavy_freight_transport_MWh'] = pop_amount_gdf['population_ratio'] * hf_sum * m
        pop_amount_gdf['railway_transport_MWh'] = pop_amount_gdf['population_ratio'] * r_sum * m
        
        pop_amount_gdf = pop_amount_gdf.drop(columns=['population_ratio', 'population_amount'], axis=1)
    
    else:
        pop_amount_gdf = gpd.read_file(pop_amount_geojson)    
        regression_df = pd.read_csv(regression_csv)

        pop_sum = pop_amount_gdf.population_amount.sum()
        pop_amount_gdf['population_ratio'] = pop_amount_gdf['population_amount'] / pop_sum
        
        k_list, b_list = list(regression_df.loc[:, 'slope']), list(regression_df.loc[:, 'intercept'])

        # Electricity demand of private passenger transport
        pop_amount_gdf['private_passenger_transport_MWh'] = pop_amount_gdf['population_amount'] * car_ownership_rate * annual_mileage * electricity_consumption / 1000

        # Electricity demand of public passenger transport
        national_pp = (k_list[0] * pop_sum + b_list[0]) / efficiency_diesel_oil_engine * efficiency_electric_engine      
        pop_amount_gdf['public_passenger_transport_MWh'] = pop_amount_gdf['population_ratio'] * national_pp * m

        # Electricity demand of light freight transport
        national_lf = (k_list[1] * pop_sum + b_list[1]) / efficiency_diesel_oil_engine * efficiency_electric_engine
        pop_amount_gdf['light_freight_transport_MWh'] = pop_amount_gdf['population_ratio'] * national_lf * m

        # Electricity demand of heavy freight transport
        national_hf = (k_list[2] * pop_sum + b_list[2]) / efficiency_diesel_oil_engine * efficiency_hydrogen_engine
        pop_amount_gdf['heavy_freight_transport_MWh'] = pop_amount_gdf['population_ratio'] * national_hf * m

        # Electricity demand of rail tranport
        r_efficiency_list = [efficiency_diesel_oil_engine, efficiency_electric_engine, efficiency_electric_engine]
        national_r = ((np.array(k_list[3: ]) * pop_sum + np.array(b_list[3: ])) / np.array(r_efficiency_list) * efficiency_electric_engine).sum()
        pop_amount_gdf['railway_transport_MWh'] = pop_amount_gdf['population_ratio'] * national_r * m
        
        pop_amount_gdf = pop_amount_gdf.drop(columns=['population_ratio', 'population_amount'], axis=1)

    value_list = ['private_passenger_transport_MWh', 'public_passenger_transport_MWh', 'light_freight_transport_MWh',
        'heavy_freight_transport_MWh', 'railway_transport_MWh']
    pop_amount_gdf['total_transport_electricity_demand'] = 0
    for i in value_list:
        pop_amount_gdf['total_transport_electricity_demand'] += pop_amount_gdf[i]

    pop_amount_gdf.to_file(output_geojson)


if __name__ == "__main__":
    calculate_other_transport_electricity_demand(
        pop_amount_geojson = snakemake.input.pop_amount_geojson, 
        regression_csv = snakemake.input.regression_csv, 
        jrc_xlsx_list = snakemake.input.jrc_xlsx_list, 
        output_geojson = snakemake.output[0]
        )