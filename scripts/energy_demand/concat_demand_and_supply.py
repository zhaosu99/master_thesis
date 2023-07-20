import geopandas as gpd
import pandas as pd


def concat_demand_and_supply(municipality_or_region, geoboundaries_geojson, demand_csv, supply_csv, output_csv):

    if municipality_or_region == 'municipality':
        merge_on = 'municipal_id'
    elif municipality_or_region == 'region':
        merge_on = 'regional_id'
    else:
        print('Input error')

    demand_df = pd.read_csv(demand_csv)
    supply_df = pd.read_csv(supply_csv)
    geoboundaries_gdf = gpd.read_file(geoboundaries_geojson)

    merged_df = supply_df.merge(demand_df, how='outer', on=merge_on).fillna(0)
    merged_gdf = geoboundaries_gdf.merge(merged_df, how='outer', on=merge_on).fillna(0)
    merged_gdf = merged_gdf.drop(['Unnamed: 0_x', 'Unnamed: 0_y'], axis=1)
    merged_gdf['index'] = merged_gdf[merge_on].apply(lambda x: float(x[9: ]))
    merged_gdf = merged_gdf.sort_values(by='index').set_index('index')
    
    result_df = merged_gdf.drop('geometry', axis=1)

    # Complete incomplete items
    def complete_all_items(input_df):
        wanted_cols = ['onshore_wind_MWh', 'offshore_wind_MWh', 'open_field_solar_PV_MWh',
            'rooftop_solar_PV_MWh', 'hydro_MWh', 'nuclear_MWh', 'biomass_MWh', 'geothermal_MWh', 
            'transport_electricity_MWh', 'industry_electricity_MWh', 'industry_heat_MWh', 'building_electricity_MWh']
        for col in wanted_cols:
            if col not in input_df.columns:
                input_df[col] = 0
        return input_df

    result_df = complete_all_items(result_df)
    
    # Calculate total energy supply
    supply_cols = ['onshore_wind_MWh', 'offshore_wind_MWh', 'open_field_solar_PV_MWh',
        'rooftop_solar_PV_MWh', 'hydro_MWh', 'nuclear_MWh', 'biomass_MWh', 'geothermal_MWh']
    result_df['total_supply_MWh'] = 0
    for i in supply_cols:
        result_df['total_supply_MWh'] += result_df[i]

    # Calculate total energy demand
    demand_cols = ['transport_electricity_MWh', 'industry_electricity_MWh', 'industry_heat_MWh', 'building_electricity_MWh']
    result_df['total_demand_MWh'] = 0
    for i in demand_cols:
        result_df['total_demand_MWh'] += result_df[i]
    
    # Calculate inventory
    result_df['inventory_MWh'] = result_df['total_supply_MWh'] - result_df['total_demand_MWh']

    result_df.to_csv(output_csv)
    

if __name__ =="__main__":
    concat_demand_and_supply(
        municipality_or_region = snakemake.params[0], 
        geoboundaries_geojson = snakemake.input[0], 
        demand_csv = snakemake.input[1], 
        supply_csv = snakemake.input[2], 
        output_csv = snakemake.output[0]
        )