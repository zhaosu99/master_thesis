import geopandas as gpd
import pandas as pd
import pycountry
from pathlib import Path


def calculate_industrial_energy_demand(country_codes_list, input_geojson_list, input_ETS_xlsx, jrc_xlsx_list, input_backup_csv, output_geojson_list):
 
    # Calculate average electricity / total energy consumption ratio
    e_sum, sum = 0, 0
    for i in jrc_xlsx_list:
        summary_sheet = pd.read_excel(open(i, 'rb'), sheet_name='Ind_Summary')
        e_sum += float(summary_sheet.loc[47, 2015])
        sum += float(summary_sheet.loc[48, 2015]) 
    electricity_ratio = e_sum / sum

    # Get the two-letter country codes of studied countries
    all_codes2_list = []
    for i in country_codes_list:
        if i == 'GRC':
            all_codes2_list.append('EL')
        elif i == "GBR":
            all_codes2_list.append('UK')
        else:
            all_codes2_list.append(pycountry.countries.get(alpha_3=i).alpha_2)

    # Preprocess raw ETS data to geodataframe
    input_df = pd.read_excel(input_ETS_xlsx).loc[:, ['Country', 'Latitude', 'Longitude', 'Subsector', 'Emissions_ETS_2014']]
    input_df = input_df.dropna(subset=['Country'])
    country_name_list = list(input_df.Country)
    country_code_2 = [pycountry.countries.search_fuzzy(i)[0].alpha_2 for i in country_name_list]
    input_df['Country'] = country_code_2 
    ETS_gdf = gpd.GeoDataFrame(input_df, geometry=gpd.points_from_xy(input_df.Latitude, input_df.Longitude), crs=4326).to_crs(3035)

    # These countries have detailed data about their energy consumption in industry sector
    existing_codes2_list = [Path(i).stem[-2: ] for i in jrc_xlsx_list]

    m = 11630 # 1 ktoe = 11630 MWh 

    for n in range(len(country_codes_list)):

        code_2 = all_codes2_list[n]
        
        input_gdf = gpd.read_file(input_geojson_list[n])

        # Get the population ratio in each municipalities
        pop_sum = input_gdf.population_amount.sum()
        pop_gdf = input_gdf.copy()
        pop_gdf['population_ratio'] = pop_gdf['population_amount'] / pop_sum
        pop_df = pop_gdf.copy().loc[:, ['municipal_id', 'population_ratio']]
        
        if code_2 in existing_codes2_list:
            index = existing_codes2_list.index(code_2)
                 
            ETS_gdf_code2 = ETS_gdf.copy()
            ETS_gdf_code2 = ETS_gdf_code2[ETS_gdf_code2.Country == code_2]
            
            raw_gdf = input_gdf.copy().sjoin(ETS_gdf_code2, how='inner', predicate='intersects').dropna()
            raw_gdf = raw_gdf.loc[:, ['municipal_id', 'Subsector', 'Emissions_ETS_2014']]
            raw_gdf = raw_gdf.rename(columns={'Emissions_ETS_2014': 'emission'})

            # Get the emission in relevent municipaltis by sector
            iron = raw_gdf[raw_gdf['Subsector']=='Iron and steel'].loc[:, ['municipal_id', 'emission']]
            iron = iron[iron['emission'] != 0]
            sum = iron.emission.sum()
            iron['emission_ratio'] = iron['emission'] / sum
            iron = iron.loc[:, ['municipal_id', 'emission_ratio']]
            
            non_ferrous = raw_gdf[raw_gdf['Subsector']=='Non-ferrous metals'].loc[:, ['municipal_id', 'emission']]
            non_ferrous = non_ferrous[non_ferrous['emission'] != 0]
            sum = non_ferrous.emission.sum()
            non_ferrous['emission_ratio'] = non_ferrous['emission'] / sum
            non_ferrous = non_ferrous.loc[:, ['municipal_id', 'emission_ratio']]

            chemicals = raw_gdf[raw_gdf['Subsector']=='Chemical industry'].loc[:, ['municipal_id', 'emission']]
            chemicals = chemicals[chemicals['emission'] != 0]
            sum = chemicals.emission.sum()
            chemicals['emission_ratio'] = chemicals['emission'] / sum
            chemicals = chemicals.loc[:, ['municipal_id', 'emission_ratio']]
            
            non_metallic = raw_gdf[raw_gdf['Subsector'].isin(['Non-ferrous metals', 'Cement', 'Glass'])].loc[:, ['municipal_id', 'emission']]
            non_metallic = non_metallic[non_metallic['emission'] != 0]
            sum = non_metallic.emission.sum()
            non_metallic['emission_ratio'] = non_metallic['emission'] / sum
            non_metallic = non_metallic.loc[:, ['municipal_id', 'emission_ratio']]
            
            paper = raw_gdf[raw_gdf['Subsector']=='Paper and printing'].loc[:, ['municipal_id', 'emission']]
            paper = paper[paper['emission'] != 0]
            sum = paper.emission.sum()
            paper['emission_ratio'] = paper['emission'] / sum
            paper = paper.loc[:, ['municipal_id', 'emission_ratio']]
            
            # Get national heat and electricity demand in iron and steel sector
            isl = pd.read_excel(open(jrc_xlsx_list[index], 'rb'), sheet_name='ISI_ued')
            isl_electricity = isl.loc[4: 8, 2015].sum()
            isl_heat = float(isl.loc[3, 2015]) - isl_electricity
            iron['iron_electricity'] = iron['emission_ratio'] * isl_electricity * m
            iron['iron_heat'] = iron['emission_ratio'] * isl_heat * m
            iron = iron.loc[:, ['municipal_id', 'iron_electricity', 'iron_heat']]

            # Get national heat and electricity demand in non-ferrous metal sector
            nfm = pd.read_excel(open(jrc_xlsx_list[index], 'rb'), sheet_name='NFM_ued')
            nfm_electricity = nfm.loc[4: 8, 2015].sum() + nfm.loc[32: 36, 2015].sum() + nfm.loc[41, 2015].sum() + nfm.loc[48, 2015].sum() + \
                nfm.loc[69: 73, 2015].sum() + nfm.loc[84, 2015].sum() + nfm.loc[91, 2015].sum() + nfm.loc[108, 2015].sum() + \
                nfm.loc[111: 115, 2015].sum() + nfm.loc[127, 2015].sum() + nfm.loc[134, 2015].sum() + nfm.loc[151, 2015].sum()
            nfm_heat = - nfm_electricity + float(nfm.loc[3, 2015]) + float(nfm.loc[31, 2015]) + float(nfm.loc[68, 2015]) + float(nfm.loc[110, 2015])  
            non_ferrous['non_ferrous_electricity'] = non_ferrous['emission_ratio'] * nfm_electricity * m
            non_ferrous['non_ferrous_heat'] = non_ferrous['emission_ratio'] * nfm_heat * m
            non_ferrous = non_ferrous.loc[:, ['municipal_id', 'non_ferrous_electricity', 'non_ferrous_heat']]

            # Get national heat and electricity demand in chemical industry sector
            chi = pd.read_excel(open(jrc_xlsx_list[index], 'rb'), sheet_name='CHI_ued')
            chi_electricity = chi.loc[4: 8, 2015].sum() + chi.loc[40, 2015].sum() + chi.loc[54, 2015].sum() + \
                chi.loc[59: 63, 2015].sum() + chi.loc[80, 2015].sum() + chi.loc[88, 2015].sum() + \
                chi.loc[102: 103, 2015].sum() + chi.loc[107: 111, 2015].sum() + chi.loc[128, 2015].sum() + \
                chi.loc[136, 2015].sum() + chi.loc[150: 151, 2015].sum()
            chi_heat = - chi_electricity + float(chi.loc[3, 2015]) + float(chi.loc[58, 2015]) + float(chi.loc[106, 2015]) - float(chi.loc[13, 2015])
            chemicals['chemicals_electricity'] = chemicals['emission_ratio'] * chi_electricity * m
            chemicals['chemicals_heat'] = chemicals['emission_ratio'] * chi_heat * m
            chemicals = chemicals.loc[:, ['municipal_id', 'chemicals_electricity', 'chemicals_heat']]

            # Get national heat and electricity demand in non-metallic industry sector
            nmm = pd.read_excel(open(jrc_xlsx_list[index], 'rb'), sheet_name='NMM_ued')
            nmm_electricity = nmm.loc[4: 8, 2015].sum() + nmm.loc[46: 50, 2015].sum() + nmm.loc[96: 100, 2015].sum()
            nmm_heat = - nmm_electricity + float(nmm.loc[3, 2015]) + float(nmm.loc[45, 2015]) + float(nmm.loc[95, 2015]) 
            non_metallic['non_metallic_electricity'] = non_metallic['emission_ratio'] * nmm_electricity * m
            non_metallic['non_metallic_heat'] = non_metallic['emission_ratio'] * nmm_heat * m
            non_metallic = non_metallic.loc[:, ['municipal_id', 'non_metallic_electricity', 'non_metallic_heat']]

            # Get national heat and electricity demand in pulp and paper sector
            ppa = pd.read_excel(open(jrc_xlsx_list[index], 'rb'), sheet_name='PPA_ued')
            ppa_electricity = ppa.loc[4: 8, 2015].sum() + ppa.loc[30: 34, 2015].sum() + ppa.loc[51, 2015].sum() + \
                ppa.loc[64, 2015].sum() + ppa.loc[77, 2015].sum() + ppa.loc[79, 2015].sum()
            ppa_heat = - ppa_electricity + float(ppa.loc[3, 2015]) + float(ppa.loc[29, 2015])
            paper['paper_electricity'] = paper['emission_ratio'] * ppa_electricity * m
            paper['paper_heat'] = paper['emission_ratio'] * ppa_heat * m
            paper = paper.loc[:, ['municipal_id', 'paper_electricity', 'paper_heat']]

            # Get national heat and electricity demand in other sectors
            summary_sheet = pd.read_excel(open(jrc_xlsx_list[index], 'rb'), sheet_name='Ind_Summary')
            total_electricity = float(summary_sheet.loc[47, 2015])
            total_heat = float(summary_sheet.loc[48, 2015]) - total_electricity 
            other_electricity = total_electricity - isl_electricity - nfm_electricity - chi_electricity - nmm_electricity - ppa_electricity
            other_heat = total_heat - isl_heat - nfm_heat - chi_heat - nmm_heat - ppa_heat

            pop_df_copy = pop_df.copy()
            pop_df_copy['other_electricity'] = pop_df_copy['population_ratio'] * other_electricity * m
            pop_df_copy['other_heat'] = pop_df_copy['population_ratio'] * other_heat * m
            pop_df_copy = pop_df_copy.loc[:, ['municipal_id', 'other_electricity', 'other_heat']]

            result_df = pd.concat([iron, non_ferrous, chemicals, non_metallic, paper, pop_df_copy], ignore_index=True, sort=False)
            result_df = result_df.fillna(0)
            result_df = result_df.groupby('municipal_id').agg('sum')
            result_df = result_df.reset_index(drop=False)
            result_gdf = input_gdf.merge(result_df, on='municipal_id')
            result_gdf = result_gdf.drop('population_amount', axis=1)
            
        else:
            backup_df = pd.read_csv(input_backup_csv)
            backup_df_copy = backup_df[backup_df.nrg_bal=='FC_IND_E'].copy()
            if code_2 != 'CH':
                backup_df_copy = backup_df_copy[backup_df_copy.geo==code_2]
                industrial_energy = float(backup_df_copy.loc[backup_df_copy.TIME_PERIOD==2019, 'OBS_VALUE'])
            else:
                industrial_energy = 2870 # Unit: ktoe https://www.odyssee-mure.eu/publications/efficiency-trends-policies-profiles/switzerland.html
            
            result_gdf = pop_gdf.copy()
            result_gdf['iron_electricity'] = 0
            result_gdf['iron_heat'] = 0
            result_gdf['non_ferrous_electricity'] = 0
            result_gdf['non_ferrous_heat'] = 0
            result_gdf['chemicals_electricity'] = 0
            result_gdf['chemicals_heat'] = 0
            result_gdf['non_metallic_electricity'] = 0
            result_gdf['non_metallic_heat'] = 0
            result_gdf['paper_electricity'] = 0
            result_gdf['paper_heat'] = 0
            result_gdf['other_electricity'] = result_gdf['population_ratio'] * electricity_ratio * m * industrial_energy
            result_gdf['other_heat'] = result_gdf['population_ratio'] * (1 - electricity_ratio) * m * industrial_energy

            # Drop useless columns
            try:
                result_gdf = result_gdf.drop(['population_amount'], axis=1)
            except:
                continue
            try:
                result_gdf = result_gdf.drop(['population_ratio'], axis=1)
            except:
                continue

        # Calculate total electricity and heat demand in industry sector
        electricity_list = ['iron_electricity', 'non_ferrous_electricity', 'chemicals_electricity',
            'non_metallic_electricity', 'paper_electricity', 'other_electricity']
        result_gdf['total_industry_electricity_demand'] = 0
        for i in electricity_list:
            result_gdf['total_industry_electricity_demand'] += result_gdf[i]

        heat_list = ['iron_electricity', 'iron_heat', 'non_ferrous_heat', 'chemicals_heat', 'non_metallic_heat',
            'paper_heat', 'other_heat']
        result_gdf['total_industry_heat_demand'] = 0
        for i in heat_list:
            result_gdf['total_industry_heat_demand'] += result_gdf[i]
            
        result_gdf.to_file(output_geojson_list[n])


if __name__ == "__main__":
    calculate_industrial_energy_demand(
        country_codes_list = snakemake.params.country_codes, 
        input_geojson_list = snakemake.input.municipalities, 
        input_ETS_xlsx = snakemake.input.ETS_data, 
        jrc_xlsx_list = snakemake.input.jrc_xlsx_list, 
        input_backup_csv = snakemake.input.backup_csv, 
        output_geojson_list = snakemake.output
        )