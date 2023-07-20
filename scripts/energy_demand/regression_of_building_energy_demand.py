import pandas as pd
from pathlib import Path
import scipy.stats as st


def regression_of_building_energy(pop_csv, res_jrc_xlsx_list, com_jrc_xlsx_list, output_csv):

    # This function calculates the regression function between population and different types of energy demand
    pop_df = pd.read_csv(pop_csv)
    
    res_jrc_xlsx_list = [i for i in res_jrc_xlsx_list if Path(i).stem[-2: ] != 'UK']
    com_jrc_xlsx_list = [i for i in com_jrc_xlsx_list if Path(i).stem[-2: ] != 'UK']
    codes2_list = [Path(i).stem[-2: ] for i in res_jrc_xlsx_list]
    
    # get population of exising JRC countries
    pop_list = []
    for i in codes2_list:
        wanted_pop_df = pop_df[pop_df['geo'] == i].copy()
        try:
            wanted_value = float(wanted_pop_df.loc[lambda x: x.TIME_PERIOD == 2015, 'OBS_VALUE'])
        except:
            wanted_value = sum(list(wanted_pop_df.OBS_VALUE)) / len(list(wanted_pop_df.OBS_VALUE))
        pop_list.append(wanted_value)
    
    # get space and water heat demand for residential and tertiary/(commercial) buildings
    res_space_heat_list, res_water_heat_list, res_other_energy_list = [], [], []
    com_space_heat_list, com_water_heat_list, com_other_energy_list = [], [], []
    
    length = len(res_jrc_xlsx_list)
    for index in range(length):

        res_sum_df = pd.read_excel(open(res_jrc_xlsx_list[index], 'rb'), sheet_name='RES_summary')
        res_use_df = pd.read_excel(open(res_jrc_xlsx_list[index], 'rb'), sheet_name='RES_hh_tes')
        com_sum_df = pd.read_excel(open(com_jrc_xlsx_list[index], 'rb'), sheet_name='SER_summary')
        com_use_df = pd.read_excel(open(com_jrc_xlsx_list[index], 'rb'), sheet_name='SER_hh_tes')

        res_space_heat_list.append(float(res_use_df.loc[2, 2015]))
        res_water_heat_list.append(float(res_use_df.loc[15, 2015]))
        res_other = float(res_use_df.loc[13, 2015]) + float(res_use_df.loc[25, 2015]) + float(res_sum_df.loc[55, 2015])
        res_other_energy_list.append(res_other)

        com_space_heat_list.append(float(com_use_df.loc[2, 2015]))
        com_water_heat_list.append(float(com_use_df.loc[17, 2015]))
        com_other = float(com_use_df.loc[14, 2015]) + float(com_use_df.loc[27, 2015]) + float(com_sum_df.loc[58, 2015])
        com_other_energy_list.append(com_other)

    y_value_list = [res_space_heat_list, res_water_heat_list, res_other_energy_list, com_space_heat_list, com_water_heat_list, com_other_energy_list]

    slope_list, intercept_list, r2_value_list = [], [], []
    for n in range(len(y_value_list)):
        locals()[f'slope_{n}'], locals()[f'intercept_{n}'], locals()[f'r_value_{n}'], locals()[f'p_value_{n}'], locals()[f'std_err_{n}'] = st.linregress(pop_list, y_value_list[n])
        slope_list.append(locals()[f'slope_{n}'])
        intercept_list.append(locals()[f'intercept_{n}'])
        r2_value_list.append(locals()[f'r_value_{n}'] ** 2)

    output_df = pd.DataFrame({
        'demand_type': ['res_space_heat', 'res_water_heat', 'res_other_energy', 'com_space_heat', 'com_water_heat', 'com_other_energy'],
        'slope': slope_list,
        'intercept': intercept_list,
        'r2': r2_value_list
        })

    output_df.to_csv(output_csv)


if __name__ == '__main__':
    regression_of_building_energy(
        pop_csv = snakemake.input.pop_csv, 
        res_jrc_xlsx_list = snakemake.input.res_jrc_xlsx_list, 
        com_jrc_xlsx_list = snakemake.input.com_jrc_xlsx_list, 
        output_csv = snakemake.output[0]
        )