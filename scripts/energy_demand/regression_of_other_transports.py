import pandas as pd
import pycountry
from pathlib import Path
import scipy.stats as st


def regression_of_other_transports(pop_csv, jrc_xlsx_list, output_csv):

    # This function calculates the regression function between population and energy consumption of differnt types of engines
    pop_df = pd.read_csv(pop_csv)
    codes2_list = [Path(i).stem[-2: ] for i in jrc_xlsx_list]

    # Get 3-letter country codes of exising JRC countries
    codes3_list = []
    for i in codes2_list:
        if i == 'EL':
            i = 'GR'
        elif i == 'UK':
            i = 'GB'
        codes3_list.append(pycountry.countries.get(alpha_2=i).alpha_3)
    
    # get population of exising JRC countries
    pop_list = []
    for i in codes2_list:
        wanted_pop_df = pop_df[pop_df['geo'] == i].copy()
        try:
            wanted_value = float(wanted_pop_df.loc[lambda x: x.TIME_PERIOD == 2015, 'OBS_VALUE'])
        except:
            wanted_value = sum(list(wanted_pop_df.OBS_VALUE)) / len(list(wanted_pop_df.OBS_VALUE))
        pop_list.append(wanted_value)
    
    # Get regression function for public passenger road transport
    pp_engine_list = []
    lf_engine_list = []
    hf_engine_list = []
    r_diesel_oil_engine_list, r_electric_engine_list, r_metro_list = [], [], []
    
    for i in jrc_xlsx_list:

        road_df = pd.read_excel(open(i, 'rb'), sheet_name='TrRoad_ene')
        railway_list_df = pd.read_excel(open(i, 'rb'), sheet_name='TrRail_ene')

        pp_engine_list.append(float(road_df.loc[31, 2015]))

        lf_engine_list.append(float(road_df.loc[41, 2015]))

        hf_engine_list.append(float(road_df.loc[44, 2015])) 

        r_diesel_oil_engine_list.append((float(railway_list_df.loc[18, 2015]) + float(railway_list_df.loc[22, 2015])))
        r_electric_engine_list.append(float(railway_list_df.loc[19, 2015]) + float(railway_list_df.loc[23, 2015]))
        r_metro_list.append(float(railway_list_df.loc[16, 2015]))

    value_list = [pp_engine_list, lf_engine_list, hf_engine_list, r_diesel_oil_engine_list, r_electric_engine_list, r_metro_list]

    slope_list, intercept_list, r2_value_list = [], [], []
    for n in range(len(value_list)):
        locals()[f'slope_{n}'], locals()[f'intercept_{n}'], locals()[f'r_value_{n}'], locals()[f'p_value_{n}'], locals()[f'std_err_{n}'] = st.linregress(pop_list, value_list[n])
        slope_list.append(locals()[f'slope_{n}'])
        intercept_list.append(locals()[f'intercept_{n}'])
        r2_value_list.append(locals()[f'r_value_{n}']**2)

    output_df = pd.DataFrame({
        'transport_type': ['pp', 'lf', 'hf', 'r_diesel_oil', 'r_electric', 'r_metro'],
        'slope': slope_list,
        'intercept': intercept_list,
        'r2': r2_value_list
    })

    output_df.to_csv(output_csv)


if __name__ == '__main__':
    regression_of_other_transports(
        pop_csv = snakemake.input.pop_csv, 
        jrc_xlsx_list = snakemake.input.jrc_xlsx_list, 
        output_csv = snakemake.output[0]
        )