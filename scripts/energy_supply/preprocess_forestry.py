import pandas as pd
import pycountry


def get_annual_average_production(input_csv, country_list, column_name):

    # Calculate annual average production for a certain forestry product type
    production_list = []
    input_df = pd.read_csv(input_csv)
    for i in country_list:
        try:
            if i == 'GRC':
                alpha_2 = 'EL'
            elif i == "GBR":
                alpha_2 = 'UK'
            else:
                alpha_2 = pycountry.countries.get(alpha_3=i).alpha_2
            temp_df = input_df.loc[lambda input_df: input_df.geo==alpha_2].copy()
            production_list.append(temp_df.OBS_VALUE.mean())
        except:
            production_list.append(0)

    annual_average_df = pd.DataFrame({
        'country_codes': country_list,
        column_name: [1000 * i for i in production_list]
        })
    
    return annual_average_df


def synthesize(input_csv_f, input_csv_r, country_list, output_csv):

    # This function summarizes the national annual production of fuelwood and roundwood
    df_f = get_annual_average_production(input_csv_f, country_list, 'fuelwood_m3')
    df_r = get_annual_average_production(input_csv_r, country_list, 'roundwood_m3')

    synthesize_df = pd.DataFrame({
        'country_codes': country_list,
        'fuelwood_m3': df_f.fuelwood_m3,
        'roundwood_m3': df_r.roundwood_m3,
        })

    synthesize_df = synthesize_df.fillna(0)
    synthesize_df.to_csv(output_csv)


if __name__ == "__main__":
    synthesize(
        input_csv_f = snakemake.input.fuelwood, 
        input_csv_r = snakemake.input.roundwood, 
        country_list = snakemake.params[0], 
        output_csv = snakemake.output[0]
        )