import pandas as pd
import pycountry


def get_annual_average_production(input_csv, crop_type, country_list):

    # Calculate annual average production for a certain crop type
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
        'crop_type': crop_type,
        'annual_production_t': [1000 * i for i in production_list]
        })
    
    return annual_average_df


def get_annual_average_production_rice(input_csv, country_list):

    # Calculate annual average production for rice
    # Here a different csv data is used
    production_list = []
    input_df = pd.read_csv(input_csv)
    for i in country_list:
        try:
            temp_df = input_df.loc[lambda input_df: input_df.Code==i].copy()
            temp_df = temp_df.loc[lambda temp_df: temp_df.Year>=2011]
            production_list.append(temp_df['Rice | 00000027 || Production | 005510 || tonnes'].mean())
        except:
            production_list.append(0)

    annual_average_df = pd.DataFrame({
        'country_codes': country_list,
        'crop_type': 'rice',
        'annual_production_t': production_list
        })
    
    return annual_average_df


def synthesize(wheat_csv, rye_csv, barley_csv, oats_csv, maize_csv, rape_csv, sunflower_csv, rice_csv, fruits_and_berries_csv, vineyards_csv, olives_csv, country_list, output_csv):

    # This function summarizes the national annual production of 11 types of agricultural products
    wheat_df = get_annual_average_production(wheat_csv, 'wheat', country_list)
    rye_df = get_annual_average_production(rye_csv, 'rye', country_list)
    barley_df = get_annual_average_production(barley_csv, 'barley', country_list)
    oats_df = get_annual_average_production(oats_csv, 'oats', country_list)
    maize_df = get_annual_average_production(maize_csv, 'maize', country_list)
    rape_df = get_annual_average_production(rape_csv, 'rape', country_list)
    sunflower_df = get_annual_average_production(sunflower_csv, 'sunflower', country_list)
    rice_df = get_annual_average_production_rice(rice_csv, country_list)
    fruits_and_berries_df = get_annual_average_production(fruits_and_berries_csv, 'fruits_and_berries', country_list)
    vineyards_df = get_annual_average_production(vineyards_csv, 'vineyards', country_list)
    olives_df = get_annual_average_production(olives_csv, 'olives', country_list)

    synthesize_df = pd.DataFrame({
        'country_codes': country_list,
        'wheat_t': wheat_df.annual_production_t,
        'rye_t': rye_df.annual_production_t,
        'barley_t': barley_df.annual_production_t,
        'oats_t': oats_df.annual_production_t,
        'maize_t': maize_df.annual_production_t,
        'rape_t': rape_df.annual_production_t,
        'sunflower_t': sunflower_df.annual_production_t,
        'rice_t': rice_df.annual_production_t,
        'fruits_and_berries_t': fruits_and_berries_df.annual_production_t,
        'vineyards_t': vineyards_df.annual_production_t,
        'olives_t': olives_df.annual_production_t,
        })

    synthesize_df = synthesize_df.fillna(0)
    synthesize_df.to_csv(output_csv)


if __name__ == "__main__":
    synthesize(
        wheat_csv = snakemake.input.wheat, 
        rye_csv = snakemake.input.rye, 
        barley_csv = snakemake.input.barley, 
        oats_csv = snakemake.input.oats, 
        maize_csv = snakemake.input.maize, 
        rape_csv = snakemake.input.rape, 
        sunflower_csv = snakemake.input.sunflower, 
        rice_csv = snakemake.input.rice, 
        fruits_and_berries_csv = snakemake.input.fruits_and_berries, 
        vineyards_csv = snakemake.input.vineyards, 
        olives_csv = snakemake.input.olives,
        country_list = snakemake.params[0], 
        output_csv = snakemake.output[0]
        )