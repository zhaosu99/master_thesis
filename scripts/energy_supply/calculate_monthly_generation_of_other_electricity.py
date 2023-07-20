import pandas as pd


def calculate_monthly_generation(input_generation_csv_list, average_month_generation_csv):

    # This function calculates the national monthly electricity generation of hydro, nuclear and geothermal power plants
    csv_list = [pd.read_csv(i) for i in input_generation_csv_list]
    monthly_generation = pd.concat(csv_list)

    monthly_generation.columns = ['country', 'month', 'geothermal_MW', 'hydro_MW_1', 'hydro_MW_2', 'marine_MW', 'nuclear_MW']
    monthly_generation = monthly_generation.replace("n/e", 0)
    monthly_generation['month'] = monthly_generation['month'].apply(lambda x: x[3: 5])
    monthly_generation['hydro_MW'] = monthly_generation.apply(lambda x: float(x['hydro_MW_1']) + float(x['hydro_MW_2']), axis=1)
    monthly_generation = monthly_generation.loc[:, ['month', 'hydro_MW', 'nuclear_MW', 'geothermal_MW']]
    monthly_generation = monthly_generation.groupby('month').mean()

    # Many countries don't have all of the three above-mentioned generation techs
    # To avoid columns discrepancy, the missed values are set to 0
    techs = ['hydro_MW', 'nuclear_MW', 'geothermal_MW']
    for i in techs:
        if i not in monthly_generation.columns:
            monthly_generation[i] = 0

    # Calculate monthly generation
    A_month = ['01', '03', '05', '07', '08', '10', '12']
    B_month = ['04', '06', '09', '11']   
    monthly_generation.loc[A_month, techs] = monthly_generation.loc[A_month, techs] * 31 * 24
    monthly_generation.loc[B_month, techs] = monthly_generation.loc[B_month, techs] * 30 * 24
    monthly_generation.loc['02', techs] = monthly_generation.loc['02', techs] * 28 * 24
    monthly_generation.columns = ['hydro_MWh', 'nuclear_MWh', 'geothermal_MWh']

    monthly_generation.to_csv(average_month_generation_csv)


def create_empty_dataframe(average_month_generation_csv):

    # Some countries don't have these three types of power plants
    month_list = [str(i).zfill(2) for i in range(1, 13)]
    empty_dataframe = pd.DataFrame()
    empty_dataframe['month'] = month_list
    empty_dataframe['hydro_MWh'] = 0
    empty_dataframe['nuclear_MWh'] = 0
    empty_dataframe['geothermal_MWh'] = 0

    empty_dataframe.to_csv(average_month_generation_csv)


def process_Albania_data(ALB_seasonal_generation, average_month_generation_csv):

    # Albania has dominent hydro power generaton, which is not included in the ENTSOE dataset, here its monthly hydro generation is infered from its own national data
    alb = pd.read_excel(ALB_seasonal_generation)
    alb = list(alb.iloc[9 , 2: 18])

    spring_index = [i for i in range(0,15,4)]
    summer_index = [i for i in range(1,15,4)]
    autumn_index = [i for i in range(2,15,4)]
    winter_index = [i for i in range(3,15,4)]

    spring_average = float(sum([alb[i] for i in spring_index]) / len([alb[i] for i in spring_index]))
    summer_average = float(sum([alb[i] for i in summer_index]) / len([alb[i] for i in summer_index]))
    autumn_average = float(sum([alb[i] for i in autumn_index]) / len([alb[i] for i in autumn_index]))
    winter_average = float(sum([alb[i] for i in winter_index]) / len([alb[i] for i in winter_index]))

    month_list = [str(i).zfill(2) for i in range(1, 13)]
    empty_dataframe = pd.DataFrame()
    empty_dataframe['month'] = month_list
    empty_dataframe['hydro_MWh'] = [spring_average] * 3 + [summer_average] * 3 + [autumn_average] * 3 + [winter_average] * 3
    empty_dataframe['nuclear_MWh'] = 0
    empty_dataframe['geothermal_MWh'] = 0

    empty_dataframe.to_csv(average_month_generation_csv)


def calculate_monthly_generation_of_other_electricity(input_generation_csv_list, ALB_seasonal_generation, average_month_generation_csv):

    if len(input_generation_csv_list) > 0:
        calculate_monthly_generation(input_generation_csv_list, average_month_generation_csv)
    elif average_month_generation_csv[-7: -4] == "ALB":
        process_Albania_data(ALB_seasonal_generation, average_month_generation_csv)
    else:
        create_empty_dataframe(average_month_generation_csv)


if __name__ == "__main__":
    calculate_monthly_generation_of_other_electricity(
        input_generation_csv_list = snakemake.input.monthly_generation_input, 
        ALB_seasonal_generation = snakemake.input.ALB_seasonal_generation,
        average_month_generation_csv = snakemake.output[0]
        )