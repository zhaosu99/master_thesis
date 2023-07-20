import pandas as pd


def concat_continental_national_energy_inventory(input_csv, output_csv):

    def add_national_id_column(input_df):

        if 'regional_id' in input_df.columns:
            input_df['national_id'] = input_df['regional_id'].apply(lambda x: x[: 3])
            input_df = input_df.drop('regional_id', axis=1)
        elif 'municipal_id' in input_df.columns:
            input_df['national_id'] = input_df['municipal_id'].apply(lambda x: x[: 3])
            input_df = input_df.drop('municipal_id', axis=1)
        else:
            print('Regional_id and municipal_id are not in columns')
        return input_df
    
    result_df = add_national_id_column(pd.read_csv(input_csv))
    result_df = result_df.groupby(by=['national_id']).sum().drop(columns=['Unnamed: 0', 'index'], axis=1)

    result_df.to_csv(output_csv)


if __name__ == "__main__":
    concat_continental_national_energy_inventory(
        input_csv = snakemake.input[0], 
        output_csv = snakemake.output[0]
        )