import pandas as pd


def concat_csv_list(input_csv_list, output_csv):

    temp_list = [pd.read_csv(i) for i in input_csv_list]
    pd.concat(temp_list).reset_index(drop=True).to_csv(output_csv)


if __name__ == "__main__":
    concat_csv_list(
        input_csv_list = snakemake.input, 
        output_csv = snakemake.output[0]
        )