import pandas as pd
import numpy as np

def calculate_average_capacity_factor(input_csv_list, ratio_list, output_csv):

    #This function creates weighted average of monthly capacity factors.
    #The weightings of West-East, South-Flat and North rooftops are given in the 'ratio_list
    month_cols = [str(i).zfill(2) for i in range(1,13)]
    column_names = ["site_id"] + month_cols
    output_dataframe = pd.DataFrame(columns=column_names).reset_index(drop=True)

    file_number = len(input_csv_list)
    example_df = pd.read_csv(input_csv_list[0])
    output_dataframe["site_id"] = example_df["site_id"]

    row_number = example_df.shape[0]
    for m in range(row_number): 
        cf_list = []
        for n in range(file_number):
            cf_list.append(np.array(pd.read_csv(input_csv_list[n]).iloc[m, 1:]) * float(ratio_list[n]))
         
        output_dataframe.loc[m, month_cols] = np.sum(cf_list, axis=0)
    
    output_dataframe.to_csv(output_csv)


if __name__ == "__main__":
    calculate_average_capacity_factor(
        input_csv_list = snakemake.input, 
        ratio_list = snakemake.params, 
        output_csv = snakemake.output[0]
        )