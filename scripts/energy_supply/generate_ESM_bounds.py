import rasterio as rio
import pandas as pd


def generate_ESM_bounds(ESM_list, output_csv):

    # This function creates a list containing the bounds of ESM rasters in the ESM list
    bounds_list = []
    index_list = []
    for i in range(len(ESM_list)):
        n = tuple(int(i) for i in list(rio.open(ESM_list[i]).bounds))
        bounds_list.append(n)
        index_list.append(i)

    index = pd.MultiIndex.from_tuples(bounds_list, names=['min_x', 'min_y', 'max_x', 'max_y'])
    ESM_bounds_dataframe = pd.DataFrame({"index": index_list}, index=index)
    
    ESM_bounds_dataframe.to_csv(output_csv)


if __name__ == "__main__":
    generate_ESM_bounds(
        ESM_list = snakemake.input.ESM_list, 
        output_csv = snakemake.output[0]
        )

