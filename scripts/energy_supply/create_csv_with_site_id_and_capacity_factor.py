import numpy as np
import xarray as xr
import pandas as pd


def calculate_monthly_capacity_factors(input_nc_list, output_csv_list):

    '''
    This function calculates the average monthly capacity factor from existing renewables.ninja silumations for each month.

    The three variablles from nc files: geometry, site id, capacity factors are grouped into two series of files:
    firstly, the site id geojson file with geometry, site id (done in the "create_voronoi_polygons_with_site_id" function);
    secondly, the capacity factor csv with site id, capacity factors (this function).

    This splitting is for accelerating calculation, since the geometry, site id columns are the same for all onshore simulations
    '''

    for i in range(len(input_nc_list)):

        # Create a dataframe with only one column "site_id" and set it as the index column
        raw_data = xr.open_dataset(input_nc_list[i])
        site_id_list = [int(i) for i in raw_data.coords['site_id']]
        create_dataframe = pd.DataFrame({"site_id": site_id_list}).set_index("site_id")

        # Calculate the average monthly capacity factor from January to December
        year_list = [str(i) for i in range(2000, 2019)]
        month_list = [str(i).zfill(2) for i in range(1, 13)]
        for month in month_list:
            locals()[month], temperal_ave_list, temperal_len_list = [], [], []
            for id in site_id_list:
                for year in year_list:
                    capacity_factor_list = raw_data.sel(time="{}-{}".format(year, month), site_id=id).electricity
                    capacity_factor_list = np.array(capacity_factor_list).tolist()
                    temperal_len_list.append(len(capacity_factor_list))
                    temperal_ave_list.append(sum(capacity_factor_list) / len(capacity_factor_list))
                locals()[month].append(np.sum(np.array(temperal_len_list) * np.array(temperal_ave_list)) / sum(temperal_len_list))
            create_dataframe[month] = locals()[month]
        create_dataframe.to_csv(output_csv_list[i])


if __name__ == "__main__":
    calculate_monthly_capacity_factors(
        input_nc_list = snakemake.input, 
        output_csv_list = snakemake.output
        )
