import pandas as pd


def concat_continental_energy_inventory(municipal_or_regional, input_csv_list, output_csv):

    wanted_cols = [municipal_or_regional, 'onshore_wind_MWh', 'offshore_wind_MWh', 'open_field_solar_PV_MWh',
       'rooftop_solar_PV_MWh', 'hydro_MWh', 'nuclear_MWh', 'biomass_MWh', 'geothermal_MWh', 
       'transport_e', 'industry_e', 'industry_h', 'building_e']
    
    df_list = []
    for i in input_csv_list:
        sel_df = pd.read_csv(i)
        for col in wanted_cols:
            if col not in sel_df.columns:
                sel_df[col] = 0
        df_list.append(sel_df.loc[:, wanted_cols])
    
    pd.concat(df_list).reset_index(drop=True).to_csv(output_csv)


if __name__ == "__main__":
    concat_continental_energy_inventory(
        municipal_or_regional = snakemake.params[0], 
        input_csv_list = snakemake.input, 
        output_csv = snakemake.output[0]
        )