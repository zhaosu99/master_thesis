import pandas as pd
import pycountry


def get_wanted_df(raw_df, incgrp, building, deg_urb, geo):

    wanted_df = raw_df.copy()[(raw_df.incgrp == incgrp) & (raw_df.building == building) & \
        (raw_df.deg_urb == deg_urb) & (raw_df.geo == geo)]
    return wanted_df


def from_code3_to_code2(code3):
    if code3 == 'GRC':
        code2 = 'EL'
    elif code3 == 'GBR':
        code2 = 'UK'
    else:
        code2 = pycountry.countries.get(alpha_3=code3).alpha_2
    return code2


def calculate_sfh_ratios(input_csv, country_codes, output_csv):

    input_df = pd.read_csv(input_csv)
    sfh_ratio_list = []
    for code3 in country_codes:
        code2 = from_code3_to_code2(code3)
        # Current dataset doesn't include Bosnia and Herzegovina and Montenegro, their data are set to be the data of their neighbor Croatia and Albaina
        if code2 == 'BA':
            code2 = 'HR'
        elif code2 == 'ME':
            code2 = 'AL'

        sfh_df = get_wanted_df(input_df, 'TOTAL', 'HOUSE', 'TOTAL', code2)
        year = sfh_df.TIME_PERIOD.max()
        sfh_value = float(sfh_df.loc[lambda x: x.TIME_PERIOD == year, 'OBS_VALUE'])
        mfh_df = get_wanted_df(input_df, 'TOTAL', 'FLAT', 'TOTAL', code2)
        year = mfh_df.TIME_PERIOD.max()
        mfh_value = float(mfh_df.loc[lambda x: x.TIME_PERIOD == year, 'OBS_VALUE'])

        sfh_ratio = sfh_value / (sfh_value + mfh_value)
        sfh_ratio_list.append(sfh_ratio)
    
    pd.DataFrame({
        'code': country_codes,
        'sfh_ratio': sfh_ratio_list
        }
    ).to_csv(output_csv)


if __name__ == "__main__":
    calculate_sfh_ratios(
        input_csv = snakemake.input[0],
        country_codes = snakemake.params[0], 
        output_csv = snakemake.output[0]
        )


