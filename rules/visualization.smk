rule download_national_countries:
    message: "Download ADM 0 level for {wildcards.gid}"
    params: "https://www.geoboundaries.org/data/geoBoundaries-3_0_0/{gid}/ADM0/geoBoundaries-3_0_0-{gid}-ADM0.geojson"
    output:
        "data/geoboundaries/{gid}_countries.geojson"
    conda: "../envs/default.yaml"
    shell:
        "curl -o {output} '{params}'"


rule add_population_amount_to_countries:
    message: "Add the population to {wildcards.gid}"
    input:
        "data/geoboundaries/{gid}_countries.geojson",
        "data/geoboundaries/NLD_simplified/geoBoundaries-NLD-ADM0_simplified.geojson",
        "data/geoboundaries/NOR_simplified/geoBoundaries-NOR-ADM0_simplified.geojson",
        "temp/population_density/pop_3035.tif",
    output:
        "temp/geoboundaries/{gid}_countries_add_pop_amount.geojson"
    conda: "../envs/default.yaml"
    script:
        "../scripts/visualization/add_population_amount_to_countries.py"


rule add_population_amount_to_regions:
    message: "Add the population to each region in {wildcards.gid}"
    input:
        "temp/geoboundaries/{gid}_regions.geojson",
        "temp/population_density/pop_3035.tif"
    output:
        "temp/geoboundaries/{gid}_regions_add_pop_amount.geojson"
    conda: "../envs/default.yaml"
    script:
        "../scripts/energy_demand/add_population_amount_to_units.py"


rule plot_potential_demand_1:
    input:
        "data/others/potentials1.csv"
    output:
        "build/result/visualization/pics/potential_demand_1.jpg"
    conda: "../envs/default.yaml"
    script:
        "../scripts/visualization/plot_potential_demand_1.py"


rule plot_potential_demand_2:
    input:
        "data/others/potentials2.csv"
    output:
        "build/result/visualization/pics/potential_demand_2.jpg"
    conda: "../envs/default.yaml"
    script:
        "../scripts/visualization/plot_potential_demand_2.py"


rule preprocess_national_visualization_data:
    input:
        inventory = "build/national/continental_national_energy_inventory.csv",
        population = expand("temp/geoboundaries/{gid}_countries_add_pop_amount.geojson", gid=config['country_codes'])
    output:
        output_csv = "build/result/visualization/data/continental_national_energy_inventory.csv"
    conda: "../envs/default.yaml"
    script:
        "../scripts/visualization/preprocess_data.py"


rule preprocess_regional_visualization_data:
    input:
        inventory = "build/regional/continental_regional_energy_inventory.csv",
        population = expand("temp/geoboundaries/{gid}_regions_add_pop_amount.geojson", gid=config['country_codes'])
    output:
        output_csv = "build/result/visualization/data/continental_regional_energy_inventory.csv"
    conda: "../envs/default.yaml"
    script:
        "../scripts/visualization/preprocess_data.py"


rule preprocess_municipal_visualization_data:
    input:
        inventory = "build/municipal/continental_municipal_energy_inventory.csv",
        population = expand("temp/geoboundaries/{gid}_municipalities_add_pop_amount.geojson", gid=config['country_codes'])
    output:
        output_csv = "build/result/visualization/data/continental_municipal_energy_inventory.csv"
    conda: "../envs/default.yaml"
    script:
        "../scripts/visualization/preprocess_data.py"


rule get_municipal_boxplots:
    input:
        "build/result/visualization/data/continental_municipal_energy_inventory.csv"
    params: "municipal_id"
    output:
        "build/result/visualization/pics/municipal_boxplots.jpg"
    conda: "../envs/default.yaml"
    script:
        "../scripts/visualization/boxplot.py"


rule get_regional_boxplots:
    input:
        "build/result/visualization/data/continental_regional_energy_inventory.csv"
    params: "regional_id"
    output:
        "build/result/visualization/pics/regional_boxplots.jpg"
    conda: "../envs/default.yaml"
    script:
        "../scripts/visualization/boxplot.py"


rule plot_national_energy_structure:
    input:
        "build/result/visualization/data/continental_national_energy_inventory.csv"
    output:
        "build/result/visualization/pics/national_energy_structure.jpg"
    conda: "../envs/default.yaml"
    script:
        "../scripts/visualization/plot_national_energy_structure.py"


rule continental_municipalities_energy_inventory_visualization:
    input:
        "temp/geoboundaries/continental_municipalities.geojson",
        "build/result/visualization/data/continental_municipal_energy_inventory.csv" 
    output:
        "build/result/visualization/pics/continental_municipalities_energy_inventory.jpg",
        "build/result/visualization/pics/continental_municipalities_energy_inventory_per_capita.jpg",
    params: 'municipal_id'
    conda: "../envs/default.yaml"
    script:
        "../scripts/visualization/continental_energy_inventory_visualization.py"


rule continental_regions_energy_inventory_visualization:
    input:
        "temp/geoboundaries/continental_regions.geojson",
        "build/result/visualization/data/continental_regional_energy_inventory.csv"
    output:
        "build/result/visualization/pics/continental_regions_energy_inventory.jpg",
        "build/result/visualization/pics/continental_regions_energy_inventory_per_capita.jpg"
    conda: "../envs/default.yaml"
    script:
        "../scripts/visualization/continental_energy_inventory_visualization.py"


rule get_continental_countries:
    input:
        expand("temp/geoboundaries/{gid}_countries_add_pop_amount.geojson", gid=config['country_codes']),
    output:
        "temp/geoboundaries/continental_countries.geojson"
    conda: "../envs/default.yaml"
    script:
        "../scripts/visualization/get_continental_countries.py"


rule continental_countries_energy_inventory_visualization:
    input:
        "temp/geoboundaries/continental_countries.geojson",
        "build/result/visualization/data/continental_national_energy_inventory.csv"
    output:
        "build/result/visualization/pics/continental_countries_energy_inventory.jpg",
        "build/result/visualization/pics/continental_countries_energy_inventory_per_capita.jpg"
    conda: "../envs/default.yaml"
    script:
        "../scripts/visualization/continental_energy_inventory_visualization.py"



