from snakemake.utils import min_version
PANDOC = "pandoc --filter pantable --filter pandoc-fignos --filter pandoc-tablenos --citeproc"

configfile: "config/default.yaml"

include: "rules/energy_supply.smk"
include: "rules/energy_demand.smk"
include: "rules/visualization.smk"


min_version("7.8")

wildcard_constraints:
    gid = "({gid_list})".format(gid_list="|".join((f"({gid})" for gid in config["country_codes"])))

onsuccess:
    if "email" in config.keys():
        shell("echo "" | mail -s 'mypaper succeeded' {config[email]}")
onerror:
    if "email" in config.keys():
        shell("echo "" | mail -s 'mypaper failed' {config[email]}")

rule all:
    message: "Run entire analysis."
    input:
        "build/result/visualization/pics/continental_municipalities_energy_inventory.jpg",
        "build/result/visualization/pics/continental_municipalities_energy_inventory_per_capita.jpg",
        "build/result/visualization/pics/continental_regions_energy_inventory.jpg",
        "build/result/visualization/pics/continental_regions_energy_inventory_per_capita.jpg",
        "build/result/visualization/pics/continental_countries_energy_inventory.jpg",
        "build/result/visualization/pics/continental_countries_energy_inventory_per_capita.jpg",
        "build/result/visualization/pics/national_energy_structure.jpg",
        "build/result/visualization/pics/municipal_boxplots.jpg",
        "build/result/visualization/pics/regional_boxplots.jpg"

rule run:
    message: "Runs the demo model."
    params:
        slope = config["slope"],
        x0 = config["x0"]
    output: "build/results.pickle"
    conda: "envs/default.yaml"
    script: "scripts/model.py"


rule plot:
    message: "Visualises the demo results."
    input:
        results = rules.run.output
    output: "build/plot.png"
    conda: "envs/default.yaml"
    script: "scripts/vis.py"


def pandoc_options(wildcards):
    suffix = wildcards["suffix"]
    if suffix == "html":
        return "--self-contained --to html5"
    elif suffix == "pdf":
        return "--pdf-engine weasyprint"
    elif suffix == "docx":
        return []
    else:
        raise ValueError(f"Cannot create report with suffix {suffix}.")


rule report:
    message: "Compile report.{wildcards.suffix}."
    input:
        "report/literature.yaml",
        "report/report.md",
        "report/pandoc-metadata.yaml",
        "report/apa.csl",
        "report/reset.css",
        "report/report.css",
        rules.plot.output
    params: options = pandoc_options
    output: "build/report.{suffix}"
    wildcard_constraints:
        suffix = "((html)|(pdf)|(docx))"
    conda: "envs/report.yaml"
    shadow: "minimal"
    shell:
        """
        cd report
        ln -s ../build .
        {PANDOC} report.md  --metadata-file=pandoc-metadata.yaml {params.options} \
        -o ../build/report.{wildcards.suffix}
        """


rule dag:
     message: "Plot dependency graph of the workflow."
     output:
         dot = "build/dag.dot",
         pdf = "build/dag.pdf"
     conda: "envs/dag.yaml"
     shell:
         """
         snakemake --rulegraph > {output.dot}
         dot -Tpdf -o {output.pdf} {output.dot}
         """


rule clean: # removes all generated results
    message: "Remove all build results but keep downloaded data."
    run:
         import shutil

         shutil.rmtree("build")
         print("Data downloaded to data/ has not been cleaned.")


rule test:
    conda: "envs/test.yaml"
    output: "build/test-report.html"
    shell:
        "py.test --html={output} --self-contained-html"
