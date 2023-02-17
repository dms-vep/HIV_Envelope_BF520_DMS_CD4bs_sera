"""Top-level ``snakemake`` file that runs analysis."""


import os


configfile: "config.yaml"


# include `dms-vep-pipeline` pipeline Snakemake file
include: os.path.join(config["pipeline_path"], "pipeline.smk")

antibody_escape_pdbs = config["antibody_escape_PDBs"]

Env_chains_by_pdb = config["env_chains_by_PDB"]

rule all:
    input:
        variant_count_files,
        rules.check_adequate_variant_counts.output.passed,
        antibody_escape_files,
        (
            [config["muteffects_observed"], config["muteffects_latent"]]
            if len(func_selections)
            else []
        ),
        config["docs"],
        expand(
            os.path.join(config["escape_PDB_structures"], "{antibody}" + "_epitope_1.pdb"),
            antibody=antibody_escape_pdbs.keys(),
        )


# Arbitrary other rules should be added here
rule site_numbering_map:
    """Map sequential numbering of protein in experiments to standard reference."""
    input:
        prot=config["gene_sequence_protein"],
        reference_site_regions=config["reference_site_regions"],
    output:
        reference="results/site_numbering/numbering_reference.fa",
        alignment="results/site_numbering/alignment.fa",
        to_align="results/site_numbering/to_align.fa",
        site_numbering_map=config["site_numbering_map"],
    params:
        numbering_reference_accession=config["numbering_reference_accession"],
    log:
        os.path.join(config["logdir"], "site_numbering_map.txt"),
    conda:
        "dms-vep-pipeline/environment.yml"
    script:
        "scripts/site_numbering_map.py"
        
rule spatial_distances:
    """Get spatial distances from PDB."""
    input: 
        pdb=config["PDB"],
    output:
        csv=config["spatial_distances"],
    params:
        target_chains=["A", "B", "C"],
    log:
        log=os.path.join(config["logdir"], "spatial_distances.txt"),
    conda:
        os.path.join(config["pipeline_path"], "environment.yml")
    script:
        "scripts/spatial_distances.py"
        
rule validation_ICs:
    """Get ``polyclonal`` predicted ICs for validated mutations."""
    input:
        config["validation_ics"],
        [
            os.path.join(config["escape_dir"], f"{antibody}.pickle")
            for antibody in pd.read_csv(config["validation_ics"])["antibody"].unique()
        ],
        nb="notebooks/validation_ICs.ipynb",
    output:
        nb="results/notebooks/validation_ICs.ipynb",
    log:
        os.path.join(config["logdir"], "validationi_ICs.txt"),
    shell:
        "papermill {input.nb} {output.nb} &> {log}"

rule mutations_vs_natural_sequences:
    """Compare mutation effects vs natural sequence variation."""
    input:
        config["natural_sequence_data"],
        config["natural_sequence_counts"],
        config["muteffects_observed"],
        nb="notebooks/mutation_effects_versus_natural_frequencies.ipynb",
    output:
        nb="results/notebooks/mutation_effects_versus_natural_frequencies.ipynb",
    log:
        os.path.join(config["logdir"], "mutation_effects_versus_natural_frequencies.txt"),
    shell:
        "papermill {input.nb} {output.nb} &> {log}"

rule color_PDB_structures:
    """Assign b factor values to PDB structures based on escape"""
    input: 
        average_escape_model = rules.avg_antibody_escape.output.avg_escape,
        input_pdb_file = lambda wc: os.path.join(
            config["PDB_structures"], 
            antibody_escape_pdbs[wc.antibody]
        ),
    params: 
        env_chains = lambda wc: Env_chains_by_pdb[antibody_escape_pdbs[wc.antibody]],
        output_file_sub_name = os.path.join(config["escape_PDB_structures"], "{antibody}"),
    output:
        output_pdb_file_name = os.path.join(config["escape_PDB_structures"], "{antibody}" + "_epitope_1.pdb")
    log:
        os.path.join(config["logdir"], "antibody_escape_pdbs_{antibody}.txt"),
    script:
        "scripts/color_PDB_structures.py"
        

# Add any extra data/results files for docs with name: file
extra_data_files = {
    "sequential to reference site numbering": config["site_numbering_map"],
}
for antibody in config["antibody_escape_PDBs"]:
    extra_data_files[f"PDB with {antibody} escape values as b factors"] = os.path.join(config["escape_PDB_structures"], f"{antibody}" + "_epitope_1.pdb")

# If you add rules with "nb" output that have wildcards, specify the rule name
# and subindex titles for the wildcards as in "docs.smk" for `nb_rule_wildcards`
# and `subindex_titles`
    
# include `dms-vep-pipeline` docs building Snakemake file
include: os.path.join(config["pipeline_path"], "docs.smk")
