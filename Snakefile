rule all:
    input:
        "results/minimizer.json",


rule download_references:
    input:
        config=config["config_file"],
    output:
        reference_genome="results/references.fasta",
    shell:
        """
        python scripts/download_references.py \
            --config-file {input.config} \
            --segment-file {output.reference_genome} \
        """


rule create_minimizer:
    input:
        reference="results/references.fasta",
        config=config["config_file"],
    output:
        minimizer="results/minimizer.json",
    shell:
        """
        python scripts/create_minimizer_index.py \
            --references-fasta {input.reference} \
            --minimizer-json {output.minimizer} \
            --config-file {input.config}
        """
