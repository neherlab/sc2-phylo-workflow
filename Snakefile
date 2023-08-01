nextcladev3 = "~/Projects/nextstrain/nextclade_dev/target/release/nextclade"
genemapv3 = "~/Projects/nextstrain/nextclade_dev/target/release/genemap"

genes = ["E", "M", "N", "ORF1a", "ORF1b", "ORF3a", "ORF6", "ORF7a", "ORF7b", "ORF8", "ORF9b", "S"]


rule get_nextclade_dataset:
    output:
        directory("sars-cov-2_dataset")
    shell:
        """
        {nextcladev3} dataset get --name sars-cov-2 --output-dir {output}
        """

rule nextclade_annotation:
    input:
        dataset = directory("sars-cov-2_dataset")
    output:
        "data/annotation.json"
    shell:
        """
        {genemapv3} {input.dataset}/genemap.gff --json > {output}
        """

rule fetch_data:
    output:
        sequences = "data/sequences.fasta.xz",
        metadata = "data/metadata.tsv",
    shell:
        """
        curl https://data.nextstrain.org/files/ncov/open/global/sequences.fasta.xz -o {output.sequences}
        curl https://data.nextstrain.org/files/ncov/open/global/metadata.tsv.xz | xz -d > {output.metadata}
        """

rule fetch_config:
    output:
        auspice_config = "data/auspice_config.json",
        lat_longs = "data/lat_longs.tsv"
    shell:
        """
        curl https://raw.githubusercontent.com/nextstrain/ncov/master/defaults/auspice_config.json -o {output.auspice_config}
        curl https://raw.githubusercontent.com/nextstrain/ncov/master/defaults/lat_longs.tsv -o {output.lat_longs}
        """


rule nextclade:
    input:
        sequences = "data/sequences.fasta.xz",
        dataset = directory("sars-cov-2_dataset")
    output:
        directory("results")
    threads: 4
    shell:
        """
        {nextcladev3} run -D {input.dataset} -j {threads} --output-all {output} {input.sequences}
        """

rule prune:
    input:
        tree = "results/nextclade.nwk",
        nextclade_tsv = "results/nextclade.tsv"
    output:
        tree = "results/tree_raw.nwk"
    run:
        import pandas as pd
        from Bio import Phylo

        tree = Phylo.read(input.tree, "newick")
        metadata = pd.read_csv(input.nextclade_tsv, sep="\t", index_col='seqName')
        for leaf in tree.get_terminals():
            if leaf.name not in metadata.index or metadata.loc[leaf.name, 'qc.overallStatus'] == 'bad':
                tree.prune(leaf)

        for clade in tree.find_clades():
            clade.branch_length /= 29903

        Phylo.write(tree, output.tree, "newick", format_branch_length="%1.8f")

rule refine:
    input:
        tree = "results/tree_raw.nwk",
        metadata = "data/metadata.tsv",
        alignment = "results/nextclade.aligned.fasta"
    output:
        tree = "results/tree.nwk",
        node_data = "results/branch_length.json"
    shell:
        """
        augur refine --tree {input.tree} --metadata {input.metadata} --alignment {input.alignment} \
                --timetree --coalescent opt --date-confidence \
                --clock-rate 0.0006 --clock-std-dev 0.0002 \
                --stochastic-resolve --keep-root \
                --output-node-data {output.node_data} --output-tree {output.tree}
        """

rule ancestral:
    input:
        tree = "results/tree.nwk",
        alignment = "results/nextclade.aligned.fasta",
        annotation = "data/annotation.json"
    output:
        node_data = "results/mutations.json"
    params:
        genes = genes,
        translations = "results/nextclade_gene_%GENE.translation.fasta"
    shell:
        """
        augur ancestral --tree {input.tree} --alignment {input.alignment} \
                --root-sequence sars-cov-2_dataset/reference.fasta --annotation {input.annotation} \
                --output-node-data {output.node_data} \
                --genes {params.genes} --translations {params.translations}
        """


rule export:
    input:
        tree = "results/tree.nwk",
        node_data = ["results/branch_length.json", "results/mutations.json"],
        metadata = "data/metadata.tsv",
        auspice_config = "data/auspice_config.json",
        lat_longs = "data/lat_longs.tsv"
    output:
        "auspice/tree.json"
    shell:
        """
        augur export v2 --tree {input.tree} --node-data {input.node_data} --metadata {input.metadata} \
                --lat-longs data/lat_longs.tsv --color-by-metadata Nextstrain_clade pango_lineage \
                --auspice-config data/auspice_config.json --output {output}
        """
