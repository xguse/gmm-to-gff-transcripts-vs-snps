# Gene set comparisons report (linked SNPs vs original SNPs)

{{ metadata }}


# Comparison Tables

## Official Genes
__Counts:__

{{ official_genes_set_counts }}
__Gene Table:__

{{ official_genes_set_table }}

## Novel Transcripts

__Counts:__

{{ novel_tx_set_counts }}

__Gene Table:__
{{ novel_tx_set_table }}

# Genome Information

## Official Genes

{% for gene in genes %}
    ### {{ gene.name }}
    __Genome Region__

    ![    ]({{ gene.region }})

{% endfor %}


## Novel Transcripts

{% for tx in txs %}
    ### {{ tx.name }}
    __Genome Region__

    ![    ]({{ tx.region }})

{% endfor %}
