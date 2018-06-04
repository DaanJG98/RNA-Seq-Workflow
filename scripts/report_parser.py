# Script parsing RNA-seq-workflow result and report them

import re
from snakemake.utils import report

# return the amount of pubmed ids for a gene.
# if empty, IndexError, it returns 0
def n_pubmed_ids(item):
    try:
        return len(item.split("\n\n")[1].split(":")[1].split(","))
    except IndexError:
        return 0

def prepare_gene_info(gene_ids, ncbi_ids, gene_info, sequence_file, kegg_ids, ortholog_file):
    report_list = []
    # Open all input files and attach their content to a list bases on index
    with open(gene_ids) as id_file:
        for idx, id in enumerate(id_file):
            if not id.split():
                continue
            else:
                report_list.append(id+"---------------------"+"\n")

    with open(ncbi_ids) as ncbi_file:
        for idx, ncbi_id in enumerate(ncbi_file):
            if not ncbi_id.split():
                continue
            else:
                report_list[idx]+="NCBI gene id: "+ncbi_id+"\n"

    with open(gene_info) as gene_info_file:
        for idx, gene_info in enumerate(gene_info_file):
            splitted_line = gene_info.split("\t")
            pubmed_ids = splitted_line[4]
            description = splitted_line[5]
            report_list[idx]+="Pubmed ids: "+pubmed_ids+"\n\n"
            report_list[idx]+="Gene function: "+description+"\n"

    with open(kegg_ids) as kegg_ids:
        for idx, kegg_id in enumerate(kegg_ids):
            report_list[idx]+="KEGG ids: "+kegg_id+"\n"

    with open(ortholog_file) as ortholog_file:
        orthologs = ortholog_file.read().split(">")
        for idx, ortholog in enumerate(orthologs):
            if not ortholog.strip():
                continue
            else:
                ortholog_ids = ortholog.rstrip().split("\n")[1:] # skip first item with is the original id of the gene
                # index minus 1 because the first split creates an empty item in list orthologs
                report_list[idx-1]+="Orthologs: "+", ".join(ortholog_ids)+"\n\n"

    with open(sequence_file) as sequences:
        for idx, sequence in enumerate(sequences):
            sequence = sequence.split("\t")[1]
            refac_seq = re.sub("(.{85})", "\\1\n", sequence, 0, re.DOTALL)
            report_list[idx]+="Sequence: "+refac_seq+"\n\n"

    # collect all data in a string to show in the report, sorted by amount of pubmed ids using the n_pubmed_ids function
    report_string = "\n".join(sorted(report_list, reverse=True, key=n_pubmed_ids))

    return report_string

def prepare_pubmed_info(pubmed_file):
    with open(pubmed_file) as pubmed_file:
        pubmed_string = ""
        for line in pubmed_file:
            pubmed_string += line.replace("\t", ": ")+"\n"

    return pubmed_string

def create_report(report_string, pubmed_string, image_files):
    image_string = """"""
    for idx in range(len(image_files)):
        i = str(idx)
        image_string += ".. image:: {image_files["+i+"]}"
        if idx != len(image_files)-1:
            image_string += """
    """

    report("""
    RNA-seq gene information report
    ===================================================
    {report_string}

    PubMed
    ----------------------------------------------------
    {pubmed_string}

    """+image_string, snakemake.output[0], metadata="Authors: Daan Gilissen and Koen Rademaker (support: daangilissen@gmail.com)")

report_string = prepare_gene_info(snakemake.input[0], snakemake.input[1], snakemake.input[2], snakemake.input[3], snakemake.input[4], snakemake.input[6])
pubmed_string = prepare_pubmed_info(snakemake.input[5])

create_report(report_string, pubmed_string, snakemake.input[7:])
