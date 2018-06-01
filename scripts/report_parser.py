# Script parsing RNA-seq-workflow result and report them

import sys
import re
from snakemake.utils import report

def prepare_gene_info(gene_ids, ncbi_ids, gene_info, sequence_file, kegg_ids):
    report_list = []
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
    with open(sequence_file) as sequences:
        for idx, sequence in enumerate(sequences):
            sequence = sequence.split("\t")[1]
            refac_seq = re.sub("(.{85})", "\\1\n", sequence, 0, re.DOTALL)
            report_list[idx]+="Sequence: "+refac_seq+"\n\n"

    report_string = "\n".join(report_list)

    return report_string

def prepare_pubmed_info(pubmed_file):
    with open(pubmed_file) as pubmed_file:
        pubmed_string = ""
        for line in pubmed_file:
            pubmed_string += line.replace("\t", ": ")

    return pubmed_string

# def prepare_image_string(image_files):
#     image_string = """"""
#     for idx,image in enumerate(image_files):
#         image_string += """.. image:: {image_files[idx]}"""
#     return image_string

def create_report(report_string, pubmed_string, image_files):
    image_string = """"""
    for idx in range(len(image_files)):
        i = str(idx)
        image_string += ".. image:: {image_files["+i+"]}"
        if idx != len(image_files)-1:
            image_string += """
    """

    report("""
    Example content report
    ===================================================
    {report_string}

    PubMed
    ----------------------------------------------------
    {pubmed_string}



    """+image_string, snakemake.output[0])


report_string = prepare_gene_info(snakemake.input[0], snakemake.input[1], snakemake.input[2], snakemake.input[3], snakemake.input[4])
pubmed_string = prepare_pubmed_info(snakemake.input[5])
# image_string = prepare_image_string(snakemake.input[6:])

create_report(report_string, pubmed_string, snakemake.input[6:])
