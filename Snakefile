mail="mymail@gmail.com"
from os import listdir
# rule all:
#     input:
#     "report.html"

rule filter_ids:
    input:
        "RNA-Seq-counts.txt"
    output:
        "data/RNA-seq-ids.txt"
    run:
        with open(output[0],"w") as out:
            with open(input[0]) as f:
                for line in f:
                    if line.startswith("#") or line.startswith("ID"):
                        next(f)
                    else:
                        out.write((line.split("\t")[0]+"\n"))

rule get_orthologs:
    input:
        "data/RNA-seq-ids.txt"
    output:
        "data/RNA-seq-orthologs.txt"
    shell:
        "./scripts/get_omadb_orthologs.sh args1 < {input} {output}"

rule get_ncbi_ids:
    input:
        "data/RNA-seq-ids.txt"
    output:
        "data/RNA-seq-ncbi-ids.txt"
    shell:
        "./scripts/get_ncbi_ids.sh args1 < {input} {output}"

rule get_gene_info:
    input:
        "data/RNA-seq-ncbi-ids.txt"
    output:
        "data/RNA-seq-gene-info.txt"
    run:
        from Bio import Entrez
        Entrez.email = mail
        with open(output[0],'a') as out:
            with open(input[0]) as file:
                for id in file:
                    if not id.strip():
                        continue
                    else:
                        handle = Entrez.efetch(db="gene", id=id, retmode="xml")
                        record = Entrez.read(handle)
                        organism = record[0]['Entrezgene_source']['BioSource']['BioSource_org']['Org-ref']['Org-ref_taxname']
                        genome = record[0]['Entrezgene_gene-source']['Gene-source']['Gene-source_src-str1']
                        start = record[0]['Entrezgene_locus'][0]['Gene-commentary_seqs'][0]['Seq-loc_int']['Seq-interval']['Seq-interval_from']
                        stop = record[0]['Entrezgene_locus'][0]['Gene-commentary_seqs'][0]['Seq-loc_int']['Seq-interval']['Seq-interval_to']
                        refs = record[0]['Entrezgene_locus'][0]['Gene-commentary_products'][0]['Gene-commentary_refs']
                        desciption = record[0]['Entrezgene_prot']['Prot-ref']['Prot-ref_name'][0]

                        pubmed = []
                        for pubmed_id in refs:
                            pubmed.append(pubmed_id['Pub_pmid']['PubMedId'])
                        out.write(genome+"\t"+start+"\t"+stop+"\t"+organism+"\t"+",".join(pubmed)+"\t"+desciption+"\n")

rule get_sequence:
    input:
        "data/RNA-seq-gene-info.txt"
    output:
        "data/RNA-seq-sequences.txt"
    run:
        import sys
        from Bio import Entrez
        from Bio.SeqUtils import GC
        Entrez.email=mail

        with open(output[0], "w") as out:
            with open(input[0]) as f:
                for info in f:
                    splittedinfo = info.split("\t")
                    genome = splittedinfo[0]
                    start = int(splittedinfo[1])+1
                    stop = int(splittedinfo[2])+1

                    handle = Entrez.efetch(db="nucleotide", id=genome, rettype="fasta", seq_start=start, seq_stop=stop, retmode="xml")
                    record = Entrez.read(handle)

                    gc = GC(record[0]['TSeq_sequence'])
                    gc = format(round(gc, 2))
                    out.write(gc+"\t"+record[0]['TSeq_sequence']+"\n")

rule conv_kegg_ids:
    input:
        "data/RNA-seq-ncbi-ids.txt"
    output:
        "data/RNA-seq-conv-kegg.txt"
    shell:
        "./scripts/conv_kegg.sh args1 < {input} {output}"

rule get_kegg_ids:
    input:
        "data/RNA-seq-conv-kegg.txt"
    output:
        "data/RNA-seq-kegg-ids.txt"
    run:
        from Bio.KEGG import REST
        with open(output[0],"w") as out:
            with open(input[0],"r") as f1:
                for line in f1:
                    splittedinfo = line.split("\t")
                    query = splittedinfo[1]
                    request = REST.kegg_get(query).readlines()

                    pathways = []
                    isPathway = False

                    for KEGG_report in iter(request):
                        if "PATHWAY" in KEGG_report:
                            isPathway=True
                            pathways.append(KEGG_report.strip("PATHWAY").lstrip().rstrip())
                            continue
                        if "BRITE" in KEGG_report or "MODULE" in KEGG_report or "POSITION" in KEGG_report:
                            break
                        if isPathway:
                            pathways.append(KEGG_report.lstrip().rstrip())
                    out.write(", ".join(pathways)+"\n")

rule get_genes_per_pubmed:
    input:
        "data/RNA-seq-gene-info.txt",
        "data/RNA-seq-ids.txt"
    output:
        "data/RNA-seq-genes-per-pubmed.txt"
    run:
        with open(input[1]) as id_file:
            with open(input[0]) as pubmed_file:

                pubmed_dict = {}

                for line in id_file:
                    if not line.strip():
                        continue
                    else:
                        gene_id = line.rstrip()
                        pubmed_ids = pubmed_file.readline().split("\t")[4].rstrip().split(",")

                        for pm_id in pubmed_ids:
                            if pm_id in pubmed_dict:
                                pubmed_dict[pm_id].append(gene_id)
                            if not pm_id in pubmed_dict:
                                pubmed_dict[pm_id] = [gene_id]

                with open(output[0], "w") as out:
                    for pm in pubmed_dict:
                        out.write(pm+"\t"+", ".join(pubmed_dict[pm])+"\n")

rule create_gc_graphs:
    input:
        "data/RNA-seq-sequences.txt",
        "data/RNA-seq-ids.txt"
    output:
        dynamic("data/figures/{gene}.png")
    script:
        "scripts/create_graph.R"

def get_figures(wildcards):
    path_list = []
    with open("data/RNA-seq-ids.txt") as f:
        for id in f:
            if not id.split():
                continue
            else:
                path_list.append("data/figures/"+id.rstrip()+".png")
    return path_list

rule report:
    input:
        "data/RNA-seq-ids.txt",
        "data/RNA-seq-ncbi-ids.txt",
        "data/RNA-seq-gene-info.txt",
        "data/RNA-seq-sequences.txt",
        "data/RNA-seq-kegg-ids.txt",
        "data/RNA-seq-genes-per-pubmed.txt",
        "data/RNA-seq-orthologs.txt",
        get_figures
    output:
        "report.html"
    script:
        "scripts/report_parser.py"
