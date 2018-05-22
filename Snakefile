mail="mymail@gmail.com"
rule filter_ids:
    input:
        "RNA-Seq-counts-1.txt"
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

rule get_gene_ids:
    input:
        "data/RNA-seq-ids.txt"
    output:
        "data/RNA-seq-ncbi-ids.txt"
    run:
        import sys
        from Bio import Entrez
        Entrez.email=mail
        with open(input[0]) as f:
            for id in f:
                handle = Entrez.esearch(db="gene", term=id,retmode="xml")
                record = Entrez.read(handle)
                with open(output[0], "a") as out:
                    out.write(str(record["IdList"][0])+"\n")

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
                    handle = Entrez.efetch(db="gene", id=id, retmode="xml")
                    record = Entrez.read(handle)
                    genome = record[0]['Entrezgene_gene-source']['Gene-source']['Gene-source_src-str1']
                    start = record[0]['Entrezgene_locus'][0]['Gene-commentary_seqs'][0]['Seq-loc_int']['Seq-interval']['Seq-interval_from']
                    stop = record[0]['Entrezgene_locus'][0]['Gene-commentary_seqs'][0]['Seq-loc_int']['Seq-interval']['Seq-interval_to']
                    refs = record[0]['Entrezgene_locus'][0]['Gene-commentary_products'][0]['Gene-commentary_refs']
                    pubmed = []
                    for pubmed_id in refs:
                        pubmed.append(pubmed_id['Pub_pmid']['PubMedId'])
                    out.write(genome+"\t"+start+"\t"+stop+"\t"+",".join(pubmed)+"\n")

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

                    out.write(str(GC(record[0]['TSeq_sequence']))+"\t"+record[0]['TSeq_sequence']+"\n")

rule get_kegg_ids:
    input:
        "data/RNA-seq-ids.txt"
    output:
        "data/RNA-seq-kegg-ids.txt"
    run:
        from Bio.KEGG.REST import kegg_get

        request = kegg_get("lpl:{id}")
