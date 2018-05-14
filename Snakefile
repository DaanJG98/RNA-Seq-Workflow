rule filter_ids:
    input:
        "RNA-Seq-counts-1.txt"
    output:
        "data/RNA-seq-ids.txt"
    run:
        with open(output[0],"w") as out:
            with open(input[0]) as f:
                for line in f:
                    if not line.startswith("#"):
                        out.write(line.split("\t")[0])

rule get_gene_ids:
    input:
        "data/RNA-seq-ids.txt"
    output:
        "data/RNA-seq-ncbi-ids.txt"
    run:
        import sys
        from Bio import Entrez
        Entrez.email="mymail@gmail.com"
        id_list = []
        with open(input[0]) as f:
            for id in f:
                if not id.startswith("ID"):
                    id_list.append(str(id).rstrip())
        #print(",".join(id_list))

        for item in id_list:
            handle = Entrez.esearch(db="nucleotide", term=item)
            record = Entrez.read(handle)
            print(record['IdList'])

            with open(output[0], "w") as out:
                out.write(str(record["IdList"])+"\n")
