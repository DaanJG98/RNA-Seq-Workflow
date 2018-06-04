# RNA-Seq-Workflow
A workflow for getting information for genes based on RNA-seq data. This workflow was built using Snakemake, a pythonic workflow system, and uses Python, Shell and R scripts for its various tasks. The workflow gathers the following information for genes and returns it a human-readable report:
* NCBI ID
* Nucleotide sequence and GC-percentage
* Literature (PubMed)
* KEGG pathways
* Orthologs

## Get started

Snakemake can be run on both Linux and Windows machines, though we find it more convenient to work with a Linux machine.
If you are on a Windows platform, don't worry there are multiple solutions:
1. Setting up a Linux virtual machine with for example VirtualBox
2. Using Vagrant to house your virtual machine. We recommend this option if you don't need Linux for other purposes

In this readme we'll cover the setup for both a full Linux OS and a Vagrant Linux VM
  

### Setup using Vagrant
  #### First you need to download and install [Vagrant](https://www.vagrantup.com/downloads.html)
  
  * Download the right package for your operating system
     
  * Run the installer and check if vagrant installed correctly
        
     * `> vagrant` 
  <br/> 
   
  #### Next checkout this git repo in a conveniet directory
  
  `> git clone https://github.com/DaanJG98/RNA-Seq-Workflow.git`
  
  <br/>
  
  #### cd into this directory
  <br/>
  
  #### Create a 64-bit Ubuntu environment in this directory with:

  `> vagrant init hashicorp/precise64`

  `> vagrant up`

  To get into the virtual machine type:

  `> vagrant ssh`
 
  <br/>
 
  *Make sure you are into your Vagrant VM when performing the next steps*
  #### Installing Miniconda 3
  
  `$ wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh`
     
  `$ bash Miniconda3-latest-Linux-x86_64.sh`
     
  You have to open a new terminal inorder to use conda, make sure to logged in to your virtual environment.
  
  <br/>
  
  #### Creating the environment using the provided environment file
  
  `$ conda env create --name {your-environment-name} --file /vagrant/environment.yaml`
  
  To activate your environment:
     
  `$ source activate {your-environment-name}`
     
  To deactivate simply:
     
  `$ source deactivate`
  
  <br/>
   
  #### Now the workflow can be run by using the `snakemake` command
  
  `$ cd /vagrant` 
  
  `$ snakemake {optional parameters}`
    
  For more info about snakemake check the docs at: [snakemake docs](http://snakemake.readthedocs.io/en/stable/)
  
  <br/>
  <br/>

### Setup using a full Linux machine
  #### Start by checking out this git repo in a new directory
  
  `> git clone https://github.com/DaanJG98/RNA-Seq-Workflow.git`
  
  <br/>
  
  #### Next install Miniconda 3
  
  `$ wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh`
    
  `$ bash Miniconda3-latest-Linux-x86_64.sh`
  
  <br/>
  
  #### Create the environment with the provided environment file
  
  `$ conda env create --name {your-environment-name} --file {path-to-file}/environment.yaml`

  To activate your environment:

  `$ source activate {your-environment-name}`

  To deactivate simply:

  `$ source deactivate`
  
  For more info about snakemake check the docs at: [snakemake docs](http://snakemake.readthedocs.io/en/stable/)
  
### Data format
RNA-seq data should be provided in the following, tab separated, format:

| ID     | Condition 1 | Condition ...  |
| :----: | :---------: | :------------: |
| gene A | {value}     | ...            |
| gene B | {value}     | ...            |
| ...    | ...         | ...            |  
  
## Snakefile rules documentation

| conv_kegg_ids   |                                                                   |
| --------------- | :---------------------------------------------------------------- |
| **Input**       | RNA-seq-ncbi-ids.txt                                              |
| **Output**      | RNA-seq-conv-kegg.txt                                             |
| **Script**      | conv_kegg.sh                                                      |
| **Description** | Pass NCBI IDs to KEGG REST API and return corresponding KEGG IDs. |

| create_gc_graphs |                                                 |
| ---------------- | :---------------------------------------------- |
| **Input**        | RNA-seq-sequences.txt, RNA-seq-ids.txt          |
| **Output**       | {gene}.png                                      |
| **Script**       | create_graph.R                                  |
| **Description**  | Create graph showing GC% and AT% for each gene. |

| filter_ids      |                               |
| --------------- | :---------------------------- |
| **Input**       | RNA-Seq-counts.txt            |
| **Output**      | RNA-seq-ids.txt               |
| **Script**      | In-rule Python script         |
| **Description** | Filter IDs out of input file. |

| get_gene_info   |                                                                                             |
| --------------- | :------------------------------------------------------------------------------------------ |
| **Input**       | RNA-seq-ncbi-ids.txt                                                                        |
| **Output**      | RNA-seq-gene-info.txt                                                                       |
| **Script**      | In-rule Python script                                                                       |
| **Description** | Fetch gene data from Entrez, make selection of specific attributes and return these values. |

| get_genes_per_pubmed   |                                                 |
| ---------------------- | :---------------------------------------------- |
| **Input**              | RNA-seq-gene-info.txt, RNA-seq-ids.txt          |
| **Output**             | RNA-seq-ids.txt                                 |
| **Script**             | In-rule Python script                           |
| **Description**        | Get per PubMed article the corresponding genes. |

| get_kegg_ids    |                                                                       |
| --------------- | :-------------------------------------------------------------------- |
| **Input**       | RNA-seq-conv-kegg.txt                                                 |
| **Output**      | RNA-seq-kegg-ids.txt                                                  |
| **Script**      | In-rule Python script                                                 |
| **Description** | Pass KEGG IDs to Bio.KEGG REST and return corresponding pathway IDs.  |

| get_ncbi_ids    |                                                                  |
| --------------- | :--------------------------------------------------------------- |
| **Input**       | RNA-seq-ids.txt                                                  |
| **Output**      | RNA-seq-ncbi-ids.txt                                             |
| **Script**      | get_ncbi_ids.sh                                                  |
| **Description** | Pass IDs from input to Entrez and return corresponding NCBI IDs. |

| get_orthologs   |                                                                                   |
| --------------- | :-------------------------------------------------------------------------------- |
| **Input**       | RNA-seq-ids.txt                                                                   |
| **Output**      | RNA-seq-orthologs.txt                                                             |
| **Script**      | get_omadb_orthologs.sh                                                            |
| **Description** | Pass IDs from input to OMA Browser API and return IDs of corresponding orthologs. |

| get_sequence    |                                                                                                 |
| --------------- | :---------------------------------------------------------------------------------------------- |
| **Input**       | RNA-seq-gene-info.txt                                                                           |
| **Output**      | RNA-seq-sequences.txt                                                                           |
| **Script**      | In-rule Python script                                                                           |
| **Description** | Fetch specific sequence in genome from Entrez, calculate GC-percentage and return these values. |

| report          |                                                                 |
| --------------- | :-------------------------------------------------------------- |
| **Input**       | RNA-seq-ids.txt, RNA-seq-ncbi-ids.txt, RNA-seq-gene-info.txt, RNA-seq-sequences.txt, RNA-seq-kegg-ids.txt, RNA-seq-genes-per-pubmed.txt, RNA-seq-orthologs.txt                        |
| **Output**      | report.html                                                     |
| **Script**      | report_parser.py                                                |
| **Description** | Parse all output data from previous rules into a single report. |
