# RNA-Seq-Workflow
A workflow for getting gene information based on RNA-Seq data
This workflow is build using Snakemake, a pythonic workflow system

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
     
  You have to open a new terminal inorder to use conda
  
  <br/>
  
  #### Creating the environment using the provided environment file
  
  `$ conda env create --name {your-environment-name} --file environment.yaml`
  
  To activate your environment:
     
  `$ source activate {your-environment-name}`
     
  To deactivate simply:
     
  `$ source deactivate`
  
  <br/>
   
  #### Now the workflow can be run by using the `snakemake` command
   
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
