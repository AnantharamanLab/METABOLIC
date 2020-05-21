<img src="https://github.com/AnantharamanLab/METABOLIC/blob/master/METABOLIC.jpg" width="85%">


# METABOLIC
**MET**abolic **A**nd **B**ioge**O**chemistry ana**L**yses **I**n mi**C**robes  
Current Version: 3.0

This software enables the prediction of metabolic and biogeochemical functional trait profiles to any given genome datasets. These genome datasets can either be metagenome-assembled genomes (MAGs), single-cell amplified genomes (SAGs) or pure culture sequenced genomes. METABOLIC has two main implementations, which are METABOLIC-G and METABOLIC-C. METABOLIC-G.pl allows for generation of metabolic profiles and biogeochemical cycling diagrams of input genomes and does not require input of sequencing reads. METABOLIC-C.pl generates the same output as METABOLIC-G.pl, but as it allows for the input of metagenomic read data, it will generate information pertaining to community metabolism. It can also calculate the genome coverage. The information is parsed and diagrams for elemental/biogeochemical cycling pathways (currently Nitrogen, Carbon, Sulfur and "other") are produced.  

| Program Name|Program Description |
|---|---|
|METABOLIC-G.pl|Allows for classification of the metabolic capabilities of input genomes. |
|METABOLIC-C.pl|Allows for classification of the metabolic capabilities of input genomes, <br />calculation of genome coverage, creation of biogeochemical cycling diagrams,<br /> and visualization of community metabolic interactions and energy flow. |

If you are using this program, please consider citing our preprint, available on [BioRxiv](https://www.biorxiv.org/content/10.1101/761643v1):
```
Zhou Z, Tran P, Liu Y, Kieft K, Anantharaman K. "METABOLIC: A scalable high-throughput metabolic and biogeochemical functional trait profiler based on microbial genomes" (2019). BioRxiv doi: https://doi.org/10.1101/761643
```

## Table of Contents:
1. [Version History](#version_history)
2. [System Requirements](#system_requirements)
3. [Dependencies Overview](#dependencies_overview)
4. [Detailed Dependencies](#dependencies_detailed)
5. [Installation Instructions](#install_instructions)  
&emsp;a. [Quick Installation](#quick_install)  
6. [Running METABOLIC](#running_metabolic)  
&emsp;a. [Required and Optional Flags](#flags)  
&emsp;b. [How to Run](#running)  
7. [METABOLIC Output Descriptions](#metabolic_output)  
&emsp;a. [Outputs Overview](#output_overview)  
&emsp;b. [Outputs Detailed](#output_detailed)  
8. [Copyright](#copyright)

## <a name="version_history"></a> Version History:
v3.0 -- Feb 18, 2020 --     
* Provide an option to let the user reduce the size of Kofam Hmm profiles (only use KOs that can be found in Modules) to speed up the calculation    
* Change HMMER to v3.3 to speed up the calculation<br />

v2.0 -- Nov 5, 2019 --     
* Add more functions on visualization, add more annotations, make the software faster<br /> 

v1.3 -- Sep 5, 2019 --     
* Fix the output folder problem, the perl script could be called in another place instead of the original place<br />

v1.2 -- Sep 5, 2019 --    
* Fix the prodigal parallel run, change "working-dir" to "METABOLIC-dir"<br />

v1.1 -- Sep 4, 2019 --     
* Fix the parallel problem, change from hmmscan to hmmsearch, and update the "METABOLIC_template_and_database"<br />   

## <a name="system_requirements"></a> System Requirements:

**System Memory Requirements:**  
- Due to requirements of some of this program's dependencies, it is highly recommended that METABOLIC-C is run on a system containing at least 100 Gb of memory.
- METABOLIC-G is not as demanding as METABOLIC-C and requires significantly less memory to run.

**System Storage Requirements:**

If you are planning to use only METABOLIC-g, you don't need to install GTDB-tk.

|Necessary Databases|Approximaate System Storage Required|
|---|---|
|METABOLIC program with unzipped files|7.69 Gb|
|GTDB-Tk Reference Data|28 Gb|

## <a name="dependencies_overview"></a> Dependencies Overview:
1. **[Perl](https://www.perl.org) (>= v5.010)**
2. **[HMMER](http://hmmer.org/) (>= v3.1b2)**
3. **[Prodigal](https://github.com/hyattpd/Prodigal) (>= v2.6.3)**
4. **[Sambamba](https://github.com/biod/sambamba) (>= v0.7.0)** (only for METABOLIG-C)
5. **[BAMtools](http://github.org/pezmaster31/bamtools) (>= v2.4.0)** (only for METABOLIG-C)
6. **[CoverM](https://github.com/wwood/CoverM)** (only for METABOLIG-C)
7. **[R](https://www.r-project.org/) (>= 3.6.0)**
8. **[Diamond](https://github.com/bbuchfink/diamond)**
9. **[Samtools](https://www.htslib.org)** (only for METABOLIG-C)
10. **[Bowtie2](https://github.com/BenLangmead/bowtie2)** (only for METABOLIG-C)
11. **[Gtdb-Tk](https://ecogenomics.github.io/GTDBTk/)** (only for METABOLIG-C)

## <a name="dependencies_detailed"></a> Perl and R Dependencies Detailed Instructions:

**Perl Modules:**<br />
&emsp;&emsp;To install, use the cpan shell by entering "perl -MCPAN -e shell cpan" and then entering<br />
&emsp;&emsp;"install [Module Name]", or install by using "cpan -i [Module Name]", or by entering<br />
&emsp;&emsp;"cpanm [Module Name]".<br />
```
Example 1:
perl -MCPAN -e shell cpan
install Data::Dumper

Example 2:
cpan -i Data::Dumper

Example 3:
cpanm Data::Dumper
```
&emsp;&emsp;&emsp;&emsp;1. *Data::Dumper*<br />
&emsp;&emsp;&emsp;&emsp;2. *POSIX*<br />
&emsp;&emsp;&emsp;&emsp;3. *Excel::Writer::XLSX*<br />
&emsp;&emsp;&emsp;&emsp;4. *Getopt::Long*<br />
&emsp;&emsp;&emsp;&emsp;5. *Statistics::Descriptive*<br />
&emsp;&emsp;&emsp;&emsp;6. *Array::Split*<br />
&emsp;&emsp;&emsp;&emsp;7. *Bio::SeqIO*<br />
&emsp;&emsp;&emsp;&emsp;8. *Bio::Perl*<br />
&emsp;&emsp;&emsp;&emsp;9. *Bio::Tools::CodonTable*<br />
&emsp;&emsp;&emsp;&emsp;10. *Carp*<br />

**R Packages**  
&emsp;&emsp;To install, open the R command line interface by entering "R" into the command line, and then enter<br />
&emsp;&emsp;"install.packages("[Package Name]")".<br />

Type:

```
R
install.packages("diagram")
install.packages("forcats")
install.packages("digest")
install.packages("htmltools")
install.packages("rmarkdown")
install.packages("reprex")
install.packages("tidyverse")
install.packages("ggthemes")
install.packages("ggalluvial")
install.packages("reshape2")
install.packages("ggraph")
q()
```

**To ensure efficient and successful installation of METABOLIC, make sure that all dependencies are properly installed prior to download of the METABOLIC software.**

## <a name="install_instructions"></a> Installation Instructions:

1. Go to where you want the program to be and clone the github repository by using the following command:  
```
git clone https://github.com/AnantharamanLab/METABOLIC.git
```    
&emsp;&emsp;or click the green button "download ZIP" folder at the top of the github and unzip the downloaded file.   
&emsp;&emsp;The perl and R scripts and dependent databases should be kept in the same directory.  

*NOTE: Before following the next steps, make sure your working directory is the directory that was created by the METABOLIC download, that is, the directory containing the main scripts for METABOLIC (METABOLIC-G.pl, METABOLIC-C.pl, etc.).*

#### <a name="quick_install"></a> Quick Installation:

For steps 2-8, we provide a "run_to_setup.sh" script along with the data downloaded from the GitHub for easy setup of dependent databases. This can be run by using the following command:
```
sh run_to_setup.sh
```


## <a name="running_metabolic"></a> Running METABOLIC: 

<img src="https://github.com/AnantharamanLab/METABOLIC/blob/master/METABOLICWorkflow.png" width="75%">

#### <a name="flags"></a> All Required and Optional Flags:

To view the options that METABOLIC-c.pl and METABOLIC-g.pl have type:

```
perl METABOLIC-G.pl -help
perl METABOLIC-C.pl -help
```

* **-in-gn [required if you are starting from nucleotide fasta files]** Defines the location of the genome fasta files ending with ".fasta" to be run by this program
* **-in [required if you are starting from faa files]** Defines the location of the genome amino acid files ending with ".faa" to be run by this program
* **-m [required]** Defines the location of the all dependent databases and scripts, which are downloaded with the program or installed above (Default './')
* **-r [required]** Defines the path to a text file containing the location of paried reads
* **-t [optional]** Defines the number of threads to run the program with (Default: 20)
* **-m-cutoff [optional]** Defines the fraction of KEGG module steps present to designate a KEGG module as present (Default: 0.75)
* **-kofam-db [optional]** Defines the use of the full ("full") or reduced ("small") KOfam database by the program (Default: 'full')
* **-p [optional]** Defines the prodigal method used to annotate ORFs ("meta" or "single")
* **-o [optional]** Defines the output directory to be created by the program (Default: current directory)

1. The directory specified by the "-in-gn" flag should contain nucleotide sequences for your genomes with the file extension ".fasta". If you are supplying amino acid fasta files for each genome, these should be contained within a directory and have the file extension ".faa", and you will be using the "-in" option instead. Ensure that the fasta headers of each .fasta or .faa file is unique, and that your file names do not contains spaces.

2. The "-r" flag allows input of a text file defining the path of metagenomic reads (if running METABOLIC-C). The metagenomic reads refer to the metagenomic read datasets that you used to generate the MAGs. Sets of paired reads are entered in one line, separated by a ",". A sample for this text file is as follows:     
```
#Read pairs: 
SRR3577362_sub_1.fastq,SRR3577362_sub_2.fastq
SRR3577362_sub2_1.fastq,SRR3577362_sub2_2.fastq
```
&emsp;&emsp;*Note that the two different sets of paired reads are separated by a line return (new line)*

#### Running Test Data:

The main METABOLIC directory also contains a set of 5 genomes, which can be used to test that the program was installed correctly. Currently, test runs can only be done with METABOLIC-G, because the only supplied test files of genomes in fasta format. These genomes can be found within the directory `5_genomes_test/`, which is contained within the METABOLIC program directory.

#### <a name="running"></a> How To Run METABOLIC:

The main scripts that should be used to run the program are METABOLIC-G.pl or METABOLIC-C.pl.  

In order to run METABOLIC-G starting from **nucleotide** sequences, **AT LEAST** the following flags should be used for METABOLIC-G:
```
METABOLIC-G.pl -in-gn [path_to_genome_files] -m [path_to_METABOLIC_directory] -o [output_directory_to_be_created]
```

In order to run METABOLIC-G starting from **amino acid** sequences, **AT LEAST** the following flags should be used for METABOLIC-G:
```
METABOLIC-G.pl -in [path_to_genome_files] -m [path_to_METABOLIC_directory] -o [output_directory_to_be_created]
```

In order to run METABOLIC-C, **AT LEAST** the following flags should be used for METABOLIC-C:
```
METABOLIC-C.pl -in-gn [path_to_genome_files] -r [path_to_list_of_paired_reads] -m [path_to_METABOLIC_directory] -o [output_directory_to_be_created]
```

<br />

We also offer a supplementary script titled "METABOLIC-C-jump.pl", which allows for the running of METABOLIC-C.pl with a reduced run time. You can use METABOLIC-C-jump.pl if you have **previously run** METABOLIC-C.pl and the `intermediate_files/` directory (See "METABOLIC Output Files" below) has already been generated.

<br />
To run this script, make sure that the `intermediate_files/` directory is within a directory which will serve as the output directory [-o option] for METABOLIC-C-jump.pl. To run this script, **AT LEAST** the following flags should be used:

```
METABOLIC-C-jump.pl -in-gn [path_to_genome_files] -r [path_to_list_of_paired_reads] -m [path_to_METABOLIC_directory] -o [premade_output_directory]
```

&emsp;*Note: METABOLIC-C-jump.pl can be used to re-run METABOLIC only if errors arose **AFTER** generation of all intermediate files*

## <a name="metabolic_output"></a> METABOLIC Output Files:

#### <a name="output_overview"></a> Output Files Overview:

|  Output File  |   File Description   |  Generated by METABOLIC-C |   Generated by METABOLIC-G   |
|-------|---------|:---------:|:---------:|
|All_gene_collections.gene|All genes combined as the mapping reference|X||
|All_gene_collections.gene.scaffold.1.bt2|Mapping reference bowtie2 files |X||
|All_gene_collections.gene.scaffold.2.bt2|Mapping reference bowtie2 files |X||
|All_gene_collections.gene.scaffold.3.bt2|Mapping reference bowtie2 files |X||
|All_gene_collections.gene.scaffold.4.bt2|Mapping reference bowtie2 files |X||
|All_gene_collections.gene.scaffold.rev.1.bt2|Mapping reference bowtie2 files |X||
|All_gene_collections.gene.scaffold.rev.2.bt2|Mapping reference bowtie2 files |X||
|All_gene_collections_mapped.1.sorted.bam|Mapping result - sorted bam file |X||
|All_gene_collections_mapped.1.sorted.bam.bai|Mapping result - sorted bam index file |X||
|All_gene_collections_mapped.1.sorted.stat|Mapping result - mapping statistics |X||
|All_gene_collections_mapped.depth.txt|The gene depth of all input genes |X||
|[Each_hmm_faa/](#faa_collection)|The faa collection for each hmm file |X|X|
|intermediate_files/|The hmmsearch, peptides (MEROPS), and CAZymes (dbCAN2) running intermediate files |X|X|
|[KEGG_identifier_result/](#kegg_ident_result)|The hit and result of each genome by Kofam database |X|X|
|Metabolic_energy_flow_input.txt|The input file for metabolic energy flow network |X||
|[Metabolic_energy_flow.pdf](#energy_flow)|The resulted diagram for metabolic energy flow network |X||
|[Metabolic_network/](#metabolic_networks)|Metabolic network result (metabolic connections, including those for each taxon and for the whole community) |X||
|Metabolic_network_input.txt|The input file for drawing metabolic network diagrams |X||
|[METABOLIC_result.xlsx](#metabolic_result_table)|The resulted excel file of METABOLIC |X|X|
|R_hm_input_1.txt|The input file 1 of sequential transformation result |X||
|R_hm_input_2.txt|The input file 2 of sequential transformation result |X||
|R_input/|The input file of element cycling pathways (for each genome and the whole community) |X|X|
|[R_output/](#r_output)|The output diagrams of element cycling pathways (for each genome and the whole community) |X|X|
|[Sequential_transformation_01.pdf](#seq_diagram)|The output diagram 1 of sequential transformation result |X||
|[Sequential_transformation_02.pdf](#seq_diagram)|The output diagram 2 of sequential transformation result |X||
<br />

#### <a name="output_detailed"></a> Output Files Detailed:
<br />

- **<a name="metabolic_result_table"></a> METABOLIC result table (`METABOLIC_result.xlsx`)**

This spreadsheet has 6 sheets: 

1. "**HMMHitNum**" = Presence or absence of custom HMM profiles within each genome, the number of times the HMM profile was identified within a genome, and the scaffold on which the HMM profile was found. The sheet provides a presence/absence indicator, the number of times a protein was identified for a given genome, and the ORF(s) that represent the identified protein.
2. "**FunctionHit**" = Presence or absence of sets of proteins which were identified and displayed as separate proteins in the sheet titled "HMMHitNum". For each genome, the functions are identified as "Present" or "Absence".
3. "**KEGGModuleHit**" = Annotation of each genome with modules from the KEGG database organized by metabolic category. For each genome, the functions are identified as "Present" or "Absence".
4. "**KEGGModuleStepHit**" = Presence or absence of modules from the KEGG database within each genome separated into the steps that make up the module. For each genome, the functions are identified as "Present" or "Absence".
5. "**dbCAN2Hit**" = The dbCAN2 annotation results against all genomes (CAZy numbers and hits). For each genome, there are two distinct columns, which show the number of times a CAZy was identified and what ORF(s) represent the protein.
6. "**MEROPSHit**" = The MEROPS peptidase searching result (MEROPS peptidase numbers and hits). For each genome, there are two distinct columns, which show the number of times a peptidase was identified and what ORF(s) represent the protein.

*In all cases if you scroll down you will see what "Gn00X" colnames refer to (they are based on your fasta file names for the genomes you gave.* 


- **<a name="faa_collection"></a> Each HMM Profile Hit Amino Acid Sequence Collection (`Each_hmm_faa/`)** 

A collection of all amino acid sequences extracted from the input genome .faa files that were identified as matches to the custom HMM profiles provided by METABOLIC.

- **<a name="kegg_ident_result"></a> KEGG identifier results (`KEGG_identifier_result/`)**  

The KEGG identifier searching result - KEGG identifier numbers and hits of each genome that could be used to visualize the pathways in KEGG Mapper

- **<a name="r_output"></a> Elemental/Biogeochemical cycling pathways for each genome and a summary scheme (`R_output/`)**

*Although the `R_output/` directory is generated by both METABOLIC-G.pl and METABOLIC-C.pl, the files contained within the directory will be dependent on which script is used.*

For both programs, METABOLIC-G.pl and METABOLIC-C.pl, the `R_output/` directory will contain the following files:
```
  [GenomeName].draw_sulfur_cycle_single.PDF
  [GenomeName].draw_nitrogen_cycle_single.PDF
  [GenomeName].draw_other_cycle_single.PDF
  [GenomeName].draw_carbon_cycle_single.PDF
```
A red arrow designates presence of a pathway step and a black arrow means absence. Note the the with of the arrows does not have any significance.

If you run METABOLIC-C.pl, the software will also calculate relative gene abundances, which will allow for generation of summary diagrams for pathways at a community scale:
```
  draw_sulfur_cycle_total.PDF
  draw_other_cycle_total.PDF
  draw_nitrogen_cycle_total.PDF
  draw_carbon_cycle_total.PDF
```
Note the the with of the arrows does not have any significance.

- **<a name="seq_diagram"></a> Sequential transformation Diagrams (`Sequential_transformation_01.pdf`, `Sequential_transformation_02.pdf`)**   

For Sequential transformation diagram, we have summarized and visualized the genome number and genome coverage (relative abundance of microorganism) of the microorganisms that were putatively involved in the sequential transformation of both important inorganic elements and organic compounds. 

The resulting files are `Sequential_transformation_01.pdf` and `Sequential_transformation_02.pdf`.        

- **<a name="energy_flow"></a> Metabolic Energy Flow (`Metabolic_energy_flow.pdf`)**

For Metabolic energy flow diagram, a Sankey diagram is generated, representing the function fractions that are contributed by various microbial groups in a given community.   

The resulting file is `Metabolic_energy_flow.pdf`.

- **<a name="metabolic_networks"></a> Metabolic Network Diagrams (`Metabolic_network/`)**

For Metabolic network diagrams, diagrams representing metabolic connections of biogeochemical cycling steps at both phylum level and the whole community level will be generated.    

The resulted files are placed in the directory `Metabolic_network/`.


## <a name="copyright"></a> Copyright
Zhichao Zhou, zczhou2017@gmail.com  
Patricia Tran, ptran5@wisc.edu  
Karthik Anantharaman, karthik@bact.wisc.edu  
Anantharaman Microbiome Laboratory  
Department of Bacteriology, University of Wisconsin, Madison 
