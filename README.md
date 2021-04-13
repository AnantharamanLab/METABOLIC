<img src="https://github.com/AnantharamanLab/METABOLIC/blob/master/METABOLIC.jpg" width="85%">


# METABOLIC
**MET**abolic **A**nd **B**ioge**O**chemistry ana**L**yses **I**n mi**C**robes  
Current Version: 4.0
Tested on: Linux Ubuntu 16.04.6 LTS (GNU/Linux 4.15.0-101-generic x86_64) (June 2020)

This software enables the prediction of metabolic and biogeochemical functional trait profiles to any given genome datasets. These genome datasets can either be metagenome-assembled genomes (MAGs), single-cell amplified genomes (SAGs) or pure culture sequenced genomes. METABOLIC has two main implementations, which are METABOLIC-G and METABOLIC-C. METABOLIC-G.pl allows for generation of metabolic profiles and biogeochemical cycling diagrams of input genomes and does not require input of sequencing reads. METABOLIC-C.pl generates the same output as METABOLIC-G.pl, but as it allows for the input of metagenomic read data, it will generate information pertaining to community metabolism. It can also calculate the genome coverage. The information is parsed and diagrams for elemental/biogeochemical cycling pathways (currently Nitrogen, Carbon, Sulfur and "other") are produced.  

| Program Name|Program Description |
|---|---|
|METABOLIC-G.pl|Allows for classification of the metabolic capabilities of input genomes. |
|METABOLIC-C.pl|Allows for classification of the metabolic capabilities of input genomes, <br />calculation of genome coverage, creation of biogeochemical cycling diagrams,<br /> and visualization of community metabolic interactions and energy flow. |

If you are using this program, please consider citing our preprint, available at [BioRxiv](https://www.biorxiv.org/content/10.1101/761643v2):
```
Zhichao Zhou, Patricia Q. Tran, Adam M. Breister, Yang Liu, Kristopher Kieft, Elise S. Cowley, Ulas Karaoz, Karthik Anantharaman."METABOLIC: High-throughput profiling of microbial genomes for functional traits, biogeochemistry, and community-scale metabolic networks"
bioRxiv 761643; doi: https://doi.org/10.1101/761643
```

## Table of Contents:
1. [Version History](#version_history)

2. [System Requirements](#system_requirements)

3. [Dependencies Overview](#dependencies_overview)

4. [Detailed Dependencies](#dependencies_detailed)

5. [Installation Instructions](#install_instructions)  
  &emsp;a. [Quick Installation](#quick_install)  

  &emsp;b. [Docker or Singularity container Installation](#docker_and_singularity_container_install)  

  &emsp;c. [Making conda environment](#making_conda_environment)  

6. [Running METABOLIC](#running_metabolic)  
  &emsp;a. [Required and Optional Flags](#flags)  
  &emsp;b. [How to Run](#running)  

7. [METABOLIC Output Descriptions](#metabolic_output)  
  &emsp;a. [Outputs Overview](#output_overview)  
  &emsp;b. [Outputs Detailed](#output_detailed)  

8. [Copyright](#copyright)

## <a name="version_history"></a> Version History:
v4.0 -- Jun 22, 2020 --  
* METABOLIC now uses an R script to generate METABOLIC_result.xlsx, which fixes issues with the generation of a corrupt METABOLIC_result.xlsx file
* Test input data now includes both five nucleotide fasta files and one set of paired sequencing reads, allowing all capabilities of both METABOLIC-G.pl and METABOLIC-C.pl to be tested
* The MN-score table has been provided as one of the results by METABOLIC-C
* Updated the motif checking step for pmo/amo, dsrE/tusD, dsrH/tusB, and dsrF/tusC
* Updated the "reads-type" option allowing the use of metatranscriptomic reads to conduct community analysis
* Updated the script error of assigning module step presence to be the same with module step presence

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

If you are planning to use only METABOLIC-G, you don't need to install GTDB-tk.

|Necessary Databases|Approximate System Storage Required|
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

Each of these programs should be in the PATH so that they can be accessed regardless of location.

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
&emsp;&emsp;&emsp;&emsp;3. *Getopt::Long*<br />
&emsp;&emsp;&emsp;&emsp;4. *Statistics::Descriptive*<br />
&emsp;&emsp;&emsp;&emsp;5. *Array::Split*<br />
&emsp;&emsp;&emsp;&emsp;6. *Bio::SeqIO*<br />
&emsp;&emsp;&emsp;&emsp;7. *Bio::Perl*<br />
&emsp;&emsp;&emsp;&emsp;8. *Bio::Tools::CodonTable*<br />
&emsp;&emsp;&emsp;&emsp;9. *Carp*<br />
&emsp;&emsp;&emsp;&emsp;10. *File::Spec*<br />
&emsp;&emsp;&emsp;&emsp;11. *File::Basename*<br />
&emsp;&emsp;&emsp;&emsp;12. *Parallel::ForkManager*<br />

Note for later: The last three lines in "run_to_setup.sh" are used to download "METABOLIC_test_files.tgz" from a  Google Drive. It requires gdown. gdown can be simply installed by calling "pip install gdown".

Note for later: The last three lines in "run_to_setup.sh" are used to download "METABOLIC_test_files.tgz" from a  Google Drive. It requires gdown. gdown can be simply installed by calling "pip install gdown".

**R Packages**  
&emsp;&emsp;To install, open the R command line interface by entering "R" into the command line, and then enter<br />
&emsp;&emsp;"install.packages("[Package Name]")".<br />

```
Example:
R
install.packages("diagram")
q()
```

&emsp;&emsp;&emsp;&emsp;1. *diagram* **(v1.6.4)**<br />
&emsp;&emsp;&emsp;&emsp;2. *forcats* **(v0.5.0)**<br />
&emsp;&emsp;&emsp;&emsp;3. *digest* **(v0.6.25)**<br />
&emsp;&emsp;&emsp;&emsp;4. *htmltools* **(v0.4.0)**<br />
&emsp;&emsp;&emsp;&emsp;5. *rmarkdown* **(v2.1)**<br />
&emsp;&emsp;&emsp;&emsp;6. *reprex* **(v0.3.0)**<br />
&emsp;&emsp;&emsp;&emsp;7. *tidyverse* **(v1.3.0)**<br />
&emsp;&emsp;&emsp;&emsp;8. *ggthemes* **(v4.2.0)**<br />
&emsp;&emsp;&emsp;&emsp;9. *ggalluvial* **(v0.11.3)**<br />
&emsp;&emsp;&emsp;&emsp;10. *reshape2* **(v1.4.3)**<br />
&emsp;&emsp;&emsp;&emsp;11. *ggraph* **(v2.0.2)**<br />
&emsp;&emsp;&emsp;&emsp;12. *pdftools* **(v2.3)**<br />
&emsp;&emsp;&emsp;&emsp;13. *igraph* **(v1.2.5)**<br />
&emsp;&emsp;&emsp;&emsp;15. *tidygraph* **(v1.1.2)**<br />
&emsp;&emsp;&emsp;&emsp;16. *stringr* **(v1.4.0)**<br />
&emsp;&emsp;&emsp;&emsp;17. *plyr* **(v1.8.6)**<br />
&emsp;&emsp;&emsp;&emsp;18. *dplyr* **(v0.8.5)**<br />
&emsp;&emsp;&emsp;&emsp;19. *openxlsx* **(v4.1.4)**<br />

**To ensure efficient and successful installation of METABOLIC, make sure that all dependencies are properly installed prior to download of the METABOLIC software.**

## <a name="install_instructions"></a> Installation Instructions:

1. Go to where you want the program to be and clone the github repository by using the following command:  
```
git clone https://github.com/AnantharamanLab/METABOLIC.git
```
&emsp;&emsp;or click the green button "download ZIP" folder at the top of the github and unzip the downloaded file.   
&emsp;&emsp;The perl and R scripts and dependent databases should be kept in the same directory.  

*NOTE: Before following the next step, make sure your working directory is the directory that was created by the METABOLIC download, that is, the directory containing the main scripts for METABOLIC (METABOLIC-G.pl, METABOLIC-C.pl, etc.).*

#### <a name="quick_install"></a> Quick Installation:

We provide a "run_to_setup.sh" script along with the data downloaded from the GitHub for easy setup of dependent databases. This can be run by using the following command:  

```
sh run_to_setup.sh
```
**Notice**: The last three lines in "run_to_setup.sh" were used to download "METABOLIC_test_files.tgz" from google drive. It requires gdown. gdown can be simply installed by calling "pip install gdown". Please also refer to **[gdown](https://github.com/wkentaro/gdown)**.  

#### <a name="docker_and_singularity_container_install"></a> Docker or Singularity container Installation:

Mr. Tin Ho has initiated a project to help intalling METABOLIC in Docker or Singularity container. Here is the link **[https://github.com/tin6150/metabolic](https://github.com/tin6150/metabolic)**. 

**Notice**:  There can be discrepancies between his version of METABOLIC and our version due to upgrading from time to time.

#### <a name="making_conda_environment"></a> Making conda environment:

Dr. Daan Speth has opened an issue teaching how to successfully install conda environment for METABOLIC-G/-C. Here is the link **[https://github.com/AnantharamanLab/METABOLIC/issues/27](https://github.com/AnantharamanLab/METABOLIC/issues/27)**. 

**Notice**:  According to Dr. Michal Strejcek, conda has issues with Perl; and if you are using R 4.x, some small changes need to be done in the script. Here is the link to his solutions to these issues **[https://github.com/AnantharamanLab/METABOLIC/issues/41](https://github.com/AnantharamanLab/METABOLIC/issues/41)**.

## <a name="running_metabolic"></a> Running METABOLIC: 

<img src="https://github.com/AnantharamanLab/METABOLIC/blob/master/METABOLICWorkflow.jpg" width="90%">

#### <a name="flags"></a> All Required and Optional Flags:

To view the options that METABOLIC-C.pl and METABOLIC-G.pl have, please type:

```
perl METABOLIC-G.pl -help
perl METABOLIC-C.pl -help
```

* **-in-gn [required if you are starting from nucleotide fasta files]** Defines the location of the FOLDER containing the genome nucleotide fasta files ending with ".fasta" to be run by this program
* **-in [required if you are starting from faa files]** Defines the location of the FOLDER containing the genome amino acid files ending with ".faa" to be run by this program
* **-r [required]** Defines the path to a text file containing the location of paried reads
* **-rt [optional]** Defines the option to use "metaG" or "metaT" to indicate whether you use the metagenomic reads or metatranscriptomic reads (default: 'metaG')
* **-t [optional]** Defines the number of threads to run the program with (Default: 20)
* **-m-cutoff [optional]** Defines the fraction of KEGG module steps present to designate a KEGG module as present (Default: 0.75)
* **-kofam-db [optional]** Defines the use of the full ("full") or reduced ("small") KOfam database by the program (Default: 'full')
* **-p [optional]** Defines the prodigal method used to annotate ORFs ("meta" or "single")
* **-o [optional]** Defines the output directory to be created by the program (Default: current directory)

1. The directory specified by the "-in-gn" flag should contain nucleotide sequences for your genomes with the file extension ".fasta". If you are supplying amino acid fasta files for each genome, these should be contained within a directory and have the file extension ".faa", and you will be using the "-in" option instead. Ensure that the fasta headers of each .fasta or .faa file is unique, and that your file names do not contains spaces. If you want to use METABOLIC-C, only "fasta" files are allowed to perform the good analysis.

2. The "-r" flag allows input of a text file defining the path of metagenomic reads (if running METABOLIC-C). The metagenomic reads refer to the metagenomic read datasets that you used to generate the MAGs. Try to confirm that you are using unzipped fastq files instead of zipped files before you run METABOLIC-C. Sets of paired reads are entered in one line, separated by a ",". A sample for this text file is as follows:     
```
#Read pairs: 
SRR3577362_sub_1.fastq,SRR3577362_sub_2.fastq
SRR3577362_sub2_1.fastq,SRR3577362_sub2_2.fastq
```
&emsp;&emsp;*Note that the two different sets of paired reads are separated by a line return (new line), and please avoid empty lines in this text file otherwise software will take blank read files as inputs.*

#### Running Test Data:

The main METABOLIC directory also contains a set of 5 genomes and one set of paired metagenomic reads, which can be used to test that METABOLIC-G and METABOLIC-C were installed correctly. These genomes and reads can be found within the directory `METABOLIC_test_files/`, which is contained within the METABOLIC program directory.

METABOLIC-C.pl and METABOLIC-G.pl can be run with the test data by using the `-test true` function of METABOLIC:
```
perl METABOLIC-G.pl -test true

perl METABOLIC-C.pl -test true
```

#### <a name="running"></a> How To Run METABOLIC:

The main scripts that should be used to run the program are METABOLIC-G.pl or METABOLIC-C.pl.  

In order to run METABOLIC-G starting from **nucleotide** sequences, **AT LEAST** the following flags should be used for METABOLIC-G:
```
perl METABOLIC-G.pl -in-gn [path_to_folder_with_genome_files] -o [output_directory_to_be_created]
```

In order to run METABOLIC-G starting from **amino acid** sequences, **AT LEAST** the following flags should be used for METABOLIC-G:
```
perl METABOLIC-G.pl -in [path_to_folder_with_genome_files] -o [output_directory_to_be_created]
```

In order to run METABOLIC-C, **AT LEAST** the following flags should be used for METABOLIC-C:
```
perl METABOLIC-C.pl -in-gn [path_to_folder_with_genome_files] -r [path_to_list_of_paired_reads] -o [output_directory_to_be_created]
```

<br />


## <a name="metabolic_output"></a> METABOLIC Output Files:

#### <a name="output_overview"></a> Output Files Overview:

|  Output File  |   File Description   |  Generated by METABOLIC-C |   Generated by METABOLIC-G   |
|-------|---------|:---------:|:---------:|
|All_gene_collections_mapped.depth.txt|The gene depth of all input genes |X||
|[Each_HMM_Amino_Acid_Sequence/](#faa_collection)|The faa collection for each hmm file |X|X|
|intermediate_files/|The hmmsearch, peptides (MEROPS), and CAZymes (dbCAN2) running intermediate files |X|X|
|[KEGG_identifier_result/](#kegg_ident_result)|The hit and result of each genome by Kofam database |X|X|
|[METABOLIC_Figures/](#metabolic_figures_out)|All figures output from the running of METABOLIC |X|X|
|METABOLIC_Figures_Input/|All input files for R-generated diagrams |X|X|
|METABOLIC_result_each_spreadsheet/|TSV files representing each sheet of the created METABOLIC_result.xlsx file |X|X|
|MN-score_result/|The resulted table for MN-score |X||
|[METABOLIC_result.xlsx](#metabolic_result_table)|The resulting excel file of METABOLIC |X|X|
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


- **<a name="faa_collection"></a> Each HMM Profile Hit Amino Acid Sequence Collection (`Each_HMM_Amino_Acid_Sequence/`)** 

A collection of all amino acid sequences extracted from the input genome .faa files that were identified as matches to the custom HMM profiles provided by METABOLIC.

- **<a name="kegg_ident_result"></a> KEGG identifier results (`KEGG_identifier_result/`)**  

The KEGG identifier searching result - KEGG identifier numbers and hits of each genome that could be used to visualize the pathways in KEGG Mapper

- **<a name="metabolic_figures_out"></a> All figures generated by METABOLIC-G.pl and METABOLIC-C.pl (`METABOLIC_Figures/`)**

Both METABOLIC-G.pl and METABOLIC-C.pl will generate a folder titled `Nutrient_Cycling_Diagrams/` within the `METABOLIC_Figures/` directory, which will contain figures that represent nutrient cycling pathways for Sulfur, Nitrogen, Carbon, and other select pathways found within each genome. METABOLIC-C.pl also has the ability to generate overall community nutrient cycling pathways.

*Although the `Nutrient_Cycling_Diagrams/` directory is generated by both METABOLIC-G.pl and METABOLIC-C.pl, the files contained within the directory will be dependent on which script is used.*

For both programs, METABOLIC-G.pl and METABOLIC-C.pl, the `Nutrient_Cycling_Diagrams/` directory will contain the following files:
```
  [GenomeName].draw_sulfur_cycle_single.PDF
  [GenomeName].draw_nitrogen_cycle_single.PDF
  [GenomeName].draw_other_cycle_single.PDF
  [GenomeName].draw_carbon_cycle_single.PDF
```
*A red arrow designates presence of a pathway step and a black arrow means absence. Note the the width of the arrows does not have any significance.*

If you run METABOLIC-C.pl, the software will also calculate relative gene abundances, which will allow for generation of summary diagrams for pathways at a community scale:
```
  draw_sulfur_cycle_total.PDF
  draw_other_cycle_total.PDF
  draw_nitrogen_cycle_total.PDF
  draw_carbon_cycle_total.PDF
```
*Note the the width of the arrows does not have any significance.*

&emsp;&emsp;> Generated only by METABOLIC-C.pl are a set of figures representing metabolic handoffs within the community:  

For Sequential transformation diagram, we have summarized and visualized the genome number and genome coverage (relative abundance of microorganism) of the microorganisms that were putatively involved in the sequential transformation of both important inorganic elements and organic compounds. 

The resulting files are `Sequential_transformation_01.pdf` and `Sequential_transformation_02.pdf`.        

&emsp;&emsp;> Generated only METABOLIC-C.pl is a figure reprsenting energy flow by the community:  

For Metabolic energy flow diagram, a Sankey diagram is generated, representing the function fractions that are contributed by various microbial groups in a given community.   

The resulting file is `Metabolic_energy_flow.pdf`.

&emsp;&emsp;> METABOLIC-C.pl generates a figure reprsenting metabolic connections between different reactions that are  
&emsp;&emsp;  found within the community:  

For Metabolic network diagrams, diagrams representing metabolic connections of biogeochemical cycling steps at both phylum level and the whole community level will be generated.    

The resulted files are placed in the directory `Metabolic_network/`.

For MN-score result, the table showing the MN-score (Metabolic Networking score) will be generated ("MN-score_result.txt"). The first column indicates the MN-score for each function. The rest part of the table indicates the contribution percentage of each phylum to the corresponding function. An example was given: 

<img src="https://github.com/AnantharamanLab/METABOLIC/blob/master/MN-score_table_example.jpg" width="75%">

The resulted files are placed in the directory `MN-score result/`.



***Notice:***

If you use metatranscriptomic reads instead of metagenomic reads in METABOLIC-C, gene coverage result will be replaced by transcript coverage [normalized into Reads Per Kilobase of transcript, per Million mapped reads (RPKM)] and all the community analyses were performed based on the transcript coverage instead.  A result file of "All_gene_collections_transcript_coverage.txt" will be generated in the output directory in addition.


## <a name="copyright"></a> Copyright
METABOLIC: **MET**abolic **A**nd **B**ioge**O**chemistry ana**L**yses **I**n mi**C**robes (C) 2019    
Zhichao Zhou, zczhou2017@gmail.com  
Patricia Tran, ptran5@wisc.edu  
Karthik Anantharaman, karthik@bact.wisc.edu  
Anantharaman Microbiome Laboratory  
Department of Bacteriology, University of Wisconsin, Madison 

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see https://www.gnu.org/licenses/.
