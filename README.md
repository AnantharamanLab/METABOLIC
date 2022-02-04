<img src="https://github.com/AnantharamanLab/METABOLIC/blob/master/METABOLIC.jpg" width="85%">


# METABOLIC
**MET**abolic **A**nd **B**ioge**O**chemistry ana**L**yses **I**n mi**C**robes  
Current Version: 4.0
Tested on: Ubuntu 18.04.5 LTS (Linux 5.4.0-81-generic x86_64) (Sep 2021)

This software enables the prediction of metabolic and biogeochemical functional trait profiles to any given genome datasets. These genome datasets can either be metagenome-assembled genomes (MAGs), single-cell amplified genomes (SAGs) or isolated strain sequenced genomes. METABOLIC has two main implementations, which are METABOLIC-G and METABOLIC-C. METABOLIC-G.pl allows for generation of metabolic profiles and biogeochemical cycling diagrams of input genomes and does not require input of sequencing reads. METABOLIC-C.pl generates the same output as METABOLIC-G.pl, but as it allows for the input of metagenomic read data, it will generate information pertaining to community metabolism. It can also calculate the genome coverage. The information is parsed and diagrams for elemental/biogeochemical cycling pathways (currently Nitrogen, Carbon, Sulfur and "other") are produced.  

| Program Name|Program Description |
|---|---|
|METABOLIC-G.pl|Allows for classification of the metabolic capabilities of input genomes. |
|METABOLIC-C.pl|Allows for classification of the metabolic capabilities of input genomes, <br />calculation of genome coverage, creation of biogeochemical cycling diagrams,<br /> and visualization of community metabolic interactions and contribution to biogeochemical processes by each microbial group. |
<br>

Slides of introducing METABOLIC were provided here: (https://github.com/AnantharamanLab/METABOLIC/blob/master/METABOLIC_C-DEBI_slides.pdf) 

(A C-DEBI series meeting talk)



If you are using this program, please consider citing our preprint, available at [BioRxiv](https://www.biorxiv.org/content/10.1101/761643v2):
```
Zhichao Zhou, Patricia Q. Tran, Adam M. Breister, Yang Liu, Kristopher Kieft, Elise S. Cowley, Ulas Karaoz, Karthik Anantharaman."METABOLIC: High-throughput profiling of microbial genomes for functional traits, biogeochemistry, and community-scale metabolic networks"
bioRxiv 761643; doi: https://doi.org/10.1101/761643
```
<br>
<br>

## Installing and using METABOLIC

Please see the project home page for usage details and installation instructions: <br>
https://github.com/AnantharamanLab/METABOLIC/wiki<br>



