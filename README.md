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

Slides of introducing METABOLIC (for a C-DEBI series meeting presentation) were provided here: (https://github.com/AnantharamanLab/METABOLIC/blob/master/METABOLIC_C-DEBI_slides.pdf) 

(The carbon fixation pathway automated annotation gets updated - in Appendix)
### GTDB-Tk and the GTDB database version (update Aug 2026, tested on Ubunutu 24.04 LTS)

METABOLIC-C.pl runs `gtdbtk classify_wf` to assign taxonomy. The GTDB-Tk
version in `METABOLIC_v4.0_env.yml` (gtdbtk 2.3.2) supports GTDB reference
data up to R214.

If you are using a newer GTDB release, run GTDB-Tk yourself and pass its
output folder to METABOLIC-C.pl, using the `-gtdbtk-dir` option.


<br>
<br>

If you are using this program, please consider citing our paper, available at [Microbiome](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-021-01213-8):
```
Zhou, Z., Tran, P.Q., Breister, A.M. et al. METABOLIC: high-throughput profiling of microbial genomes for functional traits, metabolism, biogeochemistry, and community-scale functional networks. Microbiome 10, 33 (2022). https://doi.org/10.1186/s40168-021-01213-8
```
<br>
<br>

## Installing and using METABOLIC

Please see the project home page for usage details and installation instructions: <br>
https://github.com/AnantharamanLab/METABOLIC/wiki<br>



