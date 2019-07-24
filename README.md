#[Description]
METABOLIC.pl

METABOLIC  =>  METabolic And BiogeOchemistry anaLyses In miCrobes

This software gives a metabolic and biogeochemical function trait profile to given genome datasets 
[either metagenome-assembled genomes (MAGs), single-cell amplified genomes (SAGs) or pure culture sequenced genomes]. 
It also integrates the genome coverage to make element cycling pathways.

Copyright:
Zhichao Zhou, zczhou2017@gmail.com
Karthik Anantharaman, karthik@bact.wisc.edu
Anantharaman Microbiome Laboratory
Department of Bacteriology, University of Wisconsin, Madison


#[Install instruction]
1 METABOLIC requires the following programs to be added to your system path:

  1.1 Perl (>= v5.010)
  1.2 HMMER (>= v3.1b2) 
    Link: http://hmmer.org/
  1.3 prodigal (>= v2.6.3) 
    Link: https://github.com/hyattpd/Prodigal 
    Remarks: executable must be named prodigal and not prodigal.linux
  1.4 SAMtools (>= v0.1.19)
    Link: http://www.htslib.org/
  1.5 BAMtools (>= v2.4.0)
    Link: https://github.com/pezmaster31/bamtools/wiki
  1.6 CoverM
    Link: https://github.com/wwood/CoverM

You could follow the install instruction of each program; or you could also
install them via Conda and add them to your system path:
  Link: https://anaconda.org
	

2 METABOLIC requires the KofamKOALA hmm and METABOLIC hmm database 
  Link: https://www.genome.jp/tools/kofamkoala/

  2.1 Download KofamKOALA hmm database files:
    mkdir kofam_database
    cd kofam_database
    wget -c ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz
    wget -c ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz
    gzip -d ko_list.gz
    tar xzf profiles.tar.gz; rm profiles.tar.gz
    cd profiles
    cp ../../Accessory_scripts/batch_hmmpress.pl ./
    perl batch_hmmpress.pl
    
   2.2 The METABOLIC hmm database in "METABOLIC_hmm_db.tgz" contains custom hmm files, self-parsed Pfam and TIRGfam files.
       It needs to be decompressed to the folder "METABOLIC_hmm_db" and stays in the same directory of KofamKOALA hmm database
       and scripts.

3 METABOLIC uses the "METABOLIC_temp_and_db" which contains the hmm result table and KEGG database information
  Decompress the METABOLIC_temp_and_db.tgz to the folder "METABOLIC_temp_and_db" and keep it in the same directory of 
  KofamKOALA hmm database and scripts.
  
4 This software also contains "Accessory_scripts.gz", which needs to be decompressed before use.

#[Result files]
METABOLIC result table
Each hmm hit protein collection
Element cycling pathways for each genome and a summary scheme (Both files and figures).
