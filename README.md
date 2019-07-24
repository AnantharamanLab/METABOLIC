<img src="https://github.com/ChaoLab/METABOLIC/blob/master/metabolic.png" width="100%">

[Description]
-------
**METABOLIC.pl**

**METABOLIC**  =>  **MET**abolic **A**nd **B**ioge**O**chemistry ana**L**yses **I**n mi**C**robes  

This software gives a metabolic and biogeochemical function trait profile to given genome datasets  
[either metagenome-assembled genomes (MAGs), single-cell amplified genomes (SAGs) or pure culture sequenced genomes].  
It also integrates the genome coverage to make element cycling pathways.  

Copyright:  
Zhichao Zhou, zczhou2017@gmail.com  
Patricia Tran, ptran5@wisc.edu  
Karthik Anantharaman, karthik@bact.wisc.edu  
Anantharaman Microbiome Laboratory  
Department of Bacteriology, University of Wisconsin, Madison  


[Install instruction]
-------
1 METABOLIC requires the following programs to be added to your system path:  

  1.1 Perl (>= v5.010)  
  1.2 [HMMER](http://hmmer.org/) (>= v3.1b2)   
  1.3 [Prodigal](https://github.com/hyattpd/Prodigal) (>= v2.6.3)   
      Remarks: executable must be named prodigal and not prodigal.linux  
  1.4 [SAMtools](http://www.htslib.org/) (>= v0.1.19)  
  1.5 [BAMtools](https://github.com/pezmaster31/bamtools/wiki) (>= v2.4.0)   
  1.6 [CoverM](https://github.com/wwood/CoverM)
  
You could follow the install instruction of each program; or you could also  
install them via Conda and add them to your system path:  
 Link: <https://anaconda.org>
  	  
  
2 METABOLIC requires the KofamKOALA hmm and METABOLIC hmm database   
  Link: <https://www.genome.jp/tools/kofamkoala/>   
  
  2.1 Download KofamKOALA hmm database files: 
```
    mkdir kofam_database  
    cd kofam_database  
    wget -c ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz  
    wget -c ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz  
    gzip -d ko_list.gz  
    tar xzf profiles.tar.gz; rm profiles.tar.gz  
    cd profiles  
    cp ../../Accessory_scripts/batch_hmmpress.pl ./  
    perl batch_hmmpress.pl  
 ```
      
   2.2 The METABOLIC hmm database in "METABOLIC_hmm_db.tgz" contains custom hmm files, self-parsed Pfam and TIRGfam files. It needs to be decompressed to the folder "METABOLIC_hmm_db" and stays in the same directory of KofamKOALA hmm database and scripts.  
  
3 METABOLIC uses the "METABOLIC_temp_and_db" which contains the hmm result table and KEGG database information.  
Decompress the METABOLIC_temp_and_db.tgz to the folder "METABOLIC_temp_and_db" and keep it in the same directory of KofamKOALA hmm database and scripts.     
4 This software also contains "Accessory_scripts.gz", which needs to be decompressed before use.  
  
[Result files]  
-------
METABOLIC result table  
Each hmm hit protein collection  
Element cycling pathways for each genome and a summary scheme (Both files and figures).  

[Sample inputs]  
-------
1 The genome files should be ended with ".fasta"; The genome protein files should be ended with ".faa";  
2 The "-o" option requires inputting a txt file to show the address of metagenomic reads [which is/are the metagenomic dataset(s) that you used to reconstructed the MAGs]. A sample for this txt is like:   
```
#Reads pair name with absolute pathway: 
/slowdata/Reads/METABOLIC_test_reads/SRR3577362_sub_1.fastq,/slowdata/Reads/METABOLIC_test_reads/SRR3577362_sub_2.fastq
/slowdata/Reads/METABOLIC_test_reads/SRR3577362_sub2_1.fastq,/slowdata/Reads/METABOLIC_test_reads/SRR3577362_sub2_2.fastq
```
The paired reads are in the same line. 

[Software work flow]  
-------
<img src="https://github.com/ChaoLab/METABOLIC/blob/master/Software_work_flow.jpg" width="75%">
