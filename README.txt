[Install instruction]
1 METABOLIC requires the following programs to be added to your system path:

  1.1 HMMER (>=v3.1b2) 
    Link: http://hmmer.org/
  1.2 prodigal (>=v2.6.3) 
    Link: https://github.com/hyattpd/Prodigal 
    Remarks: executable must be named prodigal and not prodigal.linux
  1.3 MetaBAT2 (>=v2.12.1)
    Link: https://bitbucket.org/berkeleylab/metabat

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
    cp ../../batch_hmmpress.pl ./
    perl batch_hmmpress.pl
    
   2.2 The METABOLIC hmm database in "METABOLIC_hmm_db.tgz" contains custom hmm files, self-parsed Pfam and TIRGfam files.
       It needs to be unpressed to the folder "METABOLIC_hmm_db" and stays in the same directory of KofamKOALA hmm database
       and scripts.

3 METABOLIC uses the "METABOLIC_temp_and_db" which contains the hmm result table and KEGG database information
  Unpress the METABOLIC_temp_and_db.tgz to the folder "METABOLIC_temp_and_db" and keep it in the same directory of 
  KofamKOALA hmm database and scripts.
        
