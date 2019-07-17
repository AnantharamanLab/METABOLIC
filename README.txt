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
	

2 METABOLIC requires the KofamKOALA hmm database
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
