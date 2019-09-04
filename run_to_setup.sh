tar zxvf Accessory_scripts.tgz
tar zxvf METABOLIC_hmm_db.tgz
tar zxvf METABOLIC_temp_and_db.tgz
tar zxvf Motif.tgz.tgz
mkdir kofam_database  
cd kofam_database  
wget -c ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz  
wget -c ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz  
gzip -d ko_list.gz  
tar xzf profiles.tar.gz; rm profiles.tar.gz  
cd profiles  
cp ../../Accessory_scripts/batch_hmmpress.pl ./  
perl batch_hmmpress.pl
cd ../
