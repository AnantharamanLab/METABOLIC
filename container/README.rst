
Starting the Metabolic v4.0 container via Docker
================================================

Interactive run (note that content in container are ephemeral.  saving kinda work if you use ``docker commit ...``::

	cd METABOLIC; mkdir temp; tar xfz ~/Downloads/5_genomes_test.tgz
	docker run -it -e DISPLAY=$DISPLAY -v /tmp/.X11-unix:/tmp/.X11-unix -v "$PWD":/tmp/home --user=$(id -u):$(id -g)  tin6150/metabolic:4.0
	cd /opt/METABOLIC
	perl METABOLIC-G.pl -help
	perl /opt/METABOLIC/METABOLIC-G.pl -in-gn /tmp/home/5_genomes_test/Genome_files -o /tmp/home/metabolic_out 

Non interactive, scriptable run::

 
	docker pull tin6150/metabolic:4.0 
	docker run  -v "$(pwd)":/tmp/home --entrypoint "perl" tin6150/metabolic:4.0 /opt/METABOLIC/METABOLIC-G.pl -t 34 -in-gn /tmp/home/5_genomes_test/Genome_files -o /tmp/home/metabolic_out 
	# Output will be in ./metabolic_out



Building Docker container
=========================

Metabolic in docker container is build as 3 cascading parts:

1. Dockerfile.base: this first layer is a Debian Linux with R 3.6.2 (and CRAN libs) and number of .deb packages to satisfy dependencies (Hmmer, etc)

2. Dockerfile.perl: this add BioPerl and other CPAN libraries

3. Dockerfile.metabolic: this final layer add miniconda, sambamba, and METABOLIC software itself.
    - run_setup.sh has been run, software is installed under /opt/METABOLIC .

4. gtdbtk binary was installed, but not the database.  
   To use it, see "DB for gtdbtk" below.

5. /usr/bin/xpdf and /usr/bin/scp are available to view resulting pdf or copy data out of the container, if needed.


Build Commands
==============

::

		docker build -t tin6150/base4metabolic  -f Dockerfile.base       .  | tee Dockerfile.base.log 
		docker build -t tin6150/perl4metabolic  -f Dockerfile.perl       .  | tee Dockerfile.perl.log 
		docker build -t tin6150/metabolic       -f Dockerfile.metabolic  .  | tee Dockerfile.log 
		docker build -t tin6150/metabolic:4.0   -f Dockerfile.metabolic  .  | tee Dockerfile.log 
		docker login --username tin6150
		docker image push tin6150/metabolic:4.0  


Debug runs/tests
================

::

        docker run  -it -v $HOME:/home/tin tin6150/metabolic
        docker exec -it pluto_amp bash                 # additional terminal into existing running container

        # testing intermediary container use:
        docker run  -it -v $HOME:/home/tin tin6150/base4metabolic
        docker run  -it -v $HOME:/home/tin tin6150/perl4metabolic

        # checking PERL5LIB @INC
        env -i perl -V    # ignores the PERL5LIB env var
        env    perl -V
        # both should return the same output, 
        # but if root's env got inherited, clear it with something like export PERL5LIB=''


DB for gtdbtk 
=============

gtdbtk maybe optional.  when running it, may need a DB.  setup as:: 

	GTDBTK_DATA_PATH=/tmp/GTDBTK_DATA
	cd $GTDBTK_DATA_PATH
	wget https://data.ace.uq.edu.au/public/gtdb/data/releases/release89/89.0/gtdbtk_r89_data.tar.gz
	tar xzf gtdbtk_r89_data.tar.gz
	# See https://github.com/Ecogenomics/GTDBTk for links to newer db

	docker run  -v /tmp:/tmp --entrypoint "perl" tin6150/metabolic:4.0 /opt/METABOLIC/METABOLIC-G.pl -t 34 -in-gn /tmp//5_genomes_test/Genome_files -o /tmp/metabolic_out 




.. # vim: tabstop=4 noexpandtab paste background=dark
