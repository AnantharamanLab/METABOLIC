# Dockerfile for creating Metabolic in a container, part 2: perl and cpan packages

# Metabolic source: https://github.com/AnantharamanLab/METABOLIC
# Dockerfile github: https://github.com/tin6150/metabolic/blob/master/Dockerfile.perl


FROM tin6150/base4metabolic
MAINTAINER Tin (at) LBL.gov

ARG TZ="America/Los_Angeles"
ARG DEBIAN_FRONTEND=noninteractive

RUN touch    _TOP_DIR_OF_CONTAINER_  ;\
    echo "====================================== " | tee -a _TOP_DIR_OF_CONTAINER_  ;\
    echo "Begin Dockerfile.perl build process at " | tee -a _TOP_DIR_OF_CONTAINER_  ;\
    hostname | tee -a       _TOP_DIR_OF_CONTAINER_        ;\
    date     | tee -a       _TOP_DIR_OF_CONTAINER_        ;\
    apt-get -y --quiet update                             ;\
    apt-get -y --force-yes --quiet install perl-base perl-doc perl-doc-html libterm-shellui-perl cpanoutdated hmmer libbio-samtools-perl ;\
    apt-get -y --quiet install git                        ;\
    test -d /opt/gitrepo  || mkdir -p /opt/gitrepo        ;\
    cd      /opt/gitrepo  ;\
    test -d /opt/gitrepo/metabolic  || git clone https://github.com/tin6150/metabolic.git  ;\
    cd      /opt/gitrepo/metabolic &&  git pull && cd -   ;\
    apt-get -y --quiet clean                              ;\
    cd      /  ;\

    echo '==================================================================' ;\
    echo "installing perl/cpan packages"  | tee -a _TOP_DIR_OF_CONTAINER_     ;\  
    date | TZ=PST8PDT tee -a                       _TOP_DIR_OF_CONTAINER_     ;\  
    echo '==================================================================' ;\
    # TBD, future adopt tin6150/bioperl use of external shell scripts to do more complete install
    export PERL_MM_USE_DEFAULT=1                                              ;\
    # cpan -f is force, -i is install.  Bio::Perl is a beast and won't install :(
    PERL_MM_USE_DEFAULT=1 cpan -fi Data::Dumper                               ;\
    PERL_MM_USE_DEFAULT=1 perl -MCPAN -e cpann "install POSIX qw(strftime)"   ;\  
    PERL_MM_USE_DEFAULT=1 cpan -fi Excel::Writer::XLSX ;\    
    PERL_MM_USE_DEFAULT=1 cpan -fi Getopt::Long ;\    
    PERL_MM_USE_DEFAULT=1 cpan -fi Statistics::Descriptive ;\
    PERL_MM_USE_DEFAULT=1 perl -MCPAN -e cpann "Array::Split qw(split_by split_into)" ;\
    PERL_MM_USE_DEFAULT=1 cpan -fi Bio::SeqIO ;\
    # found these dependencies when tried to run METABOLIC-G.pl
    PERL_MM_USE_DEFAULT=1 cpan -fi Array/Split.pm ;\
    PERL_MM_USE_DEFAULT=1 cpan -fi Data/OptList.pm  ;\  
    PERL_MM_USE_DEFAULT=1 cpan -fi Parallel/ForkManager.pm  ;\  
    # 4.0 (5acb686) addition:
    PERL_MM_USE_DEFAULT=1 cpan -fi POSIX                               ;\

    echo '==================================================================' ;\
    echo '==================================================================' ;\
    echo "installing Bio::Perl packages"  | tee -a _TOP_DIR_OF_CONTAINER_     ;\  
    echo '==================================================================' ;\
    PERL_MM_USE_DEFAULT=1 cpan -fi Bio::Perl  ;\  
    echo $? > bioperl.exit.code ;\
    echo '==================================================================' ;\
    echo "done install Bio::Perl packages"  | tee -a _TOP_DIR_OF_CONTAINER_     ;\  
    echo '==================================================================' ;\
    echo '==================================================================' ;\
    PERL_MM_USE_DEFAULT=1 cpan -fi Bio::Tools::CodonTable ;\
    PERL_MM_USE_DEFAULT=1 cpan -fi Carp ;\
    PERL_MM_USE_DEFAULT=1 cpan -fi Text::Levenshtein::XS Text::Levenshtein::Damerau::XS Text::Levenshtein Text::Levenshtein::Damerau::PP ;\
    echo $? > cpan.exit.code ;\
    cpan -a > cpan.list.out ;\

    echo 'Ending.  Last line of RUN block in Dockerbuild without continuation :)'


RUN     cd / \
  && touch _TOP_DIR_OF_CONTAINER_  \
  && TZ=PST8PDT date  >> _TOP_DIR_OF_CONTAINER_  \
  && echo  "Dockerfile.perl 2021.0424"  >> _TOP_DIR_OF_CONTAINER_   

# ENV TZ could be changed/overwritten by container's /etc/csh.cshrc
ENV DOCKERFILE Dockerfile.perl # DOES overwrite parent def of ENV DOCKERFILE
ENV TEST_DOCKER_ENV_2   Can_use_ADD_to_make_ENV_avail_in_build_process
ENV TEST_DOCKER_ENV_REF https://vsupalov.com/docker-arg-env-variable-guide/#setting-env-values
ENV DOCKER_MANTABOLIC_PERL "CPAN packages, including bioperl"

ENTRYPOINT [ "/bin/bash" ]
