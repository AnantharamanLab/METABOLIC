#!/usr/bin/perl

###########################

# METABOLIC-C.pl

# METABOLIC  =>  METabolic And BiogeOchemistry anaLyses In miCrobes

# This software gives a metabolic and biogeochemical function trait profile to given genome datasets 
# [either metagenome-assembled genomes (MAGs), single-cell amplified genomes (SAGs) or pure culture sequenced genomes]. 
# It also integrates the genome coverage to make element cycling pathways.
# METABOLIC-C.pl is specifically for users who have metagenomic reads and want to include them in the community analysis. 

# Written by Zhichao Zhou, zczhou2017@gmail.com 
# July, 2019

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

###########################

use 5.010;
use strict;
use warnings;

##modules
use Data::Dumper;
use POSIX qw(strftime);
use Getopt::Long;
use Statistics::Descriptive;
use Array::Split qw(split_by split_into);
use Parallel::ForkManager;
use File::Spec;
use File::Basename;

=head1 DESCRIPTION

        Takes a folder containing genome files to generate a profile of their metablic and biogeochemical functions
        
=head1 USAGE       

        perl METABOLIC-C.pl -t 40 -m-cutoff 0.75  -in Genome_proteins -kofam-db full -r omic_reads_parameters.txt -o METABOLIC_out
		(When you also want to calculate genome coverages, you would have to add genome files with the same basename and the extention as ".fasta" in this folder)
        
        perl METABOLIC-C.pl -t 40 -m-cutoff 0.75  -in-gn Genome_files -kofam-db full -r omic_reads_parameters.txt -o METABOLIC_out
		
		perl METABOLIC-G.pl -test true 
		use the 5 genomes to test the METABOLIC-C script
        
=head1 OPTIONS

	-t         or -cpu            [integer] The cpu numbers to run the hmmsearch (default: 20)
        -m-cutoff  or -module-cutoff  [float]   The cutoff value to assign the presence of a specific KEGG module (KEGG module step present numbers / KEGG module step total number) (default: 0.75) 
        -in                           [string]  The folder pf given genome faa files [should also give the genome fasta files and genone gene files if the (meta)genome/(meta)transciptome datasets are included]
        -in-gn                        [string]  The folder of given genome fasta files (Prodigal will be used to annotate your genomes)
        -kofam-db                     [string]  to use the "small" size or "full" size of KOfam database in METABOLIC (default: 'full')
	-p         or prodigal-method [string]  "meta" or "single" for prodigal to annotate the orf
        -r         or omic-reads      [string]  The file which indicates the address of omic reads
	    -rt        or reads-type      [string]  to use "metaG" or "metaT" to indicate whether you use the metagenomic reads or metatranscriptomic reads (default: 'metaG')
	-o         or output          [string]  The METABOLIC output folder (default: current address)
		
	
=head1 INSTRUCTIONS

	GitHub: https://github.com/AnantharamanLab/METABOLIC

=head1 OUTPUT

        Prodigal-annotated protein and gene files will be given in the input folder.
        The METABOLIC result table will be generated.
	Biogeochemical pathway diagrams will be generated.
		
=head1 COPYRIGHT

	Zhichao Zhou, zczhou2017@gmail.com
	Patricia Tran, ptran5@wisc.edu
	Karthik Anantharaman, karthik@bact.wisc.edu
	Anantharaman Microbiome Laboratory
	Department of Bacteriology, University of Wisconsin, Madison
=cut

#Intake the address of METABOLIC directory:
my $METABOLIC_dir = dirname(File::Spec->rel2abs(__FILE__));

#The options 
# option variables with default value
my $cpu_numbers = 20; #parallel running cpu numbers
my $module_cutoff = 0.75; # The cutoff value to assign the existence of a module
my $input_protein_folder;  #input microbial genome protein files
my $input_genome_folder; #input microbial genome fasta files
my $omic_reads_parameters; #The address of omic reads
my $prodigal_method = "meta"; #the prodigal method to annotate orfs
my $kofam_db_size = "full"; #the full kofam size
my $omic_reads_type = "metaG"; # metagenomic reads
my $output = `pwd`; # the output folder 
my $version="METABOLIC-C.pl v4.0";
my $test = "false";

GetOptions(
	'cpu|t=i' => \$cpu_numbers,
	'module-cutoff|m-cutoff=f' => \$module_cutoff,
	'in=s' => \$input_protein_folder,
	'in-gn=s' => \$input_genome_folder,
	'prodigal-method|p=s' => \$prodigal_method,
	'omic-reads|r=s' => \$omic_reads_parameters,
	'reads-type|rt=s' => \$omic_reads_type,
	'kofam-db=s' => \$kofam_db_size,
	'output|o=s' => \$output,
	'help|h' => sub{system('perldoc', $0); exit;},
	'v|version'=>sub{print $version."\n"; exit;},
	'test=s' => \$test
) or die("Getting options from the command line failed, please check your options");

##Pre-required files and documents
 #METABOLIC hmm database files
 my $METABOLIC_hmm_db_address = "$METABOLIC_dir/METABOLIC_hmm_db";
 
 #KofamKOALA hmm database files, July, 3 updated 
 #link: ftp://ftp.genome.jp/pub/db/kofam/
 my $kofam_db_address = "$METABOLIC_dir/kofam_database/profiles";
 my $kofam_db_KO_list = "$METABOLIC_dir/kofam_database/ko_list";
 
 #input hmm information table as a template
 my $hmm_table_temp = "$METABOLIC_dir/METABOLIC_template_and_database/hmm_table_template.txt";
 my $hmm_table_temp_2 = "$METABOLIC_dir/METABOLIC_template_and_database/hmm_table_template_2.txt"; 
 
 #the KEGG module information
 my $ko_module_table = "$METABOLIC_dir/METABOLIC_template_and_database/ko00002.keg";
 
 #the KEGG module step db 
 my $ko_module_step_db = "$METABOLIC_dir/METABOLIC_template_and_database/kegg_module_step_db.txt";
 
 #the pathway information to draw element cycling diagrams and metabolic handoff
 my $R_pathways = "$METABOLIC_dir/METABOLIC_template_and_database/R_pathways.txt";
 my $R_mh_01 = "$METABOLIC_dir/METABOLIC_template_and_database/Sequential_transformations_01.txt";
 my $R_mh_02 = "$METABOLIC_dir/METABOLIC_template_and_database/Sequential_transformations_02.txt";
 my $R_mh_tsv = "$METABOLIC_dir/METABOLIC_template_and_database/Sequential-transformations.tsv";
 my $R_order_of_input_01 = "$METABOLIC_dir/METABOLIC_template_and_database/order_of_input_01.txt";
 my $R_order_of_input_02 = "$METABOLIC_dir/METABOLIC_template_and_database/order_of_input_02.txt";
 my $CAZy_map_address = "$METABOLIC_dir/METABOLIC_template_and_database/CAZy_map.txt";
 
 #the MN-score reaction table template
 my $MN_score_reaction_table = "$METABOLIC_dir/METABOLIC_template_and_database/MN-score_reaction_table.txt";
 
 #the motif files to validate specific protein hits
 my $motif_file = "$METABOLIC_dir/METABOLIC_template_and_database/motif.txt";
 my $motif_pair_file = "$METABOLIC_dir/METABOLIC_template_and_database/motif.pair.txt";

# the test option:
if ($test eq "true"){
	$input_genome_folder = "$METABOLIC_dir/METABOLIC_test_files/Guaymas_Basin_genome_files";
	$output = "METABOLIC_out";
	$cpu_numbers = "5";
	$omic_reads_parameters = "$METABOLIC_dir/METABOLIC_test_files/Reads_address.txt";
}

##Main Body
#The present time

`mkdir $output`;

my $datestring = strftime "%Y-%m-%d %H:%M:%S", localtime; 
my $starttime = $datestring;

#Store the hmm table template
my %Hmm_table_temp = (); # line no. => each line 
my @Hmm_table_head = (); # the head of the hmm table template
my %METABOLIC_hmm2threshold = (); #hmm file id => threshold and score_type
open IN, "$hmm_table_temp";
while (<IN>){
	chomp;
	if (!/^#/){		
		my @tmp = split (/\t/);
		$Hmm_table_temp{$tmp[0]} = $_;
		if ($tmp[5] and $tmp[5] !~ /K\d\d\d\d\d/){
			$METABOLIC_hmm2threshold{$tmp[5]} = $tmp[10];
		}
	}else{
		my $line = $_; @Hmm_table_head = split (/\t/);		
	}
}
close IN;

#Store the hmm table template 2
my %Hmm_table_temp_2 = (); # line no. => each line; 
open IN, "$hmm_table_temp_2";
while (<IN>){
	chomp;
	if (!/^#/){	
		my @tmp = split (/\t/);
		$Hmm_table_temp_2{$tmp[0]}= $_;
	}
}
close IN;

#the hash of hmm file and corresponding threshold and score_type
my %Total_hmm2threshold = (%METABOLIC_hmm2threshold, _get_kofam_db_KO_threshold($kofam_db_KO_list,$kofam_db_address)); 

#run hmmsearch
#my $absolute_path = $ENV{'PWD'}; 
#my $folder_intermediate_files = "$absolute_path/intermediate_files";

`mkdir $output/intermediate_files`;

if ($input_genome_folder){
	open OUT, ">$output/tmp_run_annotate.sh";
	open OUT2, ">$output/tmp_run_annotate.sh.2";
	open IN, "ls $input_genome_folder/*.fasta |";
	$datestring = strftime "%Y-%m-%d %H:%M:%S", localtime; 
	print "\[$datestring\] The Prodigal annotation is running...\n";
	while (<IN>){
		chomp;
		my ($gn_id) = $_ =~ /^$input_genome_folder\/(.+?)\.fasta/;
		print OUT "prodigal -i $input_genome_folder/$gn_id.fasta -a $input_genome_folder/$gn_id.faa -o $input_genome_folder/$gn_id.gff -f gff -p $prodigal_method -q\n";
		print OUT2 "perl $METABOLIC_dir/Accessory_scripts/gff2fasta_mdf.pl -g $input_genome_folder/$gn_id.gff -f $input_genome_folder/$gn_id.fasta -o $input_genome_folder/$gn_id.gene\n";
	}
	close IN;
	close OUT;
	close OUT2;
	
	_run_parallel("$output/tmp_run_annotate.sh", $cpu_numbers); `rm $output/tmp_run_annotate.sh`;
	_run_parallel("$output/tmp_run_annotate.sh.2", $cpu_numbers); `rm $output/tmp_run_annotate.sh.2`;	

	$input_protein_folder = $input_genome_folder;
	$datestring = strftime "%Y-%m-%d %H:%M:%S", localtime; 
	print "\[$datestring\] The Prodigal annotation is finished\n";
} else {
	$input_genome_folder = $input_protein_folder;
}

my %Genome_id = (); # genome id => 1
my %Seqid2Genomeid = (); # seq id => genome id 
my %Total_faa_seq = (); #Store the total faa file into a hash{$line_no}
open IN,"ls $input_protein_folder/*.faa |";
while (<IN>){
	chomp;
	my $file = $_;
	#Store faa file into a hash
	%Total_faa_seq = (%Total_faa_seq, _get_faa_seq($file));
	
	my ($gn_id) = $file =~ /^$input_protein_folder\/(.+?)\.faa/; 
	$Genome_id{$gn_id} = 1; 
	open IN_, "$file";
	while (<IN_>){
		if (/>/){
			my ($seq) = $_ =~ /^>(.+?)\s/;
			$Seqid2Genomeid{$seq} = $gn_id;
		}
	}
	close IN_;
}

open OUT, ">$output/tmp_run_hmmsearch.sh";
`cat $input_protein_folder/*.faa > $input_protein_folder/faa.total; mv $input_protein_folder/faa.total $input_protein_folder/total.faa`;

`mkdir $output/intermediate_files/Hmmsearch_Outputs`;
open IN,"ls $input_protein_folder/total.faa |";
while (<IN>){
	chomp;
	my $file = $_;
	foreach my $hmm (sort keys %Total_hmm2threshold){
		my ($threshold,$score_type) = $Total_hmm2threshold{$hmm} =~ /^(.+?)\|(.+?)$/;
		#print "$hmm\t$threshold\t$score_type\n";
		if ($score_type eq "full"){
			if ($hmm !~ /K\d\d\d\d\d/){
				print OUT "hmmsearch -T $threshold --cpu 1 --tblout $output/intermediate_files/Hmmsearch_Outputs/$hmm.total.hmmsearch_result.txt $METABOLIC_hmm_db_address/$hmm $input_protein_folder/total.faa\n";
			}else{
				print OUT "hmmsearch -T $threshold --cpu 1 --tblout $output/intermediate_files/Hmmsearch_Outputs/$hmm.total.hmmsearch_result.txt $kofam_db_address/$hmm $input_protein_folder/total.faa\n";
			}
		}else{
			if ($hmm !~ /K\d\d\d\d\d/){
				print OUT "hmmsearch --domT $threshold --cpu 1 --tblout $output/intermediate_files/Hmmsearch_Outputs/$hmm.total.hmmsearch_result.txt $METABOLIC_hmm_db_address/$hmm $input_protein_folder/total.faa\n";
			}else{
				print OUT "hmmsearch --domT $threshold --cpu 1 --tblout $output/intermediate_files/Hmmsearch_Outputs/$hmm.total.hmmsearch_result.txt $kofam_db_address/$hmm $input_protein_folder/total.faa\n";
			}
		}
	}
}
close IN;
close OUT;

$datestring = strftime "%Y-%m-%d %H:%M:%S", localtime; 
print "\[$datestring\] The hmmsearch is running with $cpu_numbers cpu threads...\n";

#parallel run hmmsearch

_run_parallel("$output/tmp_run_hmmsearch.sh", $cpu_numbers); `rm $output/tmp_run_hmmsearch.sh`;
`rm $input_protein_folder/total.faa`;

$datestring = strftime "%Y-%m-%d %H:%M:%S", localtime; 
print "\[$datestring\] The hmmsearch is finished\n";

#Store motif validation files
my %Motif = _get_motif($motif_file); #protein id => motif sequences (DsrC => GPXKXXCXXXGXPXPXXCX)
my %Motif_pair = _get_motif_pair($motif_pair_file); # dsrC => tusE

#summarize hmmsearch result and print table
my %Hmmscan_result = (); # genome_name => hmm => numbers
my %Hmmscan_hits = (); # genome_name => hmm => hits
my %Hmm_id = (); # hmm => 1 
open IN, "find $output/intermediate_files/Hmmsearch_Outputs -type f -name '*.hmmsearch_result.txt' | ";
while (<IN>){
        chomp;
        my $file_name = $_;
        my ($hmm) = $file_name =~ /^$output\/intermediate_files\/Hmmsearch_Outputs\/(.+?\.hmm)\./; 
		$Hmm_id{$hmm} = 1;
		#my ($gn_id) = $file_name =~ /\.hmm\.(.+?)\.hmmsearch_result/;
		my $gn_id = "";
        open INN, "$file_name";
        while (<INN>){
                chomp;
                if (!/^#/){
                        my $line = $_; $line =~ s/\s+/\t/g;
						my @tmp = split (/\t/,$line); $gn_id = $Seqid2Genomeid{$tmp[0]};
						my ($threshold,$score_type) = $Total_hmm2threshold{$hmm} =~ /^(.+?)\|(.+?)$/; 
						if ($score_type eq "domain"){
							if ($tmp[8] >= $threshold){
								my ($hmm_basename) = $hmm =~ /^(.+?)\.hmm/; 
								if (exists $Motif{$hmm_basename}){
									my $seq; 
									my $motif = $Motif{$hmm_basename}; $motif =~ s/X/\[ARNDCQEGHILKMFPSTWYV\]/g; 									
									my %Seq_gn = _store_seq("$input_protein_folder/$gn_id.faa"); # get the genome sequences
									$seq = $Seq_gn{">$tmp[0]"};
									if ($seq =~ /$motif/){
										if (! exists $Hmmscan_hits{$gn_id}{$hmm}){
											$Hmmscan_hits{$gn_id}{$hmm} = $tmp[0];
										}else{
											$Hmmscan_hits{$gn_id}{$hmm} .= "\,".$tmp[0];
										}
										$Hmmscan_result{$gn_id}{$hmm}++;
									}
								}elsif(exists $Motif_pair{$hmm_basename}){
									my $motif_hmm = "$METABOLIC_hmm_db_address/$hmm_basename.check.hmm";
									my $motif_anti_hmm = "$METABOLIC_hmm_db_address/$Motif_pair{$hmm_basename}.check.hmm";
									_get_1_from_input_faa("$input_protein_folder/$gn_id.faa",">$tmp[0]","$output/tmp.$hmm_basename.check.faa");								
									`hmmsearch --cpu 1 --tblout $output/tmp.$hmm_basename.check.hmmsearch_result.txt $motif_hmm $output/tmp.$hmm_basename.check.faa`;
									`hmmsearch --cpu 1 --tblout $output/tmp.$Motif_pair{$hmm_basename}.check.hmmsearch_result.txt $motif_anti_hmm $output/tmp.$hmm_basename.check.faa`;
									my $motif_check_score = _get_check_score("$output/tmp.$hmm_basename.check.hmmsearch_result.txt"); 
									my $motif_anti_check_score = _get_check_score("$output/tmp.$Motif_pair{$hmm_basename}.check.hmmsearch_result.txt");
									if ($motif_check_score >= $motif_anti_check_score and $motif_check_score != 0){
										if (! exists $Hmmscan_hits{$gn_id}{$hmm}){
											$Hmmscan_hits{$gn_id}{$hmm} = $tmp[0];
										}else{
											$Hmmscan_hits{$gn_id}{$hmm} .= "\,".$tmp[0];
										}
										$Hmmscan_result{$gn_id}{$hmm}++;
									}
									`rm $output/tmp.$hmm_basename.check.faa $output/tmp.$hmm_basename.check.hmmsearch_result.txt $output/tmp.$Motif_pair{$hmm_basename}.check.hmmsearch_result.txt`;									
								}else{ # do not have motif check step
									if (! exists $Hmmscan_hits{$gn_id}{$hmm}){
										$Hmmscan_hits{$gn_id}{$hmm} = $tmp[0];
									}else{
										$Hmmscan_hits{$gn_id}{$hmm} .= "\,".$tmp[0];
									}
									$Hmmscan_result{$gn_id}{$hmm}++;
								}
							}
						}else{
							my ($hmm_basename) = $hmm =~ /^(.+?)\.hmm/; 
							if (exists $Motif{$hmm_basename}){
								my $seq; # the protein seq
								my $motif = $Motif{$hmm_basename};  $motif =~ s/X/\[ARNDCQEGHILKMFPSTWYV\]/g; 		
								my %Seq_gn = _store_seq("$input_protein_folder/$gn_id.faa"); # get the genome sequences
								$seq = $Seq_gn{">$tmp[0]"};
								if ($seq =~ /$motif/){
									if (! exists $Hmmscan_hits{$gn_id}{$hmm}){
										$Hmmscan_hits{$gn_id}{$hmm} = $tmp[0];
									}else{
										$Hmmscan_hits{$gn_id}{$hmm} .= "\,".$tmp[0];
									}
									$Hmmscan_result{$gn_id}{$hmm}++;
								}
							}elsif(exists $Motif_pair{$hmm_basename}){
								my $motif_hmm = "$METABOLIC_hmm_db_address/$hmm_basename.check.hmm";
								my $motif_anti_hmm = "$METABOLIC_hmm_db_address/$Motif_pair{$hmm_basename}.check.hmm";
								_get_1_from_input_faa("$input_protein_folder/$gn_id.faa",">$tmp[0]","$output/tmp.$hmm_basename.check.faa");
								`hmmsearch --cpu 1 --tblout $output/tmp.$hmm_basename.check.hmmsearch_result.txt $motif_hmm $output/tmp.$hmm_basename.check.faa`;
								`hmmsearch --cpu 1 --tblout $output/tmp.$Motif_pair{$hmm_basename}.check.hmmsearch_result.txt $motif_anti_hmm $output/tmp.$hmm_basename.check.faa`;
								my $motif_check_score = _get_check_score("$output/tmp.$hmm_basename.check.hmmsearch_result.txt"); 
								my $motif_anti_check_score = _get_check_score("$output/tmp.$Motif_pair{$hmm_basename}.check.hmmsearch_result.txt");
								if ($motif_check_score >= $motif_anti_check_score and $motif_check_score != 0){
									if (! exists $Hmmscan_hits{$gn_id}{$hmm}){
										$Hmmscan_hits{$gn_id}{$hmm} = $tmp[0];
									}else{
										$Hmmscan_hits{$gn_id}{$hmm} .= "\,".$tmp[0];
									}
									$Hmmscan_result{$gn_id}{$hmm}++;
								}
								`rm $output/tmp.$hmm_basename.check.faa $output/tmp.$hmm_basename.check.hmmsearch_result.txt $output/tmp.$Motif_pair{$hmm_basename}.check.hmmsearch_result.txt`;
							}else{
								if (! exists $Hmmscan_hits{$gn_id}{$hmm}){
									$Hmmscan_hits{$gn_id}{$hmm} = $tmp[0];
								}else{
									$Hmmscan_hits{$gn_id}{$hmm} .= "\,".$tmp[0];
								}
								$Hmmscan_result{$gn_id}{$hmm}++;
							}
						}						
                }
        }
        close INN;		
}
close IN;

#print out hmm result each tsv file
`mkdir $output/METABOLIC_result_each_spreadsheet`;

#print worksheet1
open OUT, ">$output/METABOLIC_result_each_spreadsheet/METABOLIC_result_worksheet1.tsv";
#print head
my @Hmm_table_head_worksheet1 = ();
for(my $i=0; $i<=9; $i++){
	push @Hmm_table_head_worksheet1, $Hmm_table_head[($i+1)];
}
foreach my $gn_id (sort keys %Genome_id){
	push @Hmm_table_head_worksheet1, "$gn_id Hmm presence";
	push @Hmm_table_head_worksheet1, "$gn_id Hit numbers";
	push @Hmm_table_head_worksheet1, "$gn_id Hits";
}
print OUT join("\t",@Hmm_table_head_worksheet1)."\n";

#print body
foreach my $line_no (sort keys %Hmm_table_temp){
	my $row = $line_no;
	my @Hmm_table_body_worksheet1 = ();
		
		my @tmp = split(/\t/,$Hmm_table_temp{$line_no});
		my $hmm = $tmp[5];
		for(my $i=0; $i<=9; $i++){
				push @Hmm_table_body_worksheet1, $tmp[($i+1)];				
		}
		
		foreach my $gn_id (sort keys %Genome_id){
			my $hmm_presence = "Absent";
			my $hit_num = 0;
			my @Hits = ();
			if ($hmm and $hmm !~ /\,\s/){
				if ($Hmmscan_result{$gn_id}{$hmm}){
					$hmm_presence = "Present";
					push @Hits, $Hmmscan_hits{$gn_id}{$hmm}; 
					$hit_num = $Hmmscan_result{$gn_id}{$hmm};
				}else{
					push @Hits, "None";
				}
			}elsif($hmm and $hmm =~ /\,\s/){
				my @tmp = split (/\,\s/,$hmm);
				my $sum = 0; 
				for(my $i=0; $i<=$#tmp; $i++){
					if ($Hmmscan_result{$gn_id}{$tmp[$i]}){
						$sum += $Hmmscan_result{$gn_id}{$tmp[$i]};
						push @Hits, $Hmmscan_hits{$gn_id}{$tmp[$i]};
					}else{
						push @Hits, "None";
					}
				}
				if ($sum){
					$hmm_presence = "Present";
				}
				$hit_num = $sum;
			}
			
			push @Hmm_table_body_worksheet1,$hmm_presence;
			push @Hmm_table_body_worksheet1,$hit_num;
			push @Hmm_table_body_worksheet1,join("\,",@Hits);
		}
		print OUT join("\t",@Hmm_table_body_worksheet1)."\n";

}
close OUT;

#print worksheet2
open OUT, ">$output/METABOLIC_result_each_spreadsheet/METABOLIC_result_worksheet2.tsv";
#print head
my @Hmm_table_head_worksheet2 = ();
for(my $i=0; $i<=3; $i++){
	push @Hmm_table_head_worksheet2, $Hmm_table_head[($i+1)];
}
foreach my $gn_id (sort keys %Genome_id){
	push @Hmm_table_head_worksheet2, "$gn_id Function presence";
}
print OUT join("\t",@Hmm_table_head_worksheet2)."\n";

#write the main body of hmm result to worksheet2
foreach my $line_no (sort keys %Hmm_table_temp_2){
	my @Hmm_table_body_worksheet2 = ();
	my @tmp_table_2 = split(/\t/,$Hmm_table_temp_2{$line_no});
	for(my $i=0; $i<=3; $i++){
		push @Hmm_table_body_worksheet2, $tmp_table_2[($i+2)];			
	}
	my $line_no_4_table_1 = $tmp_table_2[1];
	
	if ($line_no_4_table_1 !~ /\|\|/){
		my @tmp_table_1 = split(/\t/,$Hmm_table_temp{$line_no_4_table_1});
		my $hmm = $tmp_table_1[5];	
		foreach my $gn_id (sort keys %Genome_id){
			my $hmm_presence = "Absent";
			if ($hmm and $hmm !~ /\,\s/){
				if ($Hmmscan_result{$gn_id}{$hmm}){
					$hmm_presence = "Present";
				}
			}elsif($hmm and $hmm =~ /\,\s/){
				my @tmp = split (/\,\s/,$hmm);
				my $sum = 0; my @array_hit_num = ();
				for(my $i=0; $i<=$#tmp; $i++){
					if ($Hmmscan_result{$gn_id}{$tmp[$i]}){
						$sum += $Hmmscan_result{$gn_id}{$tmp[$i]};
					}
				}
				if ($sum){
					$hmm_presence = "Present";
				}
			}
			push @Hmm_table_body_worksheet2,$hmm_presence;	
		}		
	}else{
		my @array_line_no_4_table_1 = split (/\|\|/, $line_no_4_table_1);
		my @array_hmm = ();
		foreach my $line_no_4_table_1 (@array_line_no_4_table_1){
			my @tmp_table_1 = split(/\t/,$Hmm_table_temp{$line_no_4_table_1});
			my $hmm = $tmp_table_1[5];
			if ($hmm and $hmm !~ /\,\s/){
				push @array_hmm, $hmm;
			}elsif($hmm and $hmm =~ /\,\s/){
				my @tmp = split (/\,\s/,$hmm);
				foreach my $key (@tmp){
					push @array_hmm, $key;
				}
			}
		}
		foreach my $gn_id (sort keys %Genome_id){
			my $hmm_presence = "Absent";
			my $sum = 0;
			foreach my $hmm (@array_hmm){
				if ($Hmmscan_result{$gn_id}{$hmm}){
					$sum += $Hmmscan_result{$gn_id}{$hmm};
				}
			}
			if ($sum){
				$hmm_presence = "Present";
			}
			push @Hmm_table_body_worksheet2,$hmm_presence;
		}
	}
	print OUT join("\t",@Hmm_table_body_worksheet2)."\n";
}
close OUT;

#print out each hmm faa collection
$datestring = strftime "%Y-%m-%d %H:%M:%S", localtime; 
print "\[$datestring\] Generating each hmm faa collection...\n";

`mkdir $output/Each_HMM_Amino_Acid_Sequence`;

foreach my $hmm (sort keys %Hmm_id){
	my %Hmm_faa_seq = (); #Store the faa seqs in a hmm
	foreach my $gn_id (sort keys %Hmmscan_hits){
		if ($Hmmscan_hits{$gn_id}{$hmm}){			
			my @Hits = split (/\,/,$Hmmscan_hits{$gn_id}{$hmm});
			foreach my $hit (@Hits){
				my $seq_head = ">".$gn_id."~~".$hit;
				if (exists $Total_faa_seq{$seq_head}){
					$Hmm_faa_seq{$seq_head} = $Total_faa_seq{$seq_head}; #print "$Total_faa_seq{$seq_head}\n";
				}
			}

		}
	}
	if (%Hmm_faa_seq){
		open OUT, ">$output/Each_HMM_Amino_Acid_Sequence/$hmm.collection.faa";
		foreach my $key (sort keys %Hmm_faa_seq){
			print OUT "$key\n$Hmm_faa_seq{$key}\n";
		}
		close OUT;
	}
}

$datestring = strftime "%Y-%m-%d %H:%M:%S", localtime; 
print "\[$datestring\] Each hmm faa collection has been made\n";

#Do the KEGG module calculating
$datestring = strftime "%Y-%m-%d %H:%M:%S", localtime; 
print "\[$datestring\] The KEGG module result is calculating...\n";

#store the KEGG module table
my %Cat2module = ();  #module category => modules "\t"
my $head_cat2module = (); #Central carbohydrate metabolism
open IN, "$ko_module_table ";
while (<IN>){
	chomp;
	if (/^C    /){
		$head_cat2module = $_; $head_cat2module =~ s/^C    //g;
		if (!exists $Cat2module{$head_cat2module}){
			$Cat2module{$head_cat2module} = "";
		}
	}elsif(/^D      /){
		my ($tmp) = $_ =~ /^D      (M\d\d\d\d\d)/; 
		if (!exists $Cat2module{$head_cat2module}){
			$Cat2module{$head_cat2module} = $tmp;
		}else{
			$Cat2module{$head_cat2module} .= "\t".$tmp;
		}
	}
}
close IN;

#Store the KEGG module step database
my %KEGG_module = (); # M00804+01 => 0: K13327 "\t" 1: dTDP-D-forosamine biosynthesis 2: ko id
my %KEGG_module2step_number = (); #M00804 => 3
my %KEGG_module2name = (); #M00804 => dTDP-D-forosamine biosynthesis
open IN, "$ko_module_step_db";
while (<IN>){
	chomp;
	if (!/^name/){
		my @tmp = split (/\t/);
		my (@KOs) = $tmp[1] =~ /(K\d\d\d\d\d)/g;
		$KEGG_module{$tmp[2]}[0] = join("\t",@KOs);
		$KEGG_module{$tmp[2]}[1] = $tmp[0];
		$KEGG_module{$tmp[2]}[2] = $tmp[1];
		my ($module,$step_num) = $tmp[2] =~ /^(M.+?)\+(.+?)$/;
		$KEGG_module2step_number{$module} = $step_num;
		$KEGG_module2name{$module} = $tmp[0];
	}	
}
close IN;

#the hmm to ko id hash
my %Hmm2ko = _get_hmm_2_KO_hash(%Hmm_table_temp); # like: TIGR02694.hmm => K08355.hmm

#to see whether a module step exists for a given genome
my %Module_step_result = (); # M00804+01 => genome id => 1 / 0
foreach my $m_step (sort keys %KEGG_module){
	foreach my $gn_id (sort keys %Genome_id){
		$Module_step_result{$m_step}{$gn_id} = 0; 
		foreach my $hmm (sort keys %Hmm_id){
			my $hmm_new = ""; # transfer all the hmm id to ko id
			if (exists $Hmm2ko{$hmm}){
				$hmm_new = $Hmm2ko{$hmm}; 
			}elsif (!exists $Hmm2ko{$hmm} and $hmm =~ /^K\d\d\d\d\d/){
				$hmm_new = $hmm; 
			}
			my $hmm_new_wo_ext =  "";
			if ($hmm_new){
				($hmm_new_wo_ext) = $hmm_new =~ /^(.+?)\.hmm/;		
			}
			if ($Hmmscan_result{$gn_id}{$hmm} and $KEGG_module{$m_step}[0] =~ /$hmm_new_wo_ext/ and $hmm_new_wo_ext){
				$Module_step_result{$m_step}{$gn_id} = 1; 
			}
		}
	}
}

my %Module_result = (); # M00804 => genome name => Absent / Present
foreach my $module (sort keys  %KEGG_module2step_number){
	foreach my $gn_id (sort keys %Genome_id){
		my $present_no = 0; 		
		foreach my $module_step (sort keys %Module_step_result){
			if ($module_step =~ /$module/){
				$present_no += $Module_step_result{$module_step}{$gn_id};
			}
		}
		
		my $ratio = $present_no / $KEGG_module2step_number{$module};

		if ($ratio >= $module_cutoff){
			$Module_result{$module}{$gn_id} = "Present";
		}else{
			$Module_result{$module}{$gn_id} = "Absent";
		}
		
	}	
}

#print worksheet3
#write the head of hmm result to worksheet3
open OUT, ">$output/METABOLIC_result_each_spreadsheet/METABOLIC_result_worksheet3.tsv";
#print head
my @Worksheet3_head = ();
push @Worksheet3_head, "Module ID";
push @Worksheet3_head, "Module";
push @Worksheet3_head, "Module Category";

foreach my $gn_id (sort keys %Genome_id){
	push @Worksheet3_head, "$gn_id Module presence";
}
print OUT join("\t",@Worksheet3_head)."\n";

#print the worksheet3 result
foreach my $module (sort keys %Module_result){
	my @Worksheet3_body = ();
	push @Worksheet3_body, $module;
	push @Worksheet3_body, $KEGG_module2name{$module};
	my $cat_4_module = ""; # the category name for module
	foreach my $cat (sort keys %Cat2module){
		if ($Cat2module{$cat} =~ /$module/){
			$cat_4_module = $cat;
		}
	}
	push @Worksheet3_body, $cat_4_module;
	foreach my $gn_id (sort keys %Genome_id){
		push @Worksheet3_body, $Module_result{$module}{$gn_id};
	}
	print OUT join("\t",@Worksheet3_body)."\n";
}
close OUT;

#print worksheet4
#write the head of hmm result to worksheet4
open OUT, ">$output/METABOLIC_result_each_spreadsheet/METABOLIC_result_worksheet4.tsv";
#print head
my @Worksheet4_head = ();
push @Worksheet4_head, "Module step";
push @Worksheet4_head, "Module";
push @Worksheet4_head, "KO id";
push @Worksheet4_head, "Module Category";


foreach my $gn_id (sort keys %Genome_id){
	push @Worksheet4_head, "$gn_id Module step presence";
}
print OUT join("\t",@Worksheet4_head)."\n";

#print the worksheet4 result
foreach my $module_step (sort keys %Module_step_result){
	my @Worksheet4_body = ();
	push @Worksheet4_body, $module_step;
	my ($module) = $module_step =~ /^(M.+?)\+/;
	push @Worksheet4_body, $KEGG_module2name{$module};
	push @Worksheet4_body, $KEGG_module{$module_step}[2];
	my $cat_4_module = ""; # the category name for module
	foreach my $cat (sort keys %Cat2module){
		if ($Cat2module{$cat} =~ /$module/){
			$cat_4_module = $cat;
		}
	}
	push @Worksheet4_body, $cat_4_module; 
	foreach my $gn_id (sort keys %Genome_id){
		push @Worksheet4_body, $Module_result{$module}{$gn_id};
	}
	print OUT join("\t",@Worksheet4_body)."\n";
	
}
close OUT;

$datestring = strftime "%Y-%m-%d %H:%M:%S", localtime; 
print "\[$datestring\] The KEGG identifier \(KO id\) result is calculating...\n";

#print the KEGG KO hits
#my %Hmmscan_result = (); # genome_name => hmm => numbers
#my %Hmmscan_hits = (); # genome_name => hmm => hits
my %Hmmscan_result_for_KO = ();  # new gn id => KO => numbers
my %Hmmscan_hits_for_KO = (); # new gn id => KO => hits
my %New_hmmid = (); #KOs (without extension) => 1

foreach my $genome_name (sort keys %Hmmscan_result){
	my $gn = $genome_name; 
	foreach my $hmm (sort keys %Hmm_id){
		my $hmm_new = ""; # transfer all the hmm id to ko id
		if (exists $Hmm2ko{$hmm}){
			$hmm_new = $Hmm2ko{$hmm}; 
		}elsif (!exists $Hmm2ko{$hmm} and $hmm =~ /^K\d\d\d\d\d/){
			$hmm_new = $hmm; 
		}
		my $hmm_new_wo_ext =  "";
		if ($hmm_new){
			($hmm_new_wo_ext) = $hmm_new =~ /^(.+?)\.hmm/; $New_hmmid{$hmm_new_wo_ext} = 1;
			$Hmmscan_result_for_KO{$gn}{$hmm_new_wo_ext} = $Hmmscan_result{$genome_name}{$hmm};
			$Hmmscan_hits_for_KO{$gn}{$hmm_new_wo_ext} = $Hmmscan_hits{$genome_name}{$hmm};
		}
		
	}
}

`mkdir $output/KEGG_identifier_result`;
foreach my $gn_id (sort keys %Genome_id){
	open OUT1, ">$output/KEGG_identifier_result/$gn_id.result.txt";
	open OUT2, ">$output/KEGG_identifier_result/$gn_id.hits.txt";
	foreach my $hmmid (sort keys %New_hmmid){	
		my $new_gn_id = $gn_id;
		my $result = "";
		if ($Hmmscan_result_for_KO{$new_gn_id}{$hmmid}){
			$result = $Hmmscan_result_for_KO{$new_gn_id}{$hmmid}; 
		}
		print OUT1 "$hmmid\t$result\n";
		my $hits = "";
		if ($Hmmscan_hits_for_KO{$new_gn_id}{$hmmid}){
			$hits = $Hmmscan_hits_for_KO{$new_gn_id}{$hmmid}; 
		}
		print OUT2 "$hmmid\t$hits\n";
	}
	close OUT1;
	close OUT2;
}

$datestring = strftime "%Y-%m-%d %H:%M:%S", localtime; 
print "\[$datestring\] The KEGG identifier \(KO id\) seaching result is finished\n";


#run the dbCAN 
$datestring = strftime "%Y-%m-%d %H:%M:%S", localtime; 
print "\[$datestring\] Searching CAZymes by dbCAN2...\n";

`mkdir $output/intermediate_files/dbCAN2_Files`;
open OUT, ">$output/tmp_run_dbCAN2.sh";
open IN,"ls $input_protein_folder/*.faa |";
while (<IN>){
	chomp;
	my $file = $_;
	my ($gn_id) = $file =~ /^$input_protein_folder\/(.+?)\.faa/;
	print OUT "hmmscan --domtblout $output/intermediate_files/dbCAN2_Files/$gn_id.dbCAN2.out.dm --cpu 1 $METABOLIC_dir/dbCAN2/dbCAN-fam-HMMs.txt $file > $output/intermediate_files/dbCAN2_Files/$gn_id.dbCAN2.out;";
	print OUT "python $METABOLIC_dir/Accessory_scripts/hmmscan-parser-dbCANmeta.py $output/intermediate_files/dbCAN2_Files/$gn_id.dbCAN2.out.dm > $output/intermediate_files/dbCAN2_Files/$gn_id.dbCAN2.out.dm.ps\n";
}
close IN;
close OUT;

#parallel run dbCAN2
_run_parallel("$output/tmp_run_dbCAN2.sh", $cpu_numbers); `rm $output/tmp_run_dbCAN2.sh`;

my %dbCANout = (); #genome => hmmid => number
my %dbCANout2 = (); #genome => hmmid => hits
my %Hmm_dbCAN2_id = (); # hmm => 1 
open IN, "ls $output/intermediate_files/dbCAN2_Files/*.dbCAN2.out.dm.ps |";
while (<IN>)
{
	my $file = $_;
	my ($gn_id) = $file =~ /^$output\/intermediate_files\/dbCAN2_Files\/(.+?)\.dbCAN2\.out\.dm\.ps/;
    open INN, "$file";
	while (<INN>){
	   if (/^GH|^PL/){
               my @tmp = split(/\t/,$_);
               my ($hmmid) = $tmp[0] =~ /(\S+?)\.hmm/; 
			   my ($hmmid_p1,$hmmid_p2) = $hmmid =~ /^(\D+?)(\d+)/;
			   my $num=(sprintf "%03d", $hmmid_p2);
			   $hmmid = $hmmid_p1.$num;
			   $Hmm_dbCAN2_id{$hmmid} = 1; 
               my ($name) = $tmp[2];
               #my @tmp2 = split (/\|/, $name); $Genome_name{$tmp2[0]} = 1;
               $dbCANout{$gn_id}{$hmmid}++;
               if (!exists $dbCANout2{$gn_id}{$hmmid}){
                       $dbCANout2{$gn_id}{$hmmid} = $name;
               }else{
                       $dbCANout2{$gn_id}{$hmmid} .= "\;".$name;
               }
       }
	}
	close INN;
}
close IN;

$datestring = strftime "%Y-%m-%d %H:%M:%S", localtime; 
print "\[$datestring\] dbCAN2 searching is done\n";

#print worksheet5
#write the head of hmm result to worksheet5
open OUT, ">$output/METABOLIC_result_each_spreadsheet/METABOLIC_result_worksheet5.tsv";
#print head
my @Worksheet5_head = ();
push @Worksheet5_head, "CAZyme ID";

foreach my $gn_id (sort keys %Genome_id){
	push @Worksheet5_head, "$gn_id Hit numbers";
	push @Worksheet5_head, "$gn_id Hits";
}
print OUT join("\t",@Worksheet5_head)."\n";

#print the worksheet5 result
foreach my $hmmid (sort keys %Hmm_dbCAN2_id){
	my @Worksheet5_body = ();
	push @Worksheet5_body,$hmmid;	
	foreach my $gn_id (sort keys %Genome_id){
		my $hit_num = 0;
		if ($dbCANout{$gn_id}{$hmmid}){
			$hit_num = $dbCANout{$gn_id}{$hmmid};	
		}
		push @Worksheet5_body,$hit_num;	
		my $hits = "None";
		if ($dbCANout2{$gn_id}{$hmmid}){
			$hits = $dbCANout2{$gn_id}{$hmmid}; 
		}
		push @Worksheet5_body,$hits;	
	}
	print OUT join("\t",@Worksheet5_body)."\n";
}
close OUT;

#run the MEROPS
$datestring = strftime "%Y-%m-%d %H:%M:%S", localtime; 
print "\[$datestring\] Searching MEROPS peptidase...\n";

`mkdir $output/intermediate_files/MEROPS_Files`;
open OUT, ">$output/tmp_run_MEROPS.sh";
open IN,"ls $input_protein_folder/*.faa |";
while (<IN>){
	chomp;
	my $file = $_;
	my ($gn_id) = $file =~ /^$input_protein_folder\/(.+?)\.faa/;
	print OUT "diamond blastp -d $METABOLIC_dir/MEROPS/pepunit.db -q $file -o $output/intermediate_files/MEROPS_Files/$gn_id.MEROPSout.m8 -k 1 -e 1e-10 --query-cover 80 --id 50 --quiet -p 1 > /dev/null\n";
}
close IN;
close OUT;

#parallel run the MEROPS
_run_parallel("$output/tmp_run_MEROPS.sh", $cpu_numbers); `rm $output/tmp_run_MEROPS.sh`;

my %MEROPS_map; #MER id => all line
open IN, "$METABOLIC_dir/MEROPS/pepunit.lib";
while (<IN>){
	chomp;
	if (/>/){
		$_ =~ tr/\015//d;
		my ($mer_id) = $_ =~ /^>(.+?)\s/;
		$MEROPS_map{$mer_id} = $_; 
	}
}
close IN;

my %MEROPSout = (); #genome => hmmid => number
my %MEROPSout2 = (); #genome => hmmid => hits
my %MEROPSid = (); # merops_id => 1 
open IN, "ls $output/intermediate_files/MEROPS_Files/*.MEROPSout.m8 |";
while (<IN>)
{
	my $file = $_;
	my ($gn_id) = $file =~ /^$output\/intermediate_files\/MEROPS_Files\/(.+?)\.MEROPSout\.m8/;
    open INN, "$file";
	while (<INN>){
	   my @tmp = split(/\t/,$_); 
	   my ($meropsid) = $MEROPS_map{$tmp[1]} =~ /\#(.+?)\#/; $MEROPSid{$meropsid} = 1;
	   my ($name) = $tmp[0];
	   $MEROPSout{$gn_id}{$meropsid}++;
	   if (!exists $MEROPSout2{$gn_id}{$meropsid}){
                $MEROPSout2{$gn_id}{$meropsid} = $name;
        }else{
                $MEROPSout2{$gn_id}{$meropsid} .= "\;".$name;
        }
	}
	close INN;
}
close IN;

$datestring = strftime "%Y-%m-%d %H:%M:%S", localtime; 
print "\[$datestring\] MEROPS peptidase searching is done\n";

#print worksheet6
#write the head of hmm result to worksheet6
open OUT, ">$output/METABOLIC_result_each_spreadsheet/METABOLIC_result_worksheet6.tsv";
#print head
my @Worksheet6_head = ();
push @Worksheet6_head, "MEROPS peptidase ID";

foreach my $gn_id (sort keys %Genome_id){
	push @Worksheet6_head, "$gn_id Hit numbers";
	push @Worksheet6_head, "$gn_id Hits";
}
print OUT join("\t",@Worksheet6_head)."\n";

#print the worksheet6 result
foreach my $hmmid (sort keys %MEROPSid){
	my @Worksheet6_body = ();
	push @Worksheet6_body, $hmmid;
	foreach my $gn_id (sort keys %Genome_id){
		my $hit_num = 0;
		if ($MEROPSout{$gn_id}{$hmmid}){
			$hit_num = $MEROPSout{$gn_id}{$hmmid};	
		}
		push @Worksheet6_body, $hit_num;	
		my $hits = "None";
		if ($MEROPSout2{$gn_id}{$hmmid}){
			$hits = $MEROPSout2{$gn_id}{$hmmid}; 
		}
		push @Worksheet6_body, $hits;	
	}
	print OUT join("\t",@Worksheet6_body)."\n";
}
close OUT;

`Rscript $METABOLIC_dir/create_excel_spreadsheet.R $output/METABOLIC_result_each_spreadsheet/ > /dev/null`;
`mv METABOLIC_result.xlsx $output/METABOLIC_result.xlsx`;

$datestring = strftime "%Y-%m-%d %H:%M:%S", localtime; 
print "\[$datestring\] METABOLIC table has been generated\n";

#Draw element cycling diagrams
$datestring = strftime "%Y-%m-%d %H:%M:%S", localtime; 
print "\[$datestring\] Drawing element cycling diagrams...\n";

#Store R pathways
my %R_pathways = (); #step => hmms
my %R_hmm_ids = ();
open IN, "$R_pathways";
while (<IN>){
	chomp;
	my @tmp = split (/\t/);
	$R_pathways{$tmp[0]} = $tmp[1];
	if ($tmp[1] !~ /\;/){
		my @tmp2 = split (/\,/,$tmp[1]);
		foreach my $key (@tmp2){
			$R_hmm_ids{$key} = 1;
		}
	}elsif ($tmp[1] =~ /\;/){
		my @tmp2 = split (/\;/,$tmp[1]);
		foreach my $key (@tmp2){
			my @tmp3 = split (/\,/,$key);
			foreach my $key2 (@tmp3){
				if ($key2 !~ /NO/){
					$R_hmm_ids{$key2} = 1;
				}
			}
		}
	}
}
close IN;

`mkdir $output/METABOLIC_Figures_Input`;
`mkdir $output/METABOLIC_Figures_Input/Nutrient_Cycling_Diagram_Input`;
`mkdir $output/METABOLIC_Figures`;

#Get each R pathway input files
my %Total_R_input = (); #pathway => gn => 1 or 0
foreach my $gn (sort keys %Hmmscan_result){
	my %R_input = (); #for each input file
	foreach my $key (sort keys %R_pathways){
		$R_input{$key} = 0; $Total_R_input{$key}{$gn} = 0;
		my $hmms = $R_pathways{$key};
		if ($hmms !~ /\;/){
			foreach my $hmm_id (sort keys %R_hmm_ids){
				if ($hmms =~ /$hmm_id/ and $Hmmscan_result{$gn}{$hmm_id}){
					$R_input{$key} = 1; $Total_R_input{$key}{$gn} = 1;
				}
			}
		}elsif ($hmms =~ /\;/){
			my ($hmms_1,$hmms_2) = $hmms =~ /^(.+?)\;(.+?)$/;
			if ($hmms_2 !~ /NO/){
				my $logic1 = 0; my $logic2 = 0;
				foreach my $hmm_id (sort keys %R_hmm_ids){
					if ($hmms_1 =~ /$hmm_id/ and $Hmmscan_result{$gn}{$hmm_id}){
						$logic1 = 1;	
					}
					if ($hmms_2 =~ /$hmm_id/ and $Hmmscan_result{$gn}{$hmm_id}){
						$logic2 = 1;	
					}
				}
				if ($logic1 and $logic2){
					$R_input{$key} = 1; $Total_R_input{$key}{$gn} = 1;
				}
			}elsif ($hmms_2 =~ /NO/){
				my $logic1 = 0; my $logic2 = 0;
				foreach my $hmm_id (sort keys %R_hmm_ids){
					if ($hmms_1 =~ /$hmm_id/ and $Hmmscan_result{$gn}{$hmm_id}){
						$logic1 = 1;	
					}			
					if ($hmms_2 !~ /$hmm_id/ and $Hmmscan_result{$gn}{$hmm_id}){
						$logic2 = 1;	
					}		
				}
				if ($logic1 and $logic2){
					$R_input{$key} = 1; $Total_R_input{$key}{$gn} = 1;
				}
			}
		}
	}
	
	open OUT, ">$output/METABOLIC_Figures_Input/Nutrient_Cycling_Diagram_Input/$gn.R_input.txt";
	foreach my $key (sort keys %R_input){
		print OUT "$key\t$R_input{$key}\n";
	}	
	close OUT;
}

my %Genome_cov_constant = ();
#The genome coverage: genome id => coverage value
if ($omic_reads_parameters){
	my %Genome_cov = _get_Genome_coverge($omic_reads_parameters,$input_genome_folder);
	%Genome_cov_constant = %Genome_cov;
	
	my %Total_R_input_2 = (); #pathway => genome numbers \t genome coverage percentage
	foreach my $pth (sort keys %Total_R_input){
		my $gn_no = 0; my $gn_cov_percentage = 0;
		foreach my $gn (sort keys %Hmmscan_result){
			if ($Total_R_input{$pth}{$gn}){
			$gn_no += $Total_R_input{$pth}{$gn};
			if ($Genome_cov{$gn}){
				$gn_cov_percentage += $Genome_cov{$gn};
				}
			}
		}
		$Total_R_input_2{$pth} = "$gn_no\t$gn_cov_percentage";
	}
	
	open OUT, ">$output/METABOLIC_Figures_Input/Nutrient_Cycling_Diagram_Input/Total.R_input.txt"; 
	foreach my $key (sort keys %Total_R_input_2){
		print OUT "$key\t$Total_R_input_2{$key}\n";
	}
	close OUT;
}

`Rscript $METABOLIC_dir/draw_biogeochemical_cycles.R $output/METABOLIC_Figures_Input/Nutrient_Cycling_Diagram_Input $output/Output TRUE > /dev/null`;
`mv $output/Output/draw_biogeochem_cycles $output/METABOLIC_Figures/Nutrient_Cycling_Diagrams; rm -r $output/Output`;

$datestring = strftime "%Y-%m-%d %H:%M:%S", localtime; 
print "\[$datestring\] Drawing element cycling diagrams finished\n";

my $endtime = "";
if ($omic_reads_parameters){

#Draw metabolic handoff diagrams
$datestring = strftime "%Y-%m-%d %H:%M:%S", localtime; 
print "\[$datestring\] Drawing metabolic handoff diagrams...\n";

#Store metabolic handoff steps 1
my %R_mh_01 = (); #step => KOs
my %R_mh_01_hmm_ids = ();
my %Letter2reaction = (); # D -> S2-=>S0
open IN, "$R_mh_01";
while (<IN>){
	chomp;
	if (/\:/ or /\+/){
		my @tmp = split (/\t/);
		
		if ($tmp[0] =~ /\:/){
			my ($step,$reaction) = $tmp[0] =~ /^(\D+?)\:(.+?)$/;
			$R_mh_01{$step} = $tmp[2]; $Letter2reaction{$step} = $reaction; 
			my @tmp2 = split (/\,/,$tmp[2]);
			foreach my $key (@tmp2){
				$R_mh_01_hmm_ids{$key} = 1;
			}						
		}else{
			$R_mh_01{$tmp[0]} = $tmp[2];
		}
	}
}
close IN;

my %R_mh_01_summary = (); #step => gn => 1 or 0
foreach my $gn (sort keys %Hmmscan_result){
	foreach my $key (sort keys %R_mh_01){
		$R_mh_01_summary{$key}{$gn} = 0;
		if ($key !~ /\+/){
			my $hmms = $R_mh_01{$key};
			foreach my $hmm_id (sort keys %R_mh_01_hmm_ids){
				if ($hmms =~ /$hmm_id/ and $Hmmscan_result{$gn}{$hmm_id}){
					$R_mh_01_summary{$key}{$gn} = 1;
				}
			}
		}else{
			my $hmms_total = $R_mh_01{$key};
			my @tmp_hmm = split (/\;/,$hmms_total);
			my $count = 0; my $tmp_hmm_scalar = scalar @tmp_hmm;
			foreach my $hmms (@tmp_hmm){
				my $logic = 0;
				foreach my $hmm_id (sort keys %R_mh_01_hmm_ids){
					if ($hmms =~ /$hmm_id/ and $Hmmscan_result{$gn}{$hmm_id}){
						$logic = 1 ;
					}
				}
				if ($logic){
					$count++;
				}
			}
			if ($tmp_hmm_scalar == $count){
				$R_mh_01_summary{$key}{$gn} = 1;
			}
		}
	}
}

#The genome coverage: genome id => coverage value
if ($omic_reads_parameters){
	my %Genome_cov = %Genome_cov_constant;
	
	my %Total_R_hm_input_1 = (); #step => genome numbers \t genome coverage percentage
	foreach my $step (sort keys %R_mh_01){
		my $gn_no = 0; my $gn_cov_percentage = 0;
		foreach my $gn (sort keys %Hmmscan_result){
			if ($R_mh_01_summary{$step}{$gn}){
				$gn_no += $R_mh_01_summary{$step}{$gn};
			if ($Genome_cov{$gn}){
				$gn_cov_percentage += $Genome_cov{$gn};
				}
			}
		}
		$Total_R_hm_input_1{$step} = "$gn_no\t$gn_cov_percentage";
	}

	open OUT, ">$output/METABOLIC_Figures_Input/Sequential_Transformation_input_1.txt"; 
	foreach my $key (sort keys %Total_R_hm_input_1){
		print OUT "$key\t$Total_R_hm_input_1{$key}\n";
	}
	close OUT;
}

#Store the CAZy map
my %CAZy_map = (); #GH => enzymes
open IN, "$CAZy_map_address";
while (<IN>){
	chomp;
	if (!/^Family/){
		my @tmp = split (/\t/);
		$CAZy_map{$tmp[0]} = $tmp[1];
	}
}
close IN;

#Store metabolic handoff steps 2
my %R_mh_02 = (); #step => enzymes
my %R_mh_02_enzyme_ids = ();
#my %Letter2reaction = (); # D -> S2-=>S0
open IN, "$R_mh_02";
while (<IN>){
	chomp;
	if (/\:/ or /\+/){
		my @tmp = split (/\t/);
		
		if ($tmp[0] =~ /\:/){
			my ($step,$reaction) = $tmp[0] =~ /^(\D+?)\:(.+?)$/;
			$R_mh_02{$step} = $tmp[2]; $Letter2reaction{$step} = $reaction;
			my @tmp2 = split (/\;\s/,$tmp[2]);
			foreach my $key (@tmp2){
				$R_mh_02_enzyme_ids{$key} = 1;
			}						
		}else{
			$R_mh_02{$tmp[0]} = $tmp[2];
		}	
	}
}
close IN;

my %R_mh_02_summary = (); #step => gn => 1 or 0
foreach my $step (sort keys %R_mh_02){
	foreach my $gn (sort keys %dbCANout){
		$R_mh_02_summary{$step}{$gn} = 0;
		if ($step !~ /\+/){
			my $enzymes = $R_mh_02{$step};
			my @tmp_enzymes = split (/\;/, $enzymes);
			foreach my $enzyme_id (@tmp_enzymes){
				foreach my $GH_id (sort keys %CAZy_map){
					my @tmp_enzymes_2 = split (/\;/,$CAZy_map{$GH_id});
					foreach my $enzyme_id_2 (@tmp_enzymes_2){
						if ($enzyme_id_2 eq $enzyme_id and $dbCANout{$gn}{$GH_id}){
							$R_mh_02_summary{$step}{$gn} = 1;
						}
					}
				}
			}
		}else{
			my $hmms_total = $R_mh_02{$step};
			my @tmp_hmm = split (/\|/,$hmms_total);
			my $count = 0; my $tmp_hmm_scalar = scalar @tmp_hmm;			
			foreach my $hmms (@tmp_hmm){
						my $logic = 0;
						my @tmp_enzymes = split (/\;/, $hmms);
						foreach my $enzyme_id (@tmp_enzymes){
							foreach my $GH_id (sort keys %CAZy_map){
								my @tmp_enzymes_2 = split (/\;/,$CAZy_map{$GH_id});
								foreach my $enzyme_id_2 (@tmp_enzymes_2){
									if ($enzyme_id_2 eq $enzyme_id and $dbCANout{$gn}{$GH_id}){
										$logic = 1;
									}
								}				
							}
						}
						if ($logic){
							$count++;
						}
				
			}
			if ($tmp_hmm_scalar == $count){
				$R_mh_02_summary{$step}{$gn} = 1;
			}
		}
	}
}

#The genome coverage: genome id => coverage value
if ($omic_reads_parameters){
	my %Genome_cov = %Genome_cov_constant;
	
	my %Total_R_hm_input_2 = (); #step => genome numbers \t genome coverage percentage
	foreach my $step (sort keys %R_mh_02){
		my $gn_no = 0; my $gn_cov_percentage = 0;
		foreach my $gn (sort keys %Hmmscan_result){
			if ($R_mh_02_summary{$step}{$gn}){
				$gn_no += $R_mh_02_summary{$step}{$gn};
			if ($Genome_cov{$gn}){
				$gn_cov_percentage += $Genome_cov{$gn};
				}
			}
		}
		$Total_R_hm_input_2{$step} = "$gn_no\t$gn_cov_percentage";
	}

	open OUT, ">$output/METABOLIC_Figures_Input/Sequential_Transformation_input_2.txt"; 
	foreach my $key (sort keys %Total_R_hm_input_2){
		print OUT "$key\t$Total_R_hm_input_2{$key}\n";
	}
	close OUT;
}

`mkdir $output/newdir`;
`Rscript $METABOLIC_dir/draw_sequential_reaction.R $output/METABOLIC_Figures_Input/Sequential_Transformation_input_1.txt $output/METABOLIC_Figures_Input/Sequential_Transformation_input_2.txt $R_mh_tsv $R_order_of_input_01 $R_order_of_input_02 $output/newdir > /dev/null 2>/dev/null`;
`mv $output/newdir/Bar_plot/bar_plot_input_1.pdf $output/METABOLIC_Figures/Sequential_transformation_01.pdf`;
`mv $output/newdir/Bar_plot/bar_plot_input_2.pdf $output/METABOLIC_Figures/Sequential_transformation_02.pdf`;
`rm -r $output/newdir`;

$datestring = strftime "%Y-%m-%d %H:%M:%S", localtime; 
print "\[$datestring\] Drawing metabolic handoff diagrams finished\n";

$datestring = strftime "%Y-%m-%d %H:%M:%S", localtime; 
print "\[$datestring\] Drawing energy flow chart...\n";

#store the bin category
system ("gtdbtk classify_wf --cpus $cpu_numbers -x fasta --genome_dir $input_genome_folder --out_dir $output/intermediate_files/gtdbtk_Genome_files > /dev/null");

my %Bin2Cat = (); # bin => category, for instance, Acidimicrobiia_bacterium_UWMA-0264 => Actinobacteriota
my %Cat = ();
if (-e "$output/intermediate_files/gtdbtk_Genome_files/gtdbtk.bac120.summary.tsv"){
	open IN, "$output/intermediate_files/gtdbtk_Genome_files/gtdbtk.bac120.summary.tsv";
	while (<IN>){
		chomp;
		if (!/^user_genome/){
			my @tmp = split (/\t/);
			if ($tmp[1] =~ /p\_\_Proteobacteria/){
				my ($cat) = $tmp[1] =~ /\;c\_\_(.+?)\;/;
				$Bin2Cat{$tmp[0]} = $cat;
				$Cat{$cat} = 1;
			}else{
				my ($cat) = $tmp[1] =~ /\;p\_\_(.+?)\;/;
				$Bin2Cat{$tmp[0]} = $cat;
				$Cat{$cat} = 1;
			}
		}
	}
	close IN;
}

if (-e "$output/intermediate_files/gtdbtk_Genome_files/gtdbtk.ar122.summary.tsv"){
	open IN, "$output/intermediate_files/gtdbtk_Genome_files/gtdbtk.ar122.summary.tsv";
	while (<IN>){
		chomp;
		if (!/^user_genome/){
			my @tmp = split (/\t/);
			my ($cat) = $tmp[1] =~ /\;p\_\_(.+?)\;/;
			$Bin2Cat{$tmp[0]} = $cat;
			$Cat{$cat} = 1;
		}
	}
	close IN;
}

my %Hash_gn_n_pth = (); 
my %Total_R_community_coverage = (); # genome\tpathway => category \t pathway \t genome coverage percentage
if ($omic_reads_parameters){
	my %Genome_cov = %Genome_cov_constant;
	#%Total_R_input pathway => gn => 1 or 0
	foreach my $pth (sort keys %Total_R_input){
		my $gn_cov_percentage = 0;
		foreach my $gn (sort keys %Hmmscan_result){
			if ($Genome_cov{$gn} and $Total_R_input{$pth}{$gn}){
				$gn_cov_percentage = $Genome_cov{$gn};
				my $cat = $Bin2Cat{$gn};
				my $gn_n_pth = "$gn\t$pth"; $Hash_gn_n_pth{$gn_n_pth} = 1;
				$Total_R_community_coverage{$gn_n_pth} = "$cat\t$pth\t$gn_cov_percentage";
			}
		}
	}		
}

my %Total_R_community_coverage2 = (); #$genome\tpath pair => cat \t  coverage percentage average
foreach my $gn (sort keys %Hmmscan_result){
	my %Path = (); # path => 1
	foreach my $gn_n_pth (sort keys %Total_R_community_coverage){
		if ($gn_n_pth =~ /$gn\t/){
			my @tmp = split (/\t/,$gn_n_pth);
			$Path{$tmp[1]} = 1;
		}
	}
	my @Path_keys = sort keys %Path;
	for(my $i=0; $i<=$#Path_keys; $i++){
		for(my $j = $i+1; $j<=$#Path_keys; $j++){
			my $pair = "$Path_keys[$i]\t$Path_keys[$j]";
			my $coverage = 0;
			my @tmp1 = split (/\t/, $Total_R_community_coverage{"$gn\t$Path_keys[$i]"});
			my @tmp2 = split (/\t/, $Total_R_community_coverage{"$gn\t$Path_keys[$j]"});
			$coverage = ($tmp1[2] + $tmp2[2]) / 2;
			$Total_R_community_coverage2{"$gn\t$pair"} = $Bin2Cat{$gn}."\t".$coverage;
		}
	}
}

open OUT, ">$output/METABOLIC_Figures_Input/Metabolic_energy_flow_input.txt";
foreach my $gn_n_pth (sort keys %Hash_gn_n_pth){
	print OUT "$Total_R_community_coverage{$gn_n_pth}\n";
}
close OUT;

open OUT, ">$output/METABOLIC_Figures_Input/Metabolic_network_input.txt";
print OUT "#Genome\tStep1\tStep2\tTaxonomic Group\tCoverage value\(average\)\n";
foreach my $gn_n_pair (sort keys %Total_R_community_coverage2){
	print OUT "$gn_n_pair\t$Total_R_community_coverage2{$gn_n_pair}\n";
}
close OUT;

`Rscript $METABOLIC_dir/draw_metabolic_energy_flow.R $output/METABOLIC_Figures_Input/Metabolic_energy_flow_input.txt $output/Output_energy_flow > /dev/null`;
`mv $output/Output_energy_flow/Energy_plot/network.plot.pdf   $output/METABOLIC_Figures/Metabolic_energy_flow.pdf; rm -r $output/Output_energy_flow`;

`Rscript $METABOLIC_dir/draw_metabolic_network.R $output/METABOLIC_Figures_Input/Metabolic_network_input.txt $output/OutputFolder_Energy > /dev/null`;
`mv $output/OutputFolder_Energy/network_plot $output/METABOLIC_Figures/Metabolic_network_figures; rm -r $output/OutputFolder_Energy`;

$datestring = strftime "%Y-%m-%d %H:%M:%S", localtime; 
print "\[$datestring\] Drawing energy flow chart finished\n";


$datestring = strftime "%Y-%m-%d %H:%M:%S", localtime; 
print "\[$datestring\] Calculating MN-score ...\n";

###To calculate MN-score
#Store MN_score_reaction_table
my %MN_functions = (); #func => hmms
my %MN_function_hmm_ids = ();
open IN, "$MN_score_reaction_table"; # read the MN_score_reaction table
while (<IN>){
	chomp;
	if (!/^#/){
		my @tmp = split (/\t/);
		$MN_functions{$tmp[0]} = $tmp[2];
		if ($tmp[2] !~ /\;/){
			my @tmp2 = split (/\,/,$tmp[2]);
			foreach my $key (@tmp2){
				$MN_function_hmm_ids{$key} = 1;
			}
		}elsif ($tmp[2] =~ /\;/){
			my @tmp2 = split (/\;/,$tmp[2]);
			foreach my $key (@tmp2){
				my @tmp3 = split (/\,/,$key);
				foreach my $key2 (@tmp3){
					if ($key2 !~ /NO/){
						$MN_function_hmm_ids{$key2} = 1;
					}
				}
			}
		}
	}
}
close IN;

#Get MN_score_reaction_table result
my %MN_score_hash = (); #pathway => gn => 1 or 0
foreach my $gn (sort keys %Hmmscan_result){
	foreach my $key (sort keys %MN_functions){
		$MN_score_hash{$key}{$gn} = 0;
		my $hmms = $MN_functions{$key};
		if ($hmms !~ /\;/){
			foreach my $hmm_id (sort keys %MN_function_hmm_ids){
				if ($hmms =~ /$hmm_id/ and $Hmmscan_result{$gn}{$hmm_id}){
					$MN_score_hash{$key}{$gn} = 1; 
				}
			}
		}elsif ($hmms =~ /\;/){
			my ($hmms_1,$hmms_2) = $hmms =~ /^(.+?)\;(.+?)$/;
			if ($hmms_2 !~ /NO/){
				my $logic1 = 0; my $logic2 = 0;  #print "$gn\n";
				foreach my $hmm_id (sort keys %MN_function_hmm_ids){
					if ($hmms_1 =~ /$hmm_id/ and $Hmmscan_result{$gn}{$hmm_id}){
						$logic1 = 1;	
					}
					if ($hmms_2 =~ /$hmm_id/ and $Hmmscan_result{$gn}{$hmm_id}){
						$logic2 = 1;	
					}
				}
				if ($logic1 and $logic2){
					$MN_score_hash{$key}{$gn} = 1;
				}
			}elsif ($hmms_2 =~ /NO/){
				my $logic1 = 0; my $logic2 = 0;
				foreach my $hmm_id (sort keys %MN_function_hmm_ids){
					if ($hmms_1 =~ /$hmm_id/ and $Hmmscan_result{$gn}{$hmm_id}){
						$logic1 = 1;	
					}			
					if ($hmms_2 !~ /$hmm_id/ and $Hmmscan_result{$gn}{$hmm_id}){
						$logic2 = 1;	
					}		
				}
				if ($logic1 and $logic2){
					$MN_score_hash{$key}{$gn} = 1;
				}
			}
		}
	}
}

my %MN_score_community_coverage = (); # genome\tpathway => category \t pathway \t genome coverage percentage
if ($omic_reads_parameters){
	my %Genome_cov = %Genome_cov_constant;
	foreach my $pth (sort keys %MN_score_hash){
		my $gn_cov_percentage = 0; 
		foreach my $gn (sort keys %Hmmscan_result){
			if ($Genome_cov{$gn} and $MN_score_hash{$pth}{$gn}){
				$gn_cov_percentage = $Genome_cov{$gn}; 
				my $cat = $Bin2Cat{$gn};
				my $gn_n_pth = "$gn\t$pth"; 
				$MN_score_community_coverage{$gn_n_pth} = "$cat\t$pth\t$gn_cov_percentage";
				
			}
		}
	}		
}

my %MN_score_community_coverage2 = (); #$genome\tpath pair => cat \t  coverage percentage average
foreach my $gn (sort keys %Hmmscan_result){
	my %Path = (); # path => 1
	foreach my $gn_n_pth (sort keys %MN_score_community_coverage){
		if ($gn_n_pth =~ /$gn\t/){
			my @tmp = split (/\t/,$gn_n_pth);
			$Path{$tmp[1]} = 1;
		}
	}
	my @Path_keys = sort keys %Path;
	for(my $i=0; $i<=$#Path_keys; $i++){
		for(my $j = $i+1; $j<=$#Path_keys; $j++){
			my $pair = "$Path_keys[$i]\t$Path_keys[$j]";
			my $coverage = 0;
			my @tmp1 = split (/\t/, $MN_score_community_coverage{"$gn\t$Path_keys[$i]"});
			my @tmp2 = split (/\t/, $MN_score_community_coverage{"$gn\t$Path_keys[$j]"});
			$coverage = ($tmp1[2] + $tmp2[2]) / 2;
			$MN_score_community_coverage2{"$gn\t$pair"} = $Bin2Cat{$gn}."\t".$coverage;
		}
	}
}

`mkdir $output/MN-score_result`;
open OUT, ">$output/MN-score_result/MN-score_result_table_input.txt";
print OUT "#Genome\tFunc1\tFunc2\tTaxonomic Group\tCoverage value\(average\)\n";
foreach my $gn_n_pair (sort keys %MN_score_community_coverage2){
	print OUT "$gn_n_pair\t$MN_score_community_coverage2{$gn_n_pair}\n";
}
close OUT;

#read the "MN-score_result_table_input.txt" and make the "MN-score_result.txt", which is the final result of MN-score
my %Input = (); # whole line => [0]  Acidimicrobiia_bacterium_UWMA-0264	[1]  C-S-01:Organic carbon oxidation	[2] C-S-04:Acetate oxidation	[3] Actinobacteriota	[4] 0.038328883
open IN, "$output/MN-score_result/MN-score_result_table_input.txt";
while (<IN>){
	chomp;
	if (!/\#/){
		my @tmp = split (/\t/,$_);
		$Input{$_}[0] = $tmp[0];
		$Input{$_}[1] = $tmp[1];
		$Input{$_}[2] = $tmp[2];
		$Input{$_}[3] = $tmp[3];
		$Input{$_}[4] = $tmp[4];
	}
}
close IN;

my %Output1 = (); # func. => category => summed coverage
my %Output2 = (); # func. =>  summed coverage
my %Cat2 =();
foreach my $key (sort keys %Input){
	$Output1{$Input{$key}[1]}{$Input{$key}[3]} += $Input{$key}[4];
	$Output1{$Input{$key}[2]}{$Input{$key}[3]} += $Input{$key}[4];
	$Cat2{$Input{$key}[3]} = 1;
	$Output2{$Input{$key}[1]} += $Input{$key}[4];
	$Output2{$Input{$key}[2]} += $Input{$key}[4];
}

my %Output3 = (); # the contribution percentage for each function
my $sum_cov_for_output2 = 0;
foreach my $func (sort keys %Output2){
	$sum_cov_for_output2 += $Output2{$func};
}

foreach my $func (sort keys %Output2){
	my $var = ($Output2{$func} / $sum_cov_for_output2) * 100;
	$Output3{$func} = sprintf "%.1f",$var;
}

my %Output4 = (); # func. => category => percentage  
#The func. and each category contribution percentage table
foreach my $func (sort keys %Output1){
	foreach my $cat (sort keys %Cat2){
		my $var = 0;
		if ($Output2{$func} and $Output1{$func}{$cat}){
			$var = ($Output1{$func}{$cat} / $Output2{$func}) * 100;
		}
		$Output4{$func}{$cat} = sprintf "%.1f",$var;
	}
}

open OUT, ">$output/MN-score_result/MN-score_result.txt";
my $row=join("\t", sort keys %Cat2);
print OUT "Function\tMN-score for each function\t$row\n";
foreach my $tmp1 (sort keys %Output4)
{
        print OUT $tmp1."\t";
		print OUT $Output3{$tmp1}."\t";
        my @tmp = ();
        foreach my $tmp2 (sort keys %Cat2)
        {
                if (exists $Output4{$tmp1}{$tmp2})
                {
                        push @tmp, $Output4{$tmp1}{$tmp2};
                }
                else
                {
                        push @tmp,"0"
                }
        }
        print OUT join("\t",@tmp)."\n";
}
close OUT;

$datestring = strftime "%Y-%m-%d %H:%M:%S", localtime; 
print "\[$datestring\] Calculating MN-score is done\n";
$endtime = $datestring;
}

#print information about this run:
open OUT, ">$output/METABOLIC_run.log";
#print information about this run:
print OUT "$version
Run Start: $starttime
Run End: $endtime
Input Reads: $omic_reads_parameters
Reads type: $omic_reads_type
Input Genome directory (nucleotides): $input_genome_folder
Number of Threads: $cpu_numbers
Prodigal Method: $prodigal_method
KOfam DB: $kofam_db_size
Module Cutoff Value: $module_cutoff
Output directory: $output\n";
close OUT;

##subroutines
#input ko_list, return a result hash of threshold and score_type
sub _get_kofam_db_KO_threshold{
	my $list = $_[0]; 
	my $prok_list = "";
	if ($kofam_db_size eq "full"){
		$prok_list = "$_[1]/prokaryote.hal";
	}elsif ($kofam_db_size eq "small"){
		$prok_list = "$_[1]/All_Module_KO_ids.txt";
	}
	my %result = ();	
	open IN, "$list";
	while (<IN>){
		chomp;
		if (/^K/){
			my @tmp = split (/\t/);
			if ($tmp[1] eq "\-"){
				my $hmm_id = "$tmp[0]\.hmm";
				$result{$hmm_id} = "50|full";
			}else{
				my $hmm_id = "$tmp[0]\.hmm";
				$result{$hmm_id} = "$tmp[1]|$tmp[2]";
			}
		}
	}
	close IN;
	
	my %Prok_list = ();
	open IN, "$prok_list";
	while (<IN>){
		chomp;
		$Prok_list{$_} = 1;
	}	
	close IN;
	
	foreach my $key (sort keys %result){
		if (!$Prok_list{$key}){
			delete $result{$key};
		}
	}
	return (%result);
}

#input the hmm_table_temp hash, return a hmm to ko hash (like: TIGR02694.hmm => K08355.hmm)
sub _get_hmm_2_KO_hash{
	my %hash = @_;
	my %result = ();
	my %result2 = ();
	foreach my $line_no (sort keys %hash){
		my @tmp = split (/\t/, $hash{$line_no});
		my $hmm = $tmp[5];
		my $ko = $tmp[6];
		if ($hmm and $hmm !~ /\;/){
				$result{$hmm} = $ko."\.hmm";
		}elsif ($hmm){
			my @array_hmm = split (/\; /, $hmm);
			my @array_ko = split (/\; /, $ko);
			for(my $i=0; $i<=$#array_hmm; $i++){
				$result{$array_hmm[$i]} = $array_ko[$i]."\.hmm";
			}
		}
	}
	
	foreach my $hmm (sort keys %result){
		if ($result{$hmm} =~ /^K\d\d\d\d\d/){
			$result2{$hmm} = $result{$hmm};
		}
	}
	return %result2;
}

#input faa file, and output the seq hash
sub _get_faa_seq{
	my $file = $_[0]; 
	my ($file_name) = $file =~ /^.+\/(.+?)\.faa/;
	my %result = (); my $head = "";
	open _IN, "$file";
	while(<_IN>){
		chomp;
		if (/>/){
			if (/\s/){
				my ($head_old) = $_ =~ /^>(.+?)\s/; 
				$head = ">".$file_name."~~".$head_old;
			}else{
				my ($head_old) = $_ =~ /^>(.+?)$/; 
				$head = ">".$file_name."~~".$head_old;
			}	
			$result{$head} = "";
		}else{
			$result{$head} .= $_; 
			$result{$head} =~ s/\*$//g;
		}
	}
	close _IN;
	return %result;
}

sub _get_Genome_coverge{
	my $reads = $_[0];
	my $folder = $_[1]; # the input_genome_folder
	#cat all the genes
	my %Seq = (); my $head = "";
	open __IN, "ls $folder/*.gene | ";
	while (<__IN>){
		chomp;
		my $file = $_;
		my ($seq_name) = $file =~ /$folder\/(.+?)\.gene/; 
		open __INN, $file;
		while (<__INN>){
			chomp;
			if (/>/){
				my $head_old = $_;
				my $head_old_2;
				if ($head_old =~ /\s/){
					($head_old_2) = $head_old =~ />(.+?)\s/;
				}else{
					($head_old_2) = $head_old =~ />(.+?)$/;
				}
				$head = ">".$seq_name."~~".$head_old_2;
				$Seq{$head} = "";
			}else{
				$Seq{$head} .= $_;
			}
		}
		close __INN;
	}
	close __IN;
	
	open __OUT, ">$output/All_gene_collections.gene";
	foreach my $key (sort keys %Seq){
		print __OUT "$key\n$Seq{$key}\n";
	}
	close __OUT;
	
	system ("bowtie2-build $output/All_gene_collections.gene $output/All_gene_collections.gene.scaffold --quiet");
	my %Reads = (); my $i = 1;
	open __IN, "$reads";
	while (<__IN>){
		chomp;
		if (!/^#/){
			my @tmp = split (/\,/,$_);
			my $tmp_link = "";
			if ($test ne "true"){
				$tmp_link = $tmp[0]."\t".$tmp[1];
			}else{
				$tmp_link = "$METABOLIC_dir/METABOLIC_test_files/METABOLIC_test_reads/".$tmp[0]."\t"."$METABOLIC_dir/METABOLIC_test_files/METABOLIC_test_reads/".$tmp[1];
			}
			$Reads{$tmp_link} = $i;
			$i++;
		}
	}
	close __IN;
	
	my %Read_seq_numbers = (); # read pair => read seq number
	my $average_read_seq_number = 0; # The average read
	open OUT__,">$output/tmp_calculate_depth.sh";
	foreach my $key (sort keys %Reads){
		my $j = $Reads{$key};
		my @tmp = split (/\t/,$key); 
		my $seq_number = `cat $tmp[0] | wc -l`; chomp $seq_number; $Read_seq_numbers{$key} = $seq_number;
		print OUT__ "bowtie2 -x $output/All_gene_collections.gene.scaffold -1 $tmp[0] -2 $tmp[1] -S $output/All_gene_collections_mapped.$j.sam -p $cpu_numbers --quiet;";
		print OUT__ "samtools view -bS $output/All_gene_collections_mapped.$j.sam > $output/All_gene_collections_mapped.$j.bam -@ $cpu_numbers 2> /dev/null;";
		print OUT__ "mkdir $output/sambamba_tmpfiles; sambamba sort  $output/All_gene_collections_mapped.$j.bam --tmpdir $output/sambamba_tmpfiles -o $output/All_gene_collections_mapped.$j.sorted.bam 2> /dev/null;";
		print OUT__ "samtools index $output/All_gene_collections_mapped.$j.sorted.bam > /dev/null;";
		print OUT__ "samtools flagstat $output/All_gene_collections_mapped.$j.sorted.bam > $output/All_gene_collections_mapped.$j.sorted.stat 2> /dev/null\n";
		print OUT__ "rm $output/All_gene_collections_mapped.$j.sam $output/All_gene_collections_mapped.$j.bam;rm -r $output/sambamba_tmpfiles\n";
	}
	close OUT__;
	
	foreach my $key (sort keys %Read_seq_numbers){
		$average_read_seq_number += $Read_seq_numbers{$key};
	}	
	
	my @Read_seq_numbers = keys %Read_seq_numbers;
	$average_read_seq_number = $average_read_seq_number / (scalar @Read_seq_numbers) ;
	
	#parallel run calculate coverage
	my @Runs5; 
	open IN, "$output/tmp_calculate_depth.sh";
	while (<IN>){
		chomp;
		push @Runs5, $_;
	}
	close IN;
	`rm $output/tmp_calculate_depth.sh`;

	my @arrayrefs5 = split_by($cpu_numbers,@Runs5); 

	foreach my $key (@arrayrefs5){
		my $cmd; # to store the bash cmd
		my @tmp_cmd;
			$cmd = join (" &\n",@$key);
		$cmd .= "\nwait\n";
		`$cmd`;
	}
	
	system ("coverm contig --methods metabat --bam-files  $output/All_gene_collections_mapped.*.sorted.bam > $output/All_gene_collections_mapped.depth.txt 2> /dev/null");
	`rm $output/*.bt2;rm $output/All_gene_collections.gene; rm $output/*.bam; rm $output/*.sorted.stat;rm $output/*.bai`;
	
	my %h = (); # average => bin => all gene coverage values
	my @h_head = ();  
	my @h_head_num = (); 
	my %Bin = ();
	
	if ($omic_reads_type eq "metaG"){
		open __IN, "$output/All_gene_collections_mapped.depth.txt";
		while (<__IN>){
			chomp;
			if (/^contigName/){
					my @tmp = split (/\t/);@h_head = @tmp;
					for(my $i=0; $i<=$#h_head; $i++){
							if ($h_head[$i] =~ /^totalAvgDepth$/){
								push @h_head_num, $i;
							}
					}
			}else{
					my @tmp = split (/\t/);
					my ($bin) = $tmp[0] =~ /^(.+?)\~\~/;
					$Bin{$bin} = 1;
					foreach my $i (@h_head_num){
                        if (!exists $h{$h_head[$i]}{$bin}){
                                $h{$h_head[$i]}{$bin} = $tmp[$i];
                        }else{
                                $h{$h_head[$i]}{$bin} .= "\t".$tmp[$i];
                        }
					}
			}
		}
		close __IN;
	}elsif ($omic_reads_type eq "metaT"){
		open __IN, "$output/All_gene_collections_mapped.depth.txt";
		open __OUT, ">$output/All_gene_collections_transcript_coverage.txt"; # contigName => transcript coverage in RPKM
		while (<__IN>){
			chomp;
			if (/^contigName/){
				my @tmp = split (/\t/);@h_head = @tmp;
				for(my $i=0; $i<=$#h_head; $i++){
					if ($h_head[$i] =~ /^totalAvgDepth$/){
						push @h_head_num, $i;
					}	
				}
				print __OUT "contigName\tTranscript coverage in RPKM\n";
			}else{
					my @tmp = split (/\t/);
					my ($bin) = $tmp[0] =~ /^(.+?)\~\~/;
					my $geneLength = $tmp[1] / 1000; # gene length in kb
					$Bin{$bin} = 1;
					foreach my $i (@h_head_num){
                        my $transcript_coverage = $tmp[$i] * (1000000  / $average_read_seq_number) / $geneLength;
						print __OUT "$tmp[0]\t$transcript_coverage\n";
						if (!exists $h{$h_head[$i]}{$bin}){
                                $h{$h_head[$i]}{$bin} = $transcript_coverage;
                        }else{
                                $h{$h_head[$i]}{$bin} .= "\t".$transcript_coverage;
                        }
					}
			}
		}
		close __IN;
		close __OUT;
	}
	
	my %Bin2Cov = (); #bin => cov value
	my $total_cov = 0;
	foreach my $i (@h_head_num){
        foreach my $bin (sort keys %Bin){
                my @tmp = split (/\t/, $h{$h_head[$i]}{$bin});
                my $stat = Statistics::Descriptive::Full->new();
                $stat->add_data(\@tmp);
                my $mean = $stat->mean();
                $h{$h_head[$i]}{$bin} = $mean;
				$Bin2Cov{$bin} = $mean;
				$total_cov += $mean;
        }
	}
	
	my %Bin2cov_percentage = ();
	foreach my $bin (sort keys %Bin2Cov){
		my $percentage = $Bin2Cov{$bin} / $total_cov;
		$Bin2cov_percentage{$bin} = $percentage;
	}
	#system ("rm $output/All_gene_collections*");
	return %Bin2cov_percentage;
}

sub _store_seq{
	my $file = $_[0];
	my %Seq = (); my $head = "";
	open _IN, "$file";
	while (<_IN>){
		chomp;
		if (/>/){
			if (/\s/){
				($head) = $_ =~ /^(>.+?)\s/;
				$Seq{$head} = "";
			}else{
				($head) = $_ =~ /^(>.+?)$/;
				$Seq{$head} = "";
			}
		}else{
			$Seq{$head} .= $_;
		}
	}
	close _IN;
	return %Seq;
}

sub _get_1_from_input_faa{
	my $input_file = $_[0];
	my $seq = $_[1];
	my $output_file = $_[2];
	my %Seq = (); my $head = "";	
	open _IN, "$input_file";
	while (<_IN>){
		chomp;
		if (/>/){
			($head) = $_ =~ /^(>.+?)\s/;
			$Seq{$head} = "";
		}else{
			$Seq{$head} .= $_;
		}
	}
	close _IN;
	
	open OUT, ">$output_file";
	print OUT "$seq\n$Seq{$seq}\n";
	close OUT;	
}

sub _get_check_score{
	my $input_file = $_[0];
	my $score = 0;
	open _IN, "$input_file";
	while (<_IN>){
		chomp;
		if (!/^#/){
			my $line = $_;
			$line =~ s/ +/ /g;
			my @tmp = split (/ /,$line);
			if ($tmp[5] >= $score){
				$score = $tmp[5];
			}
		}
	}
	close _IN;
	return $score;
}

sub _get_motif{
	my $file = $_[0];
	my %Hash; # protein id => motif sequences (dsrC => GPXKXXCXXXGXPXPXXCX)
	open _IN, "$file";
	while (<_IN>){
		chomp;
		my @tmp = split (/\:/);
		$Hash{$tmp[0]} = $tmp[1];
	}
	close _IN;
	return %Hash;
}

sub _get_motif_pair{
	my $file = $_[0];
	my %Hash; # dsrC => tusE
	open _IN, "$file";
	while (<_IN>){
		chomp;
		my @tmp = split (/\:/);
		$Hash{$tmp[0]} = $tmp[1];
	}
	close _IN;
	return %Hash;
} 

sub _run_parallel{
	my $file = $_[0];
	my $cpu_numbers_ = $_[1];
	my @Runs; 
	open ___IN, $file;
	while (<___IN>){
		chomp;
		push @Runs, $_;
	}
	close ___IN;

	my $pm = Parallel::ForkManager->new($cpu_numbers_);
	foreach my $run (@Runs){
		my $pid = $pm->start and next;
		`$run`;
		$pm->finish;
	}
	$pm->wait_all_children;
}
