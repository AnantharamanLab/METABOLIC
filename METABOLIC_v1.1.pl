#!/usr/bin/perl

###########################

# METABOLIC.pl

# METABOLIC  =>  METabolic And BiogeOchemistry anaLyses In miCrobes

# This software gives a metabolic and biogeochemical function trait profile to given genome datasets 
# [either metagenome-assembled genomes (MAGs), single-cell amplified genomes (SAGs) or pure culture sequenced genomes]. 
# It also integrates the genome coverage to make element cycling pathways.

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
use Excel::Writer::XLSX;
use Getopt::Long;
use Statistics::Descriptive;

=head1 DESCRIPTION

        Takes a folder containing genome files to generate a profile of their metablic and biogeochemical functions
        
=head1 USAGE       

        perl METABOLIC.pl -t 40 -m-cutoff 0.75 -w . -in Genome_proteins 
		(When you also want to calculate genome coverages, you would have to add genome files with the same basename and the extention as ".fasta" in this folder)
        
        perl METABOLIC.pl -t 40 -m-cutoff 0.75 -w . -in-gn Genome_files -o omic_reads_parameters.txt
        
=head1 OPTIONS

	-t         or -cpu            [integer] The cpu numbers to run the hmmscan (default: 40)
        -m-cutoff  or -module-cutoff  [float]   The cutoff value to assign the presence of a specific KEGG module (KEGG module step present numbers / KEGG module step total number) (default: 0.75) 
        -w         or -working-dir    [string]  The working directory (default: '.') 
        -in                           [string]  The folder pf given genome faa files [should also give the genome fasta files and genone gene files if the (meta)genome/(meta)transciptome datasets are included]
        -in-gn                        [string]  The folder of given genome fasta files (Prodigal will be used to annotate your genomes)
	-p         or prodigal-method [string]  "meta" or "single" for prodigal to annotate the orf
        -o         or omic-reads      [string]  The file which indicates the address of omic reads
	
=head1 INSTRUCTIONS

	GitHub: https://github.com/AnantharamanLab/METABOLIC

=head1 OUTPUT

        Prodigal-annotated protein and gene files will be given in the input folder.
        The METABOLIC result table will be generated.
	Biogeochemical pathway diagrams will be generated.
		
=head1 COPYRIGHT

	Zhichao Zhou, zczhou2017@gmail.com
	Karthik Anantharaman, karthik@bact.wisc.edu
	Anantharaman Microbiome Laboratory
	Department of Bacteriology, University of Wisconsin, Madison
=cut

#The options 
# option variables with default value
my $cpu_numbers = 40; #parallel running cpu numbers
my $module_cutoff = 0.75; # The cutoff value to assign the existence of a module
my $working_dir = '.' ; # The working directory
my $input_protein_folder;  #input microbial genome protein files
my $input_genome_folder; #input microbial genome fasta files
my $omic_reads_parameters; #The address of omic reads
my $prodigal_method = "meta"; #the prodigal method to annotate orfs
my $version="METABOLIC.pl v1.0";

GetOptions(
	'cpu|t=i' => \$cpu_numbers,
	'module-cutoff|m-cutoff=f' => \$module_cutoff,
	'working-dir|w=s' => \$working_dir,
	'in=s' => \$input_protein_folder,
	'in-gn=s' => \$input_genome_folder,
	'prodigal-method|p=s' => \$prodigal_method,
	'omic-reads|o=s' => \$omic_reads_parameters,
	'help|h' => sub{system('perldoc', $0); exit;},
	'v|version'=>sub{print $version."\n"; exit;}
) or die("Getting options from the command line failed, please check your options");

##Pre-required files and documents
 #METABOLIC hmm database files
 my $METABOLIC_hmm_db_address = "$working_dir/METABOLIC_hmm_db";
 
 #KofamKOALA hmm database files, July, 3 updated 
 #link: ftp://ftp.genome.jp/pub/db/kofam/
 my $kofam_db_address = "$working_dir/kofam_database/profiles";
 my $kofam_db_KO_list = "$working_dir/kofam_database/ko_list";
 
 #input hmm information table as a template
 my $hmm_table_temp = "$working_dir/METABOLIC_temp_and_db/hmm_table_template.txt";
 my $hmm_table_temp_2 = "$working_dir/METABOLIC_temp_and_db/hmm_table_template_2.txt"; 
 
 #the KEGG module information
 my $ko_module_table = "$working_dir/METABOLIC_temp_and_db/ko00002.keg";
 
 #the KEGG module step db 
 my $ko_module_step_db = "$working_dir/METABOLIC_temp_and_db/kegg_module_step_db.txt";
 
 #the pathway information to draw element cycling diagrams
 my $R_pathways = "$working_dir/METABOLIC_temp_and_db/R_pathways.txt";
 
 #the motif files to validate specific protein hits
 my $motif_file = "$working_dir/METABOLIC_temp_and_db/motif.txt";
 my $motif_pair_file = "$working_dir/METABOLIC_temp_and_db/motif.pair.txt";


##Main Body
#The present time
my $datestring = strftime "%Y-%m-%d %H:%M:%S", localtime; 

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

#run hmmscan
my $absolute_path = $ENV{'PWD'}; 
#=pod
my $folder_intermediate_files = "$absolute_path/intermediate_files";
if ($folder_intermediate_files){
	`rm -r intermediate_files`;
	`mkdir intermediate_files`;
}else{
	`mkdir intermediate_files`;
}
#=cut
if ($input_genome_folder){
	open IN, "ls $input_genome_folder/*.fasta |";
	$datestring = strftime "%Y-%m-%d %H:%M:%S", localtime; 
	print "\[$datestring\] The Prodigal annotation is running...\n";
	while (<IN>){
		chomp;
		my ($gn_id) = $_ =~ /^$input_genome_folder\/(.+?)\.fasta/;
		`prodigal -i $input_genome_folder/$gn_id.fasta -a $input_genome_folder/$gn_id.faa -o $input_genome_folder/$gn_id.gff -f gff -p $prodigal_method -q `;
		`perl Accessory_scripts/gff2fasta_mdf.pl -g $input_genome_folder/$gn_id.gff -f $input_genome_folder/$gn_id.fasta -o $input_genome_folder/$gn_id.gene`;
	}
	close IN;
	$input_protein_folder = $input_genome_folder;
	$datestring = strftime "%Y-%m-%d %H:%M:%S", localtime; 
	print "\[$datestring\] The Prodigal annotation is finished\n";
}

my %Total_faa_seq = (); #Store the total faa file into a hash
my %Genome_id = (); # genome id => Gn001 like new id 
open OUT, ">tmp_run_hmmscan.sh";
open IN,"ls $input_protein_folder/*.faa |";
while (<IN>){
	chomp;
	my $file = $_;
	my ($gn_id) = $file =~ /^$input_protein_folder\/(.+?)\.faa/; $Genome_id{$gn_id} = 1; 
	#Store faa file into a hash
	%Total_faa_seq = (%Total_faa_seq, _get_faa_seq($file));
	foreach my $hmm (sort keys %Total_hmm2threshold){
		my ($threshold,$score_type) = $Total_hmm2threshold{$hmm} =~ /^(.+?)\|(.+?)$/;
		#print "$hmm\t$threshold\t$score_type\n";
		if ($score_type eq "full"){
			if ($hmm !~ /K\d\d\d\d\d/){
				print OUT "hmmscan -T $threshold --cpu 1 --tblout intermediate_files/$hmm.$gn_id.hmmscan_result.txt $METABOLIC_hmm_db_address/$hmm $input_protein_folder/$gn_id.faa\n";
			}else{
				print OUT "hmmscan -T $threshold --cpu 1 --tblout intermediate_files/$hmm.$gn_id.hmmscan_result.txt $kofam_db_address/$hmm $input_protein_folder/$gn_id.faa\n";
			}
		}else{
			if ($hmm !~ /K\d\d\d\d\d/){
				print OUT "hmmscan --domT $threshold --cpu 1 --tblout intermediate_files/$hmm.$gn_id.hmmscan_result.txt $METABOLIC_hmm_db_address/$hmm $input_protein_folder/$gn_id.faa\n";
			}else{
				print OUT "hmmscan --domT $threshold --cpu 1 --tblout intermediate_files/$hmm.$gn_id.hmmscan_result.txt $kofam_db_address/$hmm $input_protein_folder/$gn_id.faa\n";
			}
		}
	}
}
close IN;
close OUT;

$datestring = strftime "%Y-%m-%d %H:%M:%S", localtime; 
print "\[$datestring\] The hmmscan is running with $cpu_numbers cpu threads...\n";
#=pod
`cat tmp_run_hmmscan.sh | parallel -j $cpu_numbers --load 80% --noswap`;
`rm tmp_run_hmmscan.sh`;
#=cut
$datestring = strftime "%Y-%m-%d %H:%M:%S", localtime; 
print "\[$datestring\] The hmmscan is finished\n";

#Store motif validation files
my %Motif = _get_motif($motif_file); #protein id => motif sequences (DsrC => GPXKXXCXXXGXPXPXXCX)
my %Motif_pair = _get_motif_pair($motif_pair_file); # dsrC => tusE

#summarize hmmscan result and print table
my %Hmmscan_result = (); # genome_name => hmm => numbers
my %Hmmscan_hits = (); # genome_name => hmm => hits
my %Hmm_id = (); # hmm => 1 
open IN, "find intermediate_files -type f -name '*.hmmscan_result.txt' | ";
while (<IN>){
        chomp;
        my $file_name = $_;
        my ($hmm) = $file_name =~ /^intermediate_files\/(.+?\.hmm)\./; 
		$Hmm_id{$hmm} = 1;
		my ($gn_id) = $file_name =~ /\.hmm\.(.+?)\.hmmscan_result/;
        my @Hits = ();
        open INN, "$file_name";
        while (<INN>){
                chomp;
                if (!/^#/){
                        my $line = $_; $line =~ s/\s+/\t/g;
						my @tmp = split (/\t/,$line);
						my ($threshold,$score_type) = $Total_hmm2threshold{$hmm} =~ /^(.+?)\|(.+?)$/; 
						if ($score_type eq "domain"){
							if ($tmp[8] >= $threshold){
								my ($hmm_basename) = $hmm =~ /^(.+?)\.hmm/; 
								if (exists $Motif{$hmm_basename}){
									#print "$hmm_basename in $gn_id motif validaiton is effective\n";
									my $seq; # the protein seq
									my $motif = $Motif{$hmm_basename};
									my $motif_anti = "111"; # the meaningless numbers that serve as the anti motif
									if (exists $Motif_pair{$hmm_basename}){
										$motif_anti = $Motif{$Motif_pair{$hmm_basename}};
									}
									$motif =~ s/X/\[ARNDCQEGHILKMFPSTWYV\]/g; $motif_anti =~ s/X/\[ARNDCQEGHILKMFPSTWYV\]/g;
									my %Seq_gn = _store_seq("$input_protein_folder/$gn_id.faa"); # get the genome sequences
									$seq = $Seq_gn{">$tmp[2]"};
									if ($seq =~ /$motif/ and $seq !~ /$motif_anti/){
										push @Hits, $tmp[2];
										$Hmmscan_result{$gn_id}{$hmm}++; 
										#print "$hmm_basename in $gn_id has been tested, and it has passed\n";
									}else{
										#print "$hmm_basename in $gn_id has been tested, and it has failed\n";
									}
								}else{
									push @Hits, $tmp[2];
									$Hmmscan_result{$gn_id}{$hmm}++;
								}
							}
						}else{
							my ($hmm_basename) = $hmm =~ /^(.+?)\.hmm/; 
							if (exists $Motif{$hmm_basename}){
								#print "$hmm_basename in $gn_id motif validaiton is effective\n";
								my $seq; # the protein seq
								my $motif = $Motif{$hmm_basename};
								my $motif_anti = "111"; # the meaningless numbers that serve as the anti motif
								if (exists $Motif_pair{$hmm_basename}){
									$motif_anti = $Motif{$Motif_pair{$hmm_basename}};
								}
								$motif =~ s/X/\[ARNDCQEGHILKMFPSTWYV\]/g; $motif_anti =~ s/X/\[ARNDCQEGHILKMFPSTWYV\]/g;
								my %Seq_gn = _store_seq("$input_protein_folder/$gn_id.faa"); # get the genome sequences
								$seq = $Seq_gn{">$tmp[2]"};
								if ($seq =~ /$motif/ and $seq !~ /$motif_anti/){
									push @Hits, $tmp[2];
									$Hmmscan_result{$gn_id}{$hmm}++; 
									#print "$hmm_basename in $gn_id has been tested, and it has passed\n";
								}else{
									#print "$hmm_basename in $gn_id has been tested, and it has failed\n";
								}
							}else{
								push @Hits, $tmp[2];
								$Hmmscan_result{$gn_id}{$hmm}++;
							}
						}
                }
        }
        close INN;
        $Hmmscan_hits{$gn_id}{$hmm} = join("\,",@Hits);
}
close IN;

#print out hmm result excel file
my $hmm_result = Excel::Writer::XLSX->new('METABOLIC_result.xlsx');

# Add worksheets
my $worksheet1 = $hmm_result->add_worksheet("HMMHitNum"); #$worksheet1 indicates the hmm result - hmm presence, hit numbers and hits 
my $worksheet2 = $hmm_result->add_worksheet("FunctionHit"); #$worksheet2 indicates the hmm result - Function presence 
my $worksheet3 = $hmm_result->add_worksheet("KEGGModuleHit"); #$worksheet3 indicates the hmm result - module presence 
my $worksheet4 = $hmm_result->add_worksheet("KEGGModuleStepHit"); #$worksheet4 indicates the hmm result - module step presence 

#  Add and define a format
my $format_mainbody = $hmm_result->add_format(
    color => 'black', align => 'left', align  => 'vcenter', font => 'Arial', size => 10, border => 1, border_color => 'black'
);

my $format_head = $hmm_result->add_format(
    color => 'black', align => 'left', align  => 'vcenter', font => 'Arial', size => 10, border => 1, border_color => 'black', bold => 1
);

my $format_mainbody_green = $hmm_result->add_format(
    color => 'white', align => 'left', align  => 'vcenter', font => 'Arial', size => 10, border => 1, border_color => 'black', bg_color => 'green'
);

my $format_mainbody_red = $hmm_result->add_format(
    color => 'white', align => 'left', align  => 'vcenter', font => 'Arial', size => 10, border => 1, border_color => 'black', bg_color => 'red'
);

$datestring = strftime "%Y-%m-%d %H:%M:%S", localtime; 
print "\[$datestring\] The hmm hit result is calculating...\n";


#pre-set the coordinates
my $col = 0; my $row = 0; 
my $col2 = 0; my $row2 = 0; # for worksheet2
#write the head of hmm result to worksheet1
for(my $i=0; $i<=9; $i++){
	$col = $i;
	$worksheet1->write($row, $col, $Hmm_table_head[($i+1)], $format_head);
}
#write the head of hmm result to worksheet2
for(my $i=0; $i<=3; $i++){
	$col2 = $i;
	$worksheet2->write($row2, $col2, $Hmm_table_head[($i+1)], $format_head);
}

#Make new genome name and print head to worksheet1
my %Genome_id_new = ();
my $i = 1;
foreach my $gn_id (sort keys %Genome_id){	
	my $num = (sprintf "%03d", $i);
	$Genome_id_new{$gn_id} = "Gn".$num; #print "$gn_id\t$Genome_id_new{$gn_id}\n";
	$worksheet1->write($row, ($col++), "Gn$num Hmm presence", $format_head);
	$worksheet1->write($row, ($col++), "Gn$num Hit numbers", $format_head);
	$worksheet1->write($row, ($col++), "Gn$num Hits", $format_head);
	$i++;
}

#write head to worksheet2
foreach my $gn_id (sort keys %Genome_id_new){
	$worksheet2->write($row2, ($col2++), "$Genome_id_new{$gn_id} Hmm presence", $format_head);
}

#write the main body of hmm result to worksheet1
foreach my $line_no (sort keys %Hmm_table_temp){
	$row = $line_no;
	if ($row != 205){
		my @tmp = split(/\t/,$Hmm_table_temp{$line_no});
		my $hmm = $tmp[5];
		for(my $i=0; $i<=9; $i++){
				$col = $i;
				$worksheet1->write($row, $col, $tmp[($i+1)], $format_mainbody);					
		}
		
		foreach my $gn_id (sort keys %Genome_id_new){
			my $hmm_presence = "Absent";
			my $hit_num = 0;
			my @Hits = ();
			if ($hmm and $hmm !~ /\;/){
				if ($Hmmscan_result{$gn_id}{$hmm}){
					$hmm_presence = "Present";
				}
				$hit_num = $Hmmscan_result{$gn_id}{$hmm};
				if ($Hmmscan_hits{$gn_id}{$hmm}){
					push @Hits, $Hmmscan_hits{$gn_id}{$hmm}; 
				}else{
					push @Hits, "None";
				}
			}elsif($hmm){
				my @tmp = split (/\; /,$hmm);
				my $sum = 0; my @array_hit_num = ();
				for(my $i=0; $i<=$#tmp; $i++){
					if ($Hmmscan_result{$gn_id}{$tmp[$i]}){
						$sum += $Hmmscan_result{$gn_id}{$tmp[$i]};
						push @array_hit_num, $Hmmscan_result{$gn_id}{$tmp[$i]};
						push @Hits, $Hmmscan_hits{$gn_id}{$tmp[$i]};
					}else{
						push @array_hit_num, "0";
						push @Hits, "None";
					}
				}
				if ($sum){
					$hmm_presence = "Present";
				}
				$hit_num = join ("\; ",@array_hit_num);
			}
			
			if ($hmm_presence eq "Absent"){
				$worksheet1->write($row, ($col++), $hmm_presence, $format_mainbody_green);
			}else{
				$worksheet1->write($row, ($col++), $hmm_presence, $format_mainbody_red);
			}		
			if($hit_num){
				$worksheet1->write($row, ($col++), $hit_num, $format_mainbody);
			}else{
				$worksheet1->write($row, ($col++), "0", $format_mainbody);
			}
			$worksheet1->write($row, ($col++), join("\; ",@Hits), $format_mainbody);
		}
	}
}

#write the main body of hmm result to worksheet2
foreach my $line_no (sort keys %Hmm_table_temp_2){
	$row2 = $line_no;
	my @tmp_table_2 = split(/\t/,$Hmm_table_temp_2{$line_no});
	for(my $i=0; $i<=3; $i++){
		$col2 = $i;
		$worksheet2->write($row2, $col2, $tmp_table_2[($i+2)], $format_mainbody);					
	}
	my $line_no_4_table_1 = $tmp_table_2[1];
	
	if ($line_no_4_table_1 !~ /\|\|/){
		my @tmp_table_1 = split(/\t/,$Hmm_table_temp{$line_no});
		my $hmm = $tmp_table_1[5];	
		foreach my $gn_id (sort keys %Genome_id_new){
			my $hmm_presence = "Absent";
			if ($hmm and $hmm !~ /\;/){
				if ($Hmmscan_result{$gn_id}{$hmm}){
					$hmm_presence = "Present";
				}
			}elsif($hmm){
				my @tmp = split (/\; /,$hmm);
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
			if ($hmm_presence eq "Absent"){
				$worksheet2->write($row2, ($col2++), $hmm_presence, $format_mainbody_green);
			}else{
				$worksheet2->write($row2, ($col2++), $hmm_presence, $format_mainbody_red);
			}	
		}		
	}else{
		my @array_line_no_4_table_1 = split (/\|\|/, $line_no_4_table_1);
		my @array_hmm = ();
		foreach my $line_no_4_table_1 (@array_line_no_4_table_1){
			my @tmp_table_1 = split(/\t/,$Hmm_table_temp{$line_no_4_table_1});
			my $hmm = $tmp_table_1[5];
			if ($hmm and $hmm !~ /\;/){
				push @array_hmm, $hmm;
			}elsif($hmm){
				my @tmp = split (/\; /,$hmm);
				foreach my $key (@tmp){
					push @array_hmm, $key;
				}
			}
		}
		foreach my $gn_id (sort keys %Genome_id_new){
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
			
			if ($hmm_presence eq "Absent"){
				$worksheet2->write($row2, ($col2++), $hmm_presence, $format_mainbody_green);
			}else{
				$worksheet2->write($row2, ($col2++), $hmm_presence, $format_mainbody_red);
			}
		}
	
	}
}
#add genome map to worksheet1
$row = 246;
$worksheet1->write($row, 0, "Genome ID map", $format_head); 
$row++;
foreach my $gn_id (sort keys %Genome_id_new){
	$worksheet1->write($row, 0, $gn_id, $format_mainbody);
	$worksheet1->write($row, 1, $Genome_id_new{$gn_id}, $format_mainbody);
	$row++;
}
$worksheet1->set_column( 'A:E', 25 );

#add genome map to worksheet2
$row2 = 112;
$worksheet2->write($row2, 0, "Genome ID map", $format_head);
$row2++;
foreach my $gn_id (sort keys %Genome_id_new){
	$worksheet2->write($row2, 0, $gn_id, $format_mainbody);
	$worksheet2->write($row2, 1, $Genome_id_new{$gn_id}, $format_mainbody);
	$row2++;
}
$worksheet2->set_column( 'A:C', 25 );

#print out each hmm faa collection
$datestring = strftime "%Y-%m-%d %H:%M:%S", localtime; 
print "\[$datestring\] Generating each hmm faa collection...\n";
#=pod
my $folder_Each_hmm_faa = "$absolute_path/Each_hmm_faa";
if ($folder_Each_hmm_faa){
	`rm -r Each_hmm_faa`;
	`mkdir Each_hmm_faa`;
}else{
	`mkdir Each_hmm_faa`;
}
#=cut

mkdir "$working_dir/Each_hmm_faa";
foreach my $hmm (sort keys %Hmm_id){
	my %Hmm_faa_seq = (); #Store the faa seqs in a hmm
	foreach my $gn_id (sort keys %Hmmscan_hits){
		if ($Hmmscan_hits{$gn_id}{$hmm}){			
			my @Hits = split (/\,/,$Hmmscan_hits{$gn_id}{$hmm});
			foreach my $hit (@Hits){
				my $seq_head = ">".$gn_id."~~".$hit;
				if (exists $Total_faa_seq{$seq_head}){
					$Hmm_faa_seq{$seq_head} = $Total_faa_seq{$seq_head};
				}
			}

		}
	}
	if (%Hmm_faa_seq){
		open OUT, ">$working_dir/Each_hmm_faa/$hmm.collection.faa";
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
	foreach my $gn_id (sort keys %Genome_id_new){
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
	foreach my $gn_id (sort keys %Genome_id_new){
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

#pre-set the coordinates
my $col3 = 0; my $row3 = 0; # for worksheet3
my $col4 = 0; my $row4 = 0; # for worksheet4

#write the head of hmm result to worksheet3
$worksheet3->write(0, 0, "Module ID", $format_head);
$worksheet3->write(0, 1, "Module ", $format_head);
$worksheet3->write(0, 2, "Module Category", $format_head);
$col3 = 3;
foreach my $gn_id (sort keys %Genome_id_new){
	$worksheet3->write($row3, ($col3++), "$Genome_id_new{$gn_id} Hit", $format_head);
}
#write the head of hmm result to worksheet4
$worksheet4->write(0, 0, "Module step", $format_head);
$worksheet4->write(0, 1, "Module ", $format_head);
$worksheet4->write(0, 2, "KO id", $format_head);
$worksheet4->write(0, 3, "Module Category", $format_head);
$col4 = 4;
foreach my $gn_id (sort keys %Genome_id_new){
	$worksheet4->write($row4, ($col4++), "$Genome_id_new{$gn_id} Hit", $format_head);
}

$col3 = 0; $row3 = 1; # for worksheet3
$col4 = 0; $row4 = 1; # for worksheet4

#print the worksheet3 result
foreach my $module (sort keys %Module_result){
	$worksheet3->write($row3, ($col3++), $module, $format_mainbody); 
	$worksheet3->write($row3, ($col3++), $KEGG_module2name{$module}, $format_mainbody); 
	my $cat_4_module = ""; # the category name for module
	foreach my $cat (sort keys %Cat2module){
		if ($Cat2module{$cat} =~ /$module/){
			$cat_4_module = $cat;
		}
	}
	$worksheet3->write($row3, ($col3++), $cat_4_module, $format_mainbody); 
	foreach my $gn_id (sort keys %Genome_id_new){
		if ($Module_result{$module}{$gn_id} eq "Absent"){
			$worksheet3->write($row3, ($col3++), $Module_result{$module}{$gn_id}, $format_mainbody_green);
		}else{
			$worksheet3->write($row3, ($col3++), $Module_result{$module}{$gn_id}, $format_mainbody_red);
		}
	}
	$col3 = 0; $row3++;
	
}
$worksheet3->set_column( 'A:A', 10 );
$worksheet3->set_column( 'B:C', 35 );

#print the worksheet4 result
foreach my $module_step (sort keys %Module_step_result){
	$worksheet4->write($row4, ($col4++), $module_step, $format_mainbody); 
	my ($module) = $module_step =~ /^(M.+?)\+/;
	$worksheet4->write($row4, ($col4++), $KEGG_module2name{$module}, $format_mainbody); 
	$worksheet4->write($row4, ($col4++), $KEGG_module{$module_step}[2], $format_mainbody); 
	my $cat_4_module = ""; # the category name for module
	foreach my $cat (sort keys %Cat2module){
		if ($Cat2module{$cat} =~ /$module/){
			$cat_4_module = $cat;
		}
	}
	$worksheet4->write($row4, ($col4++), $cat_4_module, $format_mainbody); 
	foreach my $gn_id (sort keys %Genome_id_new){
		if (!$Module_step_result{$module_step}{$gn_id}){
			$worksheet4->write($row4, ($col4++), "Absent", $format_mainbody_green);
		}else{
			$worksheet4->write($row4, ($col4++), "Present", $format_mainbody_red);
		}
	}
	$col4 = 0; $row4++;
	
}
$worksheet4->set_column( 'A:A', 10 );
$worksheet4->set_column( 'B:D', 35 );

$hmm_result->close();

$datestring = strftime "%Y-%m-%d %H:%M:%S", localtime; 
print "\[$datestring\] METABOLIC table has been generated\n";

#Draw element cycling diagrams
$datestring = strftime "%Y-%m-%d %H:%M:%S", localtime; 
print "\[$datestring\] Drawing element cycling diagrams\n";

#Store R pathways
my %R_pathways = (); #step => KOs
my %R_hmm_ids = ();
open IN, "$R_pathways";
while (<IN>){
	chomp;
	my @tmp = split (/\t/);
	$R_pathways{$tmp[0]} = $tmp[1];
	my @tmp2 = split (/\,/,$tmp[1]);
	foreach my $key (@tmp2){
		$R_hmm_ids{$key} = 1;
	}
}
close IN;

`mkdir R_input`;
#Get each R pathway input files
my %Total_R_input = (); #pathway => gn => 1 or 0
foreach my $gn (sort keys %Hmmscan_result){
	my %R_input = (); #for each input file
	foreach my $key (sort keys %R_pathways){
		$R_input{$key} = 0; $Total_R_input{$key}{$gn} = 0;
		my $hmms = $R_pathways{$key};
		foreach my $hmm_id (sort keys %R_hmm_ids){
			if ($hmms =~ /$hmm_id/ and $Hmmscan_result{$gn}{$hmm_id}){
				$R_input{$key} = 1; $Total_R_input{$key}{$gn} = 1;
			}
		}
	}
	open OUT, ">R_input/$gn.R_input.txt";
	foreach my $key (sort keys %R_input){
		print OUT "$key\t$R_input{$key}\n";
	}	
	close OUT;
}

#The genome coverage: genome id => coverage value
if ($omic_reads_parameters){
	my %Genome_cov = _get_Genome_coverge($omic_reads_parameters,$input_genome_folder);
	
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

	open OUT, ">R_input/Total.R_input.txt"; 
	foreach my $key (sort keys %Total_R_input_2){
		print OUT "$key\t$Total_R_input_2{$key}\n";
	}
	close OUT;
}

`Rscript draw_biogeochemical_cycles.R R_input Output TURE > /dev/null`;

`mv Output/draw_biogeochem_cycles R_output; rm -r Output`;


$datestring = strftime "%Y-%m-%d %H:%M:%S", localtime; 
print "\[$datestring\] Drawing element cycling diagrams finished\n";

##subroutines

#input ko_list, return a result hash of threshold and score_type
sub _get_kofam_db_KO_threshold{
	my $list = $_[0]; my $prok_list = "$_[1]/prokaryote.hal"; 	
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
			my ($head_old) = $_ =~ /^>(.+?)\s/; 
			$head = ">".$file_name."~~".$head_old;
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
	
	open __OUT, ">All_gene_collections.gene";
	foreach my $key (sort keys %Seq){
		print __OUT "$key\n$Seq{$key}\n";
	}
	close __OUT;
	
	system ("bowtie2-build All_gene_collections.gene All_gene_collections.gene.scaffold --quiet");
	my %Reads = ();
	open __IN, "$reads";
	while (<__IN>){
		chomp;
		if (!/^#/){
			my @tmp = split (/\,/,$_);
			$Reads{$tmp[0]} = $tmp[1];
		}
	}
	close __IN;
	my @Reads1; my @Reads2; 
	foreach my $key (sort keys %Reads){
		push @Reads1,$key;
		push @Reads2,$Reads{$key};
	}
	my $reads1 = join(',',@Reads1);
	my $reads2 = join(',',@Reads2);
	system("bowtie2 -x All_gene_collections.gene.scaffold -1 $reads1 -2 $reads2 -S All_gene_collections_mapped.sam -p $cpu_numbers --quiet");
	system("samtools view -bS All_gene_collections_mapped.sam > All_gene_collections_mapped.bam -@ $cpu_numbers 2> /dev/null");
	system("bamtools sort -in All_gene_collections_mapped.bam -out All_gene_collections_mapped.sorted.bam -mem 500 > /dev/null");
	system("samtools index All_gene_collections_mapped.sorted.bam > /dev/null");
	system("samtools flagstat All_gene_collections_mapped.sorted.bam > All_gene_collections_mapped.sorted.stat 2> /dev/null");
	system("rm All_gene_collections_mapped.sam All_gene_collections_mapped.bam ");
	system ("coverm contig --methods metabat --bam-files  All_gene_collections_mapped.sorted.bam > All_gene_collections_mapped.depth.txt 2> /dev/null");
	
	my %h = (); # bam => bin => all gene coverage values
	my @h_head = ();  
	my @h_head_num = (); 
	my %Bin = ();
	open __IN, "All_gene_collections_mapped.depth.txt";
	while (<__IN>){
        chomp;
        if (/^contigName/){
                my @tmp = split (/\t/);@h_head = @tmp;
                for(my $i=0; $i<=$#h_head; $i++){
                        if ($h_head[$i] =~ /.sorted.bam$/){
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
	system ("rm All_gene_collections_mapped.depth.txt");
	
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
	system ("rm All_gene_collections*");
	return %Bin2cov_percentage;
}

sub _store_seq{
	my $file = $_[0];
	my %Seq = (); my $head = "";
	open _IN, "$file";
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
	return %Seq;
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
