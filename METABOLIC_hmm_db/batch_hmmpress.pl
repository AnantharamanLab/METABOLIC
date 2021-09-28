#!/usr/bin/perl

use strict;
use warnings;

open IN, "ls *.hmm |";
while (<IN>){
	chomp;
	`hmmpress $_`;
}
close IN;