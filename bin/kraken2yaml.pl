#!/usr/bin/env perl

use strict;
use Getopt::Long;
use Cwd;
my $usage = qq{
perl my_script.pl
  Getting help:
    [--help]


  Ouput:    
    [--outfile filename]
        The name of the output file. By default the output is the
        standard output
};

#my $outfile = "summary_mqc.yaml";
my $outfile = undef;
my $min_cov = 1;

my $help;

GetOptions(
    "help" => \$help,
    "min_cov=i" => \$min_cov,
    "outfile=s" => \$outfile);

# Print Help and exit
if ($help) {
    print $usage;
    exit(0);
}

if ($outfile) {
    open(STDOUT, ">$outfile") or die("Cannot open $outfile");
}

my %data;
my $dir = getcwd;

my $header = qq(
id: 'kraken2_reports'
section_name: 'Kraken2 Species Hits'
plot_type: 'html'
description: 'list abundant species identified in each sample. S: Species, S1: Subspecies (included in S)'
data: |\n  <dl class="dl-horizontal">
);

printf $header . "\n";


foreach my $file (glob("$dir/*.txt")) {

	my $fh = IO::File->new();
	$fh->open( $file );

	my $f = (split "/" , $file)[-1];
	my $lib = (split ".kraken" , $f)[0];
	my $entry = "<dt>Library</dt><dd><samp>$lib</samp></dd>" ;

        printf "    $entry\n";

	foreach my $line (<$fh>) {

		chomp($line);
		my ($perc,$num_frag_tax,$num_frag_ass,$rank,$tax_id,$name) = split(/\t+/,$line);

		next unless ($rank =~ /S.*/);

		$name =~ s/^\s+|\s+$//g ;
                $rank =~ s/^\s+|\s+$//g ;
                my $tab = "";
                if ($rank eq "S1") {
                	$tab = " ";
                }
		
		if ($name =~ "Severe acute respiratory syndrome coronavirus 2") {
				
			my $entry = "<dt>$perc</dt><dd><samp>$tab$rank:$name</samp></dd>" ;
                        printf "    $entry\n";
	
		} else {

			next if ($perc < $min_cov);

			my $entry = "<dt>$perc</dt><dd><samp>$tab$rank:$name</samp></dd>" ;
	        	printf "    $entry\n";
		}

	}

	close($fh);

}

printf "  </dl>\n";

