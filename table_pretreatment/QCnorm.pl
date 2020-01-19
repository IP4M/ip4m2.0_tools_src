#!/usr/bin/perl
use strict;
use warnings;

use File::Copy;
use File::Path;
use File::Basename;
use FindBin qw($Bin);
use Getopt::Long;
use Scalar::Util qw(looks_like_number);

die "!Too less arguments\n" if @ARGV < 2;

my %opts;
GetOptions (\%opts,"i=s","o=s","g=s",
);

my $html_file = $ARGV[0];
my $input = $opts{i};
my $output = $opts{o};
my $sampe2QC = $opts{g};


sub deal_NA_data{
    my $in = shift;
    my $out = shift;
    open IN,"$in";
    open OUT,">$out";
    my $h = <IN>;
    chomp $h;
    print OUT "$h\n";
    while (<IN>){
        chomp;
        my @s = split /\t/, $_;
        for(my $i = 1; $i <= $#s; $i++){
            if( check_NA( $s[$i] )){
                $s[$i] = 0;
            }
        }
        my $l = join "\t", @s;
        print OUT "$l\n";
    }
    close IN;
    close OUT;  
    
}

sub check_NA{
    my $in = shift;
    if (defined($in) && looks_like_number $in) {
        return 0;
    }else{
        return 1;   
    }
}
deal_NA_data( $input, "__noNA_data__" );
$input = "__noNA_data__";



&check_data_v2( $input, $sampe2QC );

my %name2index;
open IN,"$input";
my $head = <IN>;
chomp $head;
my @s = split /\t/, $head;
for(my $i = 1; $i <= $#s; $i++){
    $name2index{$s[$i]} = $i;
}
close IN;

my @pairs;
open IN,"$sampe2QC";
while(<IN>){
    chomp;
    my @s = split /\t/, $_;
    push @pairs, [@s];
}
close IN;



open IN,"$input";
open OUT,">$output";
$head = <IN>;
print  OUT "$head";
while(<IN>){
    chomp;
    my @s =split /\t/, $_;
    for(my $i = 0; $i <= $#pairs; $i++) {
        my $samp = $pairs[$i][0];
        my $qc = $pairs[$i][1];
        my $samp_index = $name2index{$samp};
        my $qc_index =  $name2index{$qc};
        $s[$samp_index] = ($s[$samp_index] /  $s[$qc_index]) * 10000;
    }
    for(my $i = 0; $i <= $#pairs; $i++) {
        my $qc = $pairs[$i][1];
        my $qc_index =  $name2index{$qc};
        $s[$qc_index] =  10000;
    }
    my $tmp = join "\t", @s;
    print  OUT "$tmp\n";
}

close IN;
close OUT;






get_html_link_v2($html_file, "Normalization results based on QC(pooled samples) Results", "$output", "tabH");

sub findParam{
	my $opt = shift;
	for(my $i = 0; $i <= $#ARGV; $i++){
		if ($ARGV[$i] eq $opt) {
			return $ARGV[$i+1];
		}
	}
}






use File::Basename;
use File::Copy;
#tabH, tabN, txt
sub get_html_link_v2{
	my @param = @_;
	my $html = shift @param;
	my $title = shift @param;
	my @par = @param;

    my $out_dir = "$html.files";
    mkdir $out_dir;


	open HTML, ">$html";
	print HTML "<html><head><title>$title</title></head><body><h3>Output Files:</h3><p><ul>\n";

    for(my $i = 0; $i <= $#par; $i+=2){
            my $file = $par[$i];
            my $type = $par[$i+1];
            my $name = basename $file;
            if( $type eq "tabH" ){
                get_html_table($file, "$out_dir/$name.html", "y");
            }elsif(  $type eq "tabN" ){
                get_html_table($file, "$out_dir/$name.html", "n");
            }elsif(  $type eq "txt" ){
                txt2html($file, "$out_dir/$name.html");
            }else{
                copy($file, "$out_dir/$name.html");
            }
    }

    my $basename = basename $out_dir;
	for (my $i = 0; $i <= $#par; $i+=2){
            my $file = $par[$i];
            my $name = basename $file;
	        print HTML "<li><a href=\"$basename/$name.html\" target=\"_parent\">$name</a></li>\n";
	}
	print HTML "</ul></p>\n";


}



sub get_html_table{
	my $in = shift;
	my $out = shift;
	my $isHeader = shift;
	open IN,"$in";
	open OUT,">$out";
    print OUT "<!DOCTYPE html>\n";
    print OUT "<html>\n";
	print OUT "<table border=\"1\">\n";
	if ( $isHeader eq "y") {
		my $h = <IN>;
		chomp $h;
		my @s = split /\t/, $h;
		print OUT "<tr align=\"center\">";
		for(@s){
			print OUT "<th>$_</th>";
		}
		print OUT "</tr>"
	}

	while (<IN>) {
		chomp;
		my @s = split /\t/, $_;
		print OUT "<tr align=\"center\">";
		for(@s){
			print OUT "<td>$_</td>";
		}
		print OUT "</tr>"
	}

	print OUT "</table>";
    print OUT "</html>";
}




sub txt2html{
	my $txt = shift;
	my $html = shift;
	open IN,"$txt";
	open OUT,">$html";
    print OUT "<!DOCTYPE html>\n";
    print OUT "<html>";
	print OUT "<body>\n";
	while (<IN>) {
		chomp;
		print OUT "$_<br />\n";
	}

	print OUT "</body>\n";
    print OUT "</html>\n";
	close IN;
	close OUT;
}



##check matrix and sample group data

sub check_data{
	my $matrix = shift;
	my $sample_group = shift;
	open IN,"$matrix";
	my $header = <IN>;
	chomp $header;
	my @s = split /\t/, $header;
        my %hs; 
        foreach ( @s ){ 
                $hs{$_} = 1;
        }
	my $col_num = @s;
	my $row_num = 2;
	while (<IN>){
        	chomp;
	        my @s = split /\t/, $_;
	        die "Error: col number in row $row_num is not equal to header!\n" if @s != $col_num;
	        for (my $i = 1; $i <= $#s; $i++){
	                my $col = $i + 1;
	                die "Error: element in row $row_num col $col is not numeric!\n" unless isnumeric( $s[$i] );
	        }
	        $row_num++;
	}
	close IN;
	if ( defined $sample_group ){
		my %cnt;
	        open IN,"$sample_group";
	        while (<IN>){
			chomp;
        	        my @s = split /\t/, $_;
			die "Error: sample group file each line must be 2 columns!\n" if @s != 2;
	                $cnt{$s[0]} += 1;
	                die "Error: sample name: $s[0] in sample group file is not exist in the matrix file!\n" unless defined($hs{$s[0]});
	        }
	        close IN;

		foreach ( keys %cnt){
		        if ( $cnt{$_} > 1){
		                die "Error: sample name: $_ is not uniq in the sample group file!\n";
		        }
		}
	}
}

sub isnumeric {  
        my $val = shift;  
        looks_like_number $val;
}  




sub check_data_v2{
	my $matrix = shift;
	my $sample_group = shift;
	open IN,"$matrix";
	my $header = <IN>;
	chomp $header;
	my @s = split /\t/, $header;
        my %hs; 
        foreach ( @s ){ 
                $hs{$_} = 1;
        }
	my $col_num = @s;
	my $row_num = 2;
	while (<IN>){
        	chomp;
	        my @s = split /\t/, $_;
	        die "Error: col number in row $row_num is not equal to header!\n" if @s != $col_num;
	        for (my $i = 1; $i <= $#s; $i++){
	                my $col = $i + 1;
	                die "Error: element in row $row_num col $col is not numeric!\n" unless isnumeric( $s[$i] );
	        }
	        $row_num++;
	}
	close IN;
	if ( defined $sample_group ){
		    my %cnt;
	        open IN,"$sample_group";
            my %gp;
	        while (<IN>){
			chomp;
        	        my @s = split /\t/, $_;
			die "Error: sample to QC design file each line must be 2 columns!\n" if @s != 2;
	                $cnt{$s[0]} += 1;
                    $gp{$s[1]} += 1;
	                die "Error: sample name: $s[0] in sample to QC design file is not exist in the peak table file!\n" unless defined($hs{$s[0]});
                    die "Error: QC name: $s[1] in sample to QC design file is not exist in the peak table file!\n" unless defined($hs{$s[1]});
	        }
	        close IN;

		foreach ( keys %cnt){
		        if ( $cnt{$_} > 1){
		                die "Error: sample name: $_ is not uniq in the sample group file!\n";
		        }
		}

		foreach ( keys %gp){
		        if ( defined($cnt{$_}) ){
		                die "Error: $_ both existed in sample and QC columns!\n";
		        }
		}


	}
}