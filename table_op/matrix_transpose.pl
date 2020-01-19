#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
my %opts;
GetOptions (\%opts,"i=s",  "o=s" , "html=s");
my $usage = <<"USAGE";

       Usage :perl $0 [options]
                   -i  input file
                   -o  output file
                 
USAGE
die $usage if ( !($opts{i} && $opts{o}) );

my $input = $opts{i};
my $output = $opts{o};

open IN,"$input";

my @matrix;
my $header = <IN>;
chomp $header;
my @s = split /\t/, $header;
my $col_num = @s;
push @matrix, [@s];

my $lines = 2;
while (<IN>){
	chomp;
	@s = split /\t/, $_;

	my $num = @s;
	die "format error: line $lines elements number is not equal to the header\n" if @s != $col_num;
	
	push @matrix, [@s];
	$lines++;
}
close IN;

my @tr_matrix;

for(my $j = 0; $j < $col_num; $j++){
	for(my $i = 0; $i <= $#matrix; $i++){
		$tr_matrix[$j][$i] = $matrix[$i][$j];
	}
}

open OUT,">$output";
for (my $i = 0; $i <= $#tr_matrix; $i++){
	my $l = join "\t", @{$tr_matrix[$i]};
	print OUT "$l\n";
}

close OUT;

#get_html_table($output, $opts{html}, "y");
get_html_link_v2($opts{html}, "Matrix Transpose Results", "$output", "tabH");




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



