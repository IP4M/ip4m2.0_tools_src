#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Getopt::Long qw(:config no_ignore_case bundling);
use Cwd;
use File::Basename;
use Data::Dumper;


my $html = pop @ARGV;
my $output = pop @ARGV;
my @matrices = @ARGV;

unless (scalar @matrices > 1) {
    die "\n\n\tusage: $0 matrixA matrixB ...\n\n";
}

my %matrix;
my %genes;

main: {
    open OUT,">$output";

    foreach my $matrix (@matrices) {

        &parse_matrix($matrix);
        
    }

    ## output new matrix:

    my @colnames = sort keys %matrix;
    print OUT "\t" . join("\t", @colnames) . "\n";

    foreach my $gene (sort keys %genes) {
        
        print OUT "$gene";
        foreach my $colname (@colnames) {
            my $val = $matrix{$colname}->{$gene};
            my $val_;
            unless (defined $val) {
                $val_ = "0";
            }else{
                my @aa  = @{$val};
                my $sum = 0;
                my $num = @aa;
                foreach(@aa){
                    $sum += $_;
                }
                $val_ = $sum/ $num;
            }
            print OUT "\t$val_";
        }
        print OUT "\n";
    }
    
#get_html_table($output, $html, "y");
get_html_link_v2($html, "Tables merged Results", "$output", "tabH");
    close OUT;
    exit(0);

}




####
sub parse_matrix {
    my ($matrix_file) = @_;
    
    open (my $fh, $matrix_file);
    my $header = <$fh>;
    chomp $header;
    my @pos_to_col = split(/\t/, $header);
    my $check_column_ordering_flag = 0;
    
    #foreach my $sample (@pos_to_col) {
    #    if (exists $matrix{$sample}) {
    #        die "Error, already encountered column header: $sample, cannot have redundant column names across matrices.";
    #    }
    #}

    
    while (<$fh>) {
        chomp;
        my @x = split(/\t/);
        
        unless ($check_column_ordering_flag) {
            if (scalar(@x) == scalar(@pos_to_col) + 1) {
                ## header is offset, as is acceptable by R
                ## not acceptable here.  fix it:
                unshift (@pos_to_col, "");
            }
            $check_column_ordering_flag = 1;
          
        }
        

        my $gene = $x[0];
        $genes{$gene} = 1;
        
        for (my $i = 1; $i <= $#x; $i++) {
            my $col = $pos_to_col[$i];
            my $val = $x[$i];
                
            #$matrix{$col}->{$gene} = $val;
            push @{$matrix{$col}->{$gene}}, $val;
            
        }
        
    }
   

    return %matrix;

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