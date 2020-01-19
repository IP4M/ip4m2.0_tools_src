#!/usr/bin/perl -w
use strict;
use warnings;

use Getopt::Long;
use Scalar::Util qw(looks_like_number);
use File::Copy;
use File::Copy::Recursive qw(dircopy);
use File::Basename;
use FindBin qw($Bin);
my $r_path = "$Bin/../../bin/R/bin";
$ENV{'PATH'} = "$r_path;$ENV{'PATH'}";


my %opts;
GetOptions (\%opts, "i=s",    "html=s", 
);



my $input = $opts{i};
my $html = $opts{html};
my $db = "$Bin/hmdb_info.txt";

sub get_greek_hs {
    my %hs;
    open IN,"$Bin/replace_info.txt";
    while(<IN>){
        chomp;
        my @s = split /\t/, $_;
        $hs{$s[0]} = $s[1];
    }
    close IN;
    return %hs;
}

my %gks = get_greek_hs();
my @greeks = keys %gks;

sub process_greek {
    my $string = shift;
    foreach (@greeks) {
        $string=~s/$_/$gks{$_}/g;
    }
    return $string;
}

sub process_name {
    my $string = shift;
    $string =~ s/^\s+//;
    $string =~ s/\s+$//;
    $string =~ s/\s+/ /g;
    $string =~ s/_/ /g;
    $string =~ s/\?//g;
    $string = process_greek($string);
    $string = lc($string);
    return $string;
}


my %database;
open IN, "$db";
<IN>;
while(<IN>){
    chomp;
    my @s = split /\t/, $_;
    for(my $i = 0; $i <= 11; $i++){
        unless(defined($s[$i])){
            $s[$i] = "";
        }
    }
    my $hit = join "\t", @s[1,0,9,8,4,5,7, 11];
    foreach( @s[1, 2, 3] ){
        $_ = lc($_);
        $database{$_} = $hit;
    }
    my @ss = split / \|\| /, $s[10];
    foreach( @ss ){
       $_ = lc($_);
       $database{$_} = $hit;
    }
}
close IN;


my $html_dir = dirname($html);
my $output = "$html_dir/compounds_idmapping.txt";
open IN,"$input";
open OUT, ">$output";
print OUT "Query_Name\tMatched_Name\tHMDB_ID\tPubChem_Compound\tKEGG_Compound_ID\tChemical_Formula\tAverage_Molecular_Weight\tSuper_Class\tPathways\n";
while (<IN>) {
    chomp;
    my $name = process_name($_);
    if(defined($database{$name})){
        my $line = $database{$name};
        print OUT "$_\t$line\n";
    }else{
        print OUT "$_\tNo Hit\n";
    }
}
close IN;
close OUT;


get_html_link_v2($html, "Compounds ID Mapping Results",  "$output", "tabH", );















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
                get_html_table_v2($file, "$out_dir/$name.html", "y");
            }elsif(  $type eq "tabN" ){
                get_html_table_v2($file, "$out_dir/$name.html", "n");
            }elsif(  $type eq "txt" ){
                txt2html($file, "$out_dir/$name.html");
            }elsif(  $type eq "relation" ){
                get_relation_table($file, "$out_dir/$name.html", "y");}
            else{
                copy($file, "$out_dir/$name");
            }
    }

    my $basename = basename $out_dir;
	for (my $i = 0; $i <= $#par; $i+=2){
            my $file = $par[$i];
            my $type = $par[$i+1];
            my $name = basename $file;
            if( $type eq "tabH" ){
                print HTML "<li><a href=\"$basename/$name.html\" target=\"_parent\">$name</a></li>\n";
            }elsif(  $type eq "tabN" ){
                print HTML "<li><a href=\"$basename/$name.html\" target=\"_parent\">$name</a></li>\n";
            }elsif(  $type eq "txt" ){
                print HTML "<li><a href=\"$basename/$name.html\" target=\"_parent\">$name</a></li>\n";
            }elsif(  $type eq "relation" ){
                print HTML "<li><a href=\"$basename/$name.html\" target=\"_parent\">$name</a></li>\n";
            }
            else{
                print HTML "<li><a href=\"$basename/$name\" target=\"_parent\">$name</a></li>\n";
            }
	        
	}
	print HTML "</ul></p>\n";


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



sub get_relation_table{
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
        my $first = shift @s;
        print OUT "<td><b>$first</b></td>";
		for(@s){
			print OUT "<td>$_</td>";
		}
		print OUT "</tr>"
	}

	print OUT "</table>";
    print OUT "</html>";
}



sub get_html_table_v2{
	my $in = shift;
	my $out = shift;
	my $isHeader = shift;
	open IN,"$in";
	open OUT,">$out";
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
        for(my $i = 0; $i<=$#s; $i++){
            if($i == 2){
                if($s[$i] eq ""){
                    print OUT "<td>$s[$i]</td>";
                }else{
                print OUT "<td><a href=\"http://www.hmdb.ca/metabolites/$s[$i]\" target=\"_blank\">$s[$i]<a></td>";}
            }elsif($i == 4){
                if($s[$i] eq ""){
                    print OUT "<td>$s[$i]</td>";
                }else{
                print OUT "<td><a href=\"http://www.kegg.jp/dbget-bin/www_bget?cpd:$s[$i]\" target=\"_blank\">$s[$i]<a></td>";}
            }elsif($i == 3){
                if($s[$i] eq ""){
                    print OUT "<td>$s[$i]</td>";
                }else{
                print OUT "<td><a href=\"https://pubchem.ncbi.nlm.nih.gov/compound/$s[$i]\" target=\"_blank\">$s[$i]<a></td>";}
            }else{
                print OUT "<td>$s[$i]</td>";
            }
        }
=a        
        for(@s){
			print OUT "<td>$_</td>";
		}
=cut
		print OUT "</tr>"
	}

	print OUT "</table>";
}