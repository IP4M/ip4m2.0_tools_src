#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use File::Copy;
use FindBin qw($Bin);
use Data::Dumper;
use File::Basename;
use Scalar::Util qw(looks_like_number);
use File::Copy::Recursive qw(dircopy);

my %opts;
GetOptions (\%opts, "i=s", "html=s", 
"lib1=s", 
"lib2=s", 
);

my $perl = $^X;


my $idmapping_bin = "$Bin/../pathways/idMapping_v2.pl";
my $pathway_bin = "$Bin/../pathways/keggEnrich_v2.pl";
my $enrichment_bin = "$Bin/../pathways/enrichment_v2.pl";


my $input = $opts{i};
my $lib1 = $opts{lib1};
my $lib2 = $opts{lib2};
my $html = $opts{html};

my $dir = dirname $opts{html};
my $num = 1;
my @output_files;

mkdir "01_ana";
my $cmd = "\"$perl\" \"$idmapping_bin\"  -i  \"$input\"  -html \"01_ana/outputs.html\" ";
run_cmd($cmd);
push @output_files,  (
    "01_ana/compounds_idmapping.txt", "idmap", "$num",
);
$num += 1;


mkdir "02_ana";
$cmd = "\"$perl\" \"$pathway_bin\" -html \"02_ana/outputs.html\" -i \"01_ana/compounds_idmapping.txt\" -method rbc -type hyperg  -lib \"$lib1\" ";
run_cmd($cmd);
push @output_files,  (
    "02_ana/Pathway_Result.txt", "path", "$num",
    "02_ana/Pathway_Bubbleplot.pdf", "pdf", "$num",
    "02_ana/Pathway_Barplot.pdf", "pdf", "$num",
);
$num += 1;


mkdir "03_ana";
$cmd = "\"$perl\" \"$enrichment_bin\" -html \"03_ana/outputs.html\" -i \"01_ana/compounds_idmapping.txt\"  -cmps 2 -lib \"$lib2\"";
run_cmd($cmd);
push @output_files,  (
    "03_ana/MSEA_Result.txt", "tabH", "$num",
    "03_ana/Enrichment_Barplot.pdf", "pdf", "$num",
    "03_ana/Network.pdf", "pdf", "$num",
);
copy("03_ana/MSEA_Network.json","$dir/$num.MSEA_Network.json");
$num += 1;




get_html_link_workflow($opts{html}, "Pathway and Enrichment Analysis Results", @output_files);










sub run_cmd {
        my $cmd = shift;
        print "$cmd\n";
        my $ret = system( $cmd );
        if ( $ret ){
                die "Error, died with $ret";
        }
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




sub get_html_link_workflow{
	my @param = @_;
	my $html = shift @param;
	my $title = shift @param;
	my @par = @param;

    my $out_dir = "$html.files";
    mkdir $out_dir;
    my $html_dir = dirname($html);

	open HTML, ">$html";
	print HTML "<html><head><title>$title</title></head><body><h3>Output Files:</h3><p><ul>\n";

    for(my $i = 0; $i <= $#par; $i+=3){
            my $file = $par[$i];
            my $type = $par[$i+1];
            my $num = $par[$i+2];
            my $name = basename $file;
            $name = "$num.$name";
            if( $type eq "tabH" ){
                copy($file, "$html_dir/$name");
                get_html_table($file, "$out_dir/$name.html", "y");
            }elsif(  $type eq "tabN" ){
                copy($file, "$html_dir/$name");
                get_html_table($file, "$out_dir/$name.html", "n");
            }elsif(  $type eq "txt" ){
                copy($file, "$html_dir/$name");
                txt2html($file, "$out_dir/$name.html");
            }elsif(  $type eq "pdf"){
                copy($file, "$html_dir/$name");
            }elsif(  $type eq "dir" ){
                dircopy(  $file , "$html_dir/$name") or die $!;
            }elsif(  $type eq "relation" ){
                copy($file, "$html_dir/$name");
                get_relation_table($file, "$out_dir/$name.html", "y");
            }elsif(  $type eq "idmap" ){
                copy($file, "$html_dir/$name");
                get_html_table_v2($file, "$out_dir/$name.html", "y");
            }elsif(  $type eq "path" ){
                copy($file, "$html_dir/$name");
                get_html_table_v3($file, "$out_dir/$name.html", "y");
            }
            else{
                copy($file, "$out_dir/$name");
            }
    }

    my $basename = basename $out_dir;
	for (my $i = 0; $i <= $#par; $i+=3){
            my $file = $par[$i];
            my $type = $par[$i+1];
            my $num = $par[$i+2];
            my $name = basename $file;
            $name = "$num.$name";
            if( $type eq "tabH" ){
                print HTML "<li><a href=\"$basename/$name.html\" target=\"_parent\">$name</a></li>\n";
            }elsif(  $type eq "tabN" ){
                print HTML "<li><a href=\"$basename/$name.html\" target=\"_parent\">$name</a></li>\n";
            }elsif(  $type eq "txt" ){
                print HTML "<li><a href=\"$basename/$name.html\" target=\"_parent\">$name</a></li>\n";
            }elsif( $type eq "pdf" ){
                print HTML "<li><a href=\"$name\" target=\"_parent\">$name</a></li>\n";
            }elsif(  $type eq "dir" ){
                print HTML "<li><a href=\"$name\" target=\"_parent\">$name</a></li>\n";
            }elsif(  $type eq "relation" ){
                print HTML "<li><a href=\"$basename/$name.html\" target=\"_parent\">$name</a></li>\n";
            }elsif(  $type eq "idmap" ){
                print HTML "<li><a href=\"$basename/$name.html\" target=\"_parent\">$name</a></li>\n";
            }elsif(  $type eq "path" ){
                print HTML "<li><a href=\"$basename/$name.html\" target=\"_parent\">$name</a></li>\n";
            }
            else{
                print HTML "<li><a href=\"$basename/$name\" target=\"_parent\">$name</a></li>\n";
            }
	        
	}
	print HTML "</ul></p>\n";


}






##check matrix and sample group data, group number must be 2
sub check_data_v3{
	my $matrix = shift;
	my $sample_group = shift;
    my $paired = shift;
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
        my %cnt2;
        open IN,"$sample_group";
        while (<IN>){
            chomp;
        	my @s = split /\t/, $_;
			die "Error: sample group file each line must be 2 columns!\n" if @s != 2;
	        $cnt{$s[0]} += 1;
            $cnt2{$s[1]} += 1;
	        die "Error: sample name: $s[0] in sample group file is not exist in the matrix file!\n" unless defined($hs{$s[0]});
	    }
	    close IN;

		foreach ( keys %cnt){
		     if ( $cnt{$_} > 1){
		        die "Error: sample name: $_ is not uniq in the sample group file!\n";
		     }
		}

        my @gps = keys %cnt2;
        if ( @gps != 2) {
            die "Error: groups number must be 2 in the sample group file!\n";
        }

        if( $paired eq "T"){
            if ( $cnt2{$gps[0]} != $cnt2{$gps[1]} ) {
                die "Error: for paired test, samples number must be equal in the two groups!\n";
            }
        }
	}
}


use Scalar::Util qw(looks_like_number);
sub isnumeric {  
        my $val = shift;  
        looks_like_number $val;
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





sub get_html_table_v3{
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
            if($i == $#s){
                print OUT "<td><a href=\"$s[$i]\" target=\"_blank\">link<a></td>";}
            else{
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