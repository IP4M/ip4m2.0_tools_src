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
"lib=s",    "type=s", "method=s",  
);



my $input = $opts{i};
my $html = $opts{html};
my $lib = $opts{lib};
my $type = $opts{type};
my $method = $opts{method};
my $pathway_libs = "$Bin/mylibs";


process_idmapping( $input, "compounds_idmapping.txt" );
cmp_color( $input, "color.txt" );

my $lib_path = "$Bin/metabo_base_with_id_set.r";
copy( $lib_path, "./__libs.r__" );

open RCMD, ">cmd.r";
print RCMD <<"RSCRIPT";
library(magrittr)
library(tidyr)
library(tibble)
library(tidyverse)
library(viridis)

createWhenNoExist <- function(f){
    ! dir.exists(f) && dir.create(f)
}

source("./__libs.r__")
data = read.table("compounds_idmapping.txt", header=T, com='', quote="",sep="\t", check.names=F)

nm.mat = data.frame(query=data[,1], hmdb=data[,2], kegg=data[,5], hmdbid=data[,3])

mSet <- InitDataObjects("conc", "pathora", FALSE)
mSet <- Setup.MapData(mSet, nm.mat)
mSet <- CrossReferencing(mSet, "kegg")

fcData <- read_tsv("color.txt")

mSet <- SetMetabolomeFilter(mSet, F)

dataDir = "$pathway_libs"
mSet <- SetKEGG.PathLib(mSet, "$lib", dataDir = dataDir)
mSet <- SetMetabolomeFilter(mSet, F)
nodeImp = "$method"
method = "$type"
mSet <- CalculateOraScore(mSet, nodeImp, method, dataDir = dataDir, fcData = fcData, FALSE)

mSet <- PlotPathSummary(mSet, imgName = "Pathway_Bubbleplot.pdf", width = 7, height = 5)

data <- read_csv("Pathway_Result.csv") \%>\%
rename(Metabolite = X1)

plotData <- data \%>\%
    rename(p = `Raw P`) \%>\%
    mutate(p = as.numeric(p), Impact = as.numeric(Impact)) \%>\%
    mutate(logP = - log(p)) \%>\%
    arrange(logP) \%>\%
    slice(1 : 30) \%>\%
    mutate(Metabolite = factor(Metabolite, levels = unique(Metabolite)))

p <- ggplot(plotData, mapping = aes(x = Metabolite, y = Impact)) +
    ylab("Impact") +
    xlab("") +
    theme_bw(base_size = 8.8, base_family = "Times") +
    theme(axis.text.x = element_text(size = 10, hjust = 1, vjust = 1), legend.position = 'right',
    legend.text = element_text(size = 9), legend.title = element_text(size = 11), axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 11), axis.title.x = element_text(size = 12), panel.grid.major.x = element_blank(),
    panel.border = element_rect(size = 0.75), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(),
    ) +
    geom_col(aes(fill = logP)) +
    coord_flip() +
    scale_fill_gradient("-log(P)", low = "yellow", high = "red")

ggsave(limitsize = FALSE,"Pathway_Barplot.pdf", p, width = 6, height = 8)



RSCRIPT
close RCMD;
my $ret = system("R --restore --no-save < cmd.r");
if ( $ret ) {
        die "Error, died with $ret";
}




my $html_dir = dirname($html);

getRTable("Pathway_Result.txt", "$html_dir/Pathway_Result.txt" );
&get_html_link_v3($html, "Pathway Analysis Results",  
    "$html_dir/Pathway_Result.txt", "tabH", 
    "Pathway_Bubbleplot.pdf", "pdf", 
    "Pathway_Barplot.pdf", "pdf", 

);




sub cmp_color{
	my $in = shift;
	my $out = shift;
	open IN,"$in";
	open OUT,">$out";
    my $header = <IN>;
    print OUT "KEGG\tcol\n";

    while(<IN>){
        chomp;
        my @s = split /\t/, $_;
        if ($s[1] =~/No Hit/) {
            
        }else{
            if ( defined($s[4]) && $s[4]=~/^C/ ) {
                print OUT "$s[4]\tDeepSkyBlue\n";
            }
        }
    }

	close IN;
	close OUT;    
}




sub process_idmapping{
	my $in = shift;
	my $out = shift;
	open IN,"$in";
	open OUT,">$out";
    my $header = <IN>;
    chomp $header;
    my @s = split /\t/, $header;
    $header = join "\t", @s[0..4];
    print OUT "$header\n";

    while(<IN>){
        chomp;
        my @s = split /\t/, $_;
        if ($s[1] =~/No Hit/) {
            print OUT "$s[0]\tNA\tNA\tNA\tNA\n";
        }else{
            for(my $i = 0; $i <= 4; $i++){
                if( !defined($s[$i]) || $s[$i]=~/^\s*$/ ){
                    $s[$i] = "NA";
                }
            }
            my $line = join "\t", @s[0..4];
            print OUT "$line\n";
        }
    }

	close IN;
	close OUT;
}




sub transpose{
my $input = shift;
my $output = shift;

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
}


sub transpose2{
my $input = shift;
my $output = shift;

open IN,"$input";

my @matrix;
my $header = <IN>;
chomp $header;
my @s = split /\t/, $header;
$s[0] = "Support Vectors";
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
}


sub old_new_name_convert{

    my $file = shift;
    my $out = shift;
    my %hs;

    open IN,"pairs.txt";
    <IN>;
    while(<IN>){
        chomp;
        my @s = split /\t/, $_;
        $hs{$s[2]} = $s[1];
    }
    close IN;
    open IN,"$file";
    open OUT,">$out";
    my $h = <IN>;
    print OUT "\t$h";
    while(<IN>){
        my @s = split /\t/, $_;
        $s[0] = $hs{$s[0]};
        my $tmp  = join "\t", @s;
        print OUT "$tmp";
    }
    close IN;
    close OUT;
}


# getRTable(in, out)
sub getRTable{
    my $in = shift;
    my $out = shift;
    my $name = shift;
    open IN,"$in";
    open OUT,">$out";
    my $h = <IN>;
    if( defined($name)){
        print OUT "$name\t$h";
    }else{
        print OUT "\t$h";
    }
    while (<IN>){
          print OUT "$_";
    }
    close IN;
    close OUT;
}


sub getRTable_v2{
    my $in = shift;
    my $out = shift;
    open IN,"$in";
    open OUT,">$out";
    my $h = <IN>;
    $h=~s/holm/p_holm/;
    $h=~s/bonferroni/p_bonferroni/;
    print OUT "\t$h";
    while (<IN>){
          print OUT "$_";
    }
    close IN;
    close OUT;
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





##check matrix and sample group data

sub check_data_v4{
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
        if ( @gps < 2) {
            die "Error: groups number must be >= 2 in the sample group file!\n";
        }
	}
}
use Scalar::Util qw(looks_like_number);
sub isnumeric {  
        my $val = shift;  
        looks_like_number $val;
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




sub get_html_link_v3{
	my @param = @_;
	my $html = shift @param;
	my $title = shift @param;
	my @par = @param;

    my $out_dir = "$html.files";
    mkdir $out_dir;
    my $html_dir = dirname($html);

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
            }elsif(  $type eq "pdf"){
                copy($file, "$html_dir/$name");
            }elsif(  $type eq "dir" ){
                dircopy(  $file , "$html_dir/$name") or die $!;
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
            }elsif( $type eq "pdf" ){
                print HTML "<li><a href=\"$name\" target=\"_parent\">$name</a></li>\n";
            }elsif(  $type eq "dir" ){
                print HTML "<li><a href=\"$name\" target=\"_parent\">$name</a></li>\n";
            }elsif(  $type eq "relation" ){
                print HTML "<li><a href=\"$basename/$name.html\" target=\"_parent\">$name</a></li>\n";
            }
            else{
                print HTML "<li><a href=\"$basename/$name\" target=\"_parent\">$name</a></li>\n";
            }
	        
	}
	print HTML "</ul></p>\n";


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