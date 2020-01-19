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
GetOptions (\%opts, "i=s", "g=s",  "html=s",  "zscal=s", 
);

my $matrix = $opts{i};
my $zscal = $opts{zscal};
my $html = $opts{html};
my $group_file = $opts{g};


use Scalar::Util qw(looks_like_number);
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
deal_NA_data( $matrix, "__noNA_data__" );
$matrix = "__noNA_data__";


check_data_v4($matrix, $group_file, "F");

copy($matrix, '__matrix__');
copy($group_file, '__group__');
get_csv_group_file($group_file, '__group2__');

open RCMD, ">cmd.r";
print RCMD <<"RSCRIPT";

library(ropls)
library(tibble)
library(magrittr)
library(dplyr)
library(stringr)
library(readr)
library(ggplot2)
library(ggrepel)
library(tidyr)
zscal = $zscal
xMN <- read.table(quote="", "__matrix__", check.names = FALSE, header = TRUE, row.names = 1, sep = "\t", comment.char = "")
tmpDf = xMN
samDF <- read.table(quote="", "__group__", check.names = FALSE, header = F, row.names = 1,  com = '', sep = "\t")
colnames(samDF) = c("Group")

samDF\$Group = as.character(samDF\$Group)
samDF\$Group = as.factor(samDF\$Group)

sampleNames = rownames(samDF)
sampleSize = length(sampleNames)
xMN = xMN[rownames(samDF)]
xMN = t(xMN)
yMCN <- matrix(samDF[, "Group"], ncol = 1, dimnames = list(rownames(xMN), "Group"))

if( zscal){
    xMN = scale(xMN)
}
crossValI <- min(nrow(xMN), 7)

comp <- length(unique(yMCN))
if (comp < 3) {
    comp <- 3
}
plsdaRs <- opls(x = xMN, y = yMCN, predI = comp, plotL = F, crossvalI = crossValI)

modelDf <- plsdaRs\@modelDF \%>\%
    rownames_to_column("pcName") \%>\%
    slice(1 : 5) \%>\%
    mutate(pcName = str_replace(pcName, "p", "C")) \%>\%
    column_to_rownames("pcName")


plotData <- plsdaRs\@scoreMN \%>\%
    as.data.frame() %>%
    rownames_to_column("SampleID") \%>\%
    select("SampleID", num_range("p", 1 : 5)) \%>\%
    rename_at(vars(- c("SampleID")), function(x){
        str_replace(x, "p", "C")
    }) \%>\%
    column_to_rownames("SampleID")


write.csv(plotData, "PLSDA_Score.csv")
write.table(plotData, "PLSDA_Score.txt", quote = FALSE, sep = "\t")
write.csv(modelDf, "PLSDA_R2X_R2Y_Q2.csv", quote = FALSE)
write.table(modelDf, "PLSDA_R2X_R2Y_Q2.txt", quote = FALSE, sep = "\t")


sampleInfo <- read.csv("__group2__", header = T, stringsAsFactors = F) \%>\%
    select(c("SampleID", "ClassNote")) \%>\%
    mutate(ClassNote = as.character(ClassNote))
plotData <- read_csv("PLSDA_Score.csv") \%>\%
    rename(SampleID = X1) \%>\%
    inner_join(sampleInfo, by = c("SampleID")) \%>\%
    mutate(ClassNote = factor(ClassNote, levels = unique(ClassNote)))

modelDf <- read.csv("PLSDA_R2X_R2Y_Q2.csv", header = T, row.names = 1)

p <- ggplot(plotData, mapping = aes(x = C1, y = C2, color = ClassNote, label = SampleID, fill = ClassNote)) +
    xlab(paste0("Component 1 (", modelDf["C1", "R2X"] * 100 , "%)")) +
    ylab(paste0("Component 2 (", modelDf["C2", "R2X"] * 100 , "%)")) +
    ggtitle("") +
    theme_bw(base_size = 8.8, base_family = "Times") +
    theme(axis.text.x = element_text(size = 9, vjust = 0.5),
    axis.text.y = element_text(size = 8.8), legend.position = 'right',
    axis.title.y = element_text(size = 12), legend.margin = margin(t = 0.3, b = 0.1, unit = 'cm'),
    legend.text = element_text(size = 6), axis.title.x = element_text(size = 12),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 13)
    ) +
#0 line
    geom_vline(aes(xintercept = 0), colour = "#BEBEBE", linetype = "solid") +
    geom_hline(aes(yintercept = 0), colour = "#BEBEBE", linetype = "solid") +
#point
    geom_point(aes(colour = ClassNote), size = 4, stroke = 0) +
    stat_ellipse(aes(fill = ClassNote), colour = NA, size = 0.3, level = 0.95, type = "norm",
    geom = "polygon", alpha = 0.2, show.legend = F) +
    geom_text_repel(segment.size=0.2,size = 2, family = "Times")

ggsave("PLSDA_Score_2D_Label.pdf", p, width = 5, height = 4)


p <- ggplot(plotData, mapping = aes(x = C1, y = C2, color = ClassNote)) +
    xlab(paste0("Component 1 (", modelDf["C1", "R2X"] * 100 , "%)")) +
    ylab(paste0("Component 2 (", modelDf["C2", "R2X"] * 100 , "%)")) +
    ggtitle("") +
    theme_bw(base_size = 8.8, base_family = "Times") +
    theme(axis.text.x = element_text(size = 9, vjust = 0.5),
    axis.text.y = element_text(size = 8.8), legend.position = 'right',
    axis.title.y = element_text(size = 12), legend.margin = margin(t = 0.3, b = 0.1, unit = 'cm'),
    legend.text = element_text(size = 6), axis.title.x = element_text(size = 12),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 13)
    ) +
#0 line
    geom_vline(aes(xintercept = 0), colour = "#BEBEBE", linetype = "solid") +
    geom_hline(aes(yintercept = 0), colour = "#BEBEBE", linetype = "solid") +
#point
    geom_point(aes(colour = ClassNote), size = 4, stroke = 0) +
    stat_ellipse(aes(fill = ClassNote), colour = NA, size = 0.3, level = 0.95, type = "norm",
    geom = "polygon", alpha = 0.2, show.legend = F)


ggsave("PLSDA_Score_2D.pdf", p, width = 5, height = 4)


RSCRIPT
close RCMD;
my $ret = system("R --restore --no-save < cmd.r");
if ( $ret ) {
        die "Error, died with $ret";
}


my $html_dir = dirname($html);
getRTable_oplsda("PLSDA_Score.txt", "$html_dir/PLSDA_Score.txt");
getRTable_oplsda("PLSDA_R2X_R2Y_Q2.txt", "$html_dir/PLSDA_R2X_R2Y_Q2.txt");

&get_html_link_v3($html, "PLS-DA Analysis Results", 
    "$html_dir/PLSDA_Score.txt", "tabH", 
    "$html_dir/PLSDA_R2X_R2Y_Q2.txt", "tabH", 

    "PLSDA_Score_2D_Label.pdf", "pdf", 
    "PLSDA_Score_2D.pdf", "pdf", 
);





sub getRTable_oplsda{
    my $in = shift;
    my $out = shift;
    open IN,"$in";
    open OUT,">$out";
    my $h = <IN>;
    print OUT "\t$h";
    while (<IN>){
          print OUT "$_";
    }
    close IN;
    close OUT;
}



# getRTable(in, out)
sub getRTable_v10{
    my $in = shift;
    my $out = shift;

    open IN,"$in";
    open OUT,">$out";
    my $h = <IN>;
    chomp $h;
    my @s = split /\t/, $h;
    pop @s;
    $h = join "\t", @s;
    print OUT "\t$h\n";
    while (<IN>){
          chomp;
          my @s = split /\t/, $_;
          pop @s;
          my $line = join "\t", @s;
          print OUT "$line\n";
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



sub get_csv_group_file{
    my $in = shift;
    my $out = shift;
    open IN,"$in";
    open OUT,">$out";
    print OUT "SampleID,ClassNote\n";
    while (<IN>){
          chomp;
          my @s = split /\t/, $_;
          my $l = join ",", @s;
          print OUT "$l\n";
    }
    close IN;
    close OUT;
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
                get_html_table($file, "$out_dir/$name.html", "y");
            }elsif(  $type eq "tabN" ){
                get_html_table($file, "$out_dir/$name.html", "n");
            }elsif(  $type eq "txt" ){
                txt2html($file, "$out_dir/$name.html");
            }elsif(  $type eq "pdf"){
                copy($file, "$html_dir/$name");
            }elsif(  $type eq "dir" ){
                dircopy(  $file , "$html_dir/$name") or die $!;
            }
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
            }
            else{
                print HTML "<li><a href=\"$basename/$name\" target=\"_parent\">$name</a></li>\n";
            }
	        
	}
	print HTML "</ul></p>\n";


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
