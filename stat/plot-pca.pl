#!/usr/bin/perl

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

GetOptions( \%opts,"i=s",  "g=s",  "zscal=s", "html=s");
my $usage = <<"USAGE";
       Program : $0   
USAGE

die $usage if(!($opts{i}));

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




check_data_v4($matrix, $group_file);

copy($matrix, '__matrix__');
if ( defined($group_file) ) {
    get_csv_group_file($group_file, '__group__');
}

my $group;
if(defined($group_file)){
    $group = "TRUE";
}else{
    $group = "FALSE";
}


open RCMD, ">cmd.r";
print RCMD <<__EORSCRIPT__;
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
group = $group

data <-read.table("__matrix__", header=T, quote="",com='', sep="\t", row.names=1, check.names=F )
tmpDF = data
data = t(data)

if( zscal){
    data = scale(data)
}

crossValI <- min(nrow(data), 7)
pcaRs <- opls(data, plotL = F, crossvalI = crossValI)
df <- pcaRs\@modelDF \%>\% as_tibble()
if (nrow(df) < 5) {
    pcaRs <- opls(data, plotL = F, crossvalI = crossValI, predI = 5)
}

parameterData <- pcaRs\@modelDF \%>\%
    select(c(1 : 2)) \%>\%
    set_colnames(c("R2", "Cumulative R2")) \%>\%
    mutate(pcName = paste0("PC", 1 : n())) \%>\%
    filter(pcName \%in\% str_c("PC", 1 : 5)) \%>\%
    column_to_rownames("pcName") \%>\%
    t()
write.table(parameterData, "PCA_R2.txt", quote = FALSE, sep = "\t")
write.csv(parameterData,  "PCA_R2.csv", quote = FALSE)

pcData <- pcaRs\@scoreMN \%>\%
    as.data.frame() \%>\%
    rownames_to_column("SampleID") \%>\%
    select("SampleID", num_range("p", 1 : 5)) \%>\%
    rename_at(vars(- c("SampleID")), function(x){
        str_replace(x, "p", "PC")
    }) \%>\%
    column_to_rownames("SampleID")
write.table(pcData, "PCA_Score.txt", quote = FALSE, sep = "\t")
write.csv(pcData, "PCA_Score.csv", quote = FALSE)


pcDataFileName <- "PCA_Score.csv"
pcData <- read_csv(pcDataFileName) \%>\%
    rename(SampleID = X1)
parameterFileName <- "PCA_R2.csv"
parameterData <- read.csv(parameterFileName, header = T, stringsAsFactors = F, comment.char = "")

pNames <- pcData \%>\%
    column_to_rownames("SampleID") \%>\%
    colnames()
cn <- combn(pNames, 2)
if (group ) {
    sampleInfo <- read_csv("__group__") \%>\%
        select(c("SampleID", "ClassNote")) \%>\%
        mutate(ClassNote = factor(ClassNote, levels = unique(ClassNote)))
    for (i in 1 : ncol(cn)) {
        row <- cn[, i]
        p1Name <- row[1]
        p2Name <- row[2]
        p1Index <- str_replace(p1Name, "PC", "")
        p2Index <- str_replace(p2Name, "PC", "")
        pc12 <- pcData \%>\%
            select(c("SampleID", p1Name, p2Name)) \%>\%
            set_colnames(c("SampleID", "p1", "p2")) \%>\%
            left_join(sampleInfo, by = c("SampleID"))
        
        impoPc1 <- parameterData[1, str_replace(p1Name, "p", "PC")]
        impoPc1 <- round(impoPc1 * 100, 2)
        impoPc2 <- parameterData[1, str_replace(p2Name, "p", "PC")]
        impoPc2 <- round(impoPc2 * 100, 2)
        p <- ggplot(pc12, mapping = aes(x = p1, y = p2, label = SampleID, color = ClassNote)) +
            xlab(paste0("PC", p1Index, "(", impoPc1, "%)", sep = "")) +
            ylab(paste0("PC", p2Index, "(", impoPc2, "%)", sep = "")) +
            theme_bw(base_size = 8.8, base_family = "Times") +
            theme(axis.text.x = element_text(size = 9, vjust = 0.5),
            axis.text.y = element_text(size = 8.8), legend.position = 'right',
            axis.title.y = element_text(size = 11), legend.margin = margin(t = 0.3, b = 0.1, unit = 'cm'),
            legend.text = element_text(size = 6), axis.title.x = element_text(size = 11),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            ) +
             #0 line
            geom_vline(aes(xintercept = 0), colour = "#BEBEBE", linetype = "solid") +
            geom_hline(aes(yintercept = 0), colour = "#BEBEBE", linetype = "solid") +
            #point
            geom_point(aes(colour = factor(ClassNote)), size = 3, stroke = 0) +
            stat_ellipse(aes(fill = ClassNote), colour = NA, size = 0.3, level = 0.95, type = "norm",
            geom = "polygon", alpha = 0.2, show.legend = F) +
            geom_text_repel(segment.size=0.2,size = 2, family = "Times")
        fileName <- paste0("PC", p1Index, p2Index, "_Score_2D_Label.pdf")
        ggsave(fileName, p, width = 5, height = 4)
    }
    for (i in 1 : ncol(cn)) {
        row <- cn[, i]
        p1Name <- row[1]
        p2Name <- row[2]
        p1Index <- str_replace(p1Name, "PC", "")
        p2Index <- str_replace(p2Name, "PC", "")
        pc12 <- pcData \%>\%
            select(c("SampleID", p1Name, p2Name)) \%>\%
            set_colnames(c("SampleID", "p1", "p2")) \%>\%
            left_join(sampleInfo, by = c("SampleID"))
        
        impoPc1 <- parameterData[1, str_replace(p1Name, "p", "PC")]
        impoPc1 <- round(impoPc1 * 100, 2)
        impoPc2 <- parameterData[1, str_replace(p2Name, "p", "PC")]
        impoPc2 <- round(impoPc2 * 100, 2)
        p <- ggplot(pc12, mapping = aes(x = p1, y = p2, label = SampleID, color = ClassNote)) +
            xlab(paste0("PC", p1Index, "(", impoPc1, "%)", sep = "")) +
            ylab(paste0("PC", p2Index, "(", impoPc2, "%)", sep = "")) +
            theme_bw(base_size = 8.8, base_family = "Times") +
            theme(axis.text.x = element_text(size = 9, vjust = 0.5),
            axis.text.y = element_text(size = 8.8), legend.position = 'right',
            axis.title.y = element_text(size = 11), legend.margin = margin(t = 0.3, b = 0.1, unit = 'cm'),
            legend.text = element_text(size = 6), axis.title.x = element_text(size = 11),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            ) +
             #0 line
            geom_vline(aes(xintercept = 0), colour = "#BEBEBE", linetype = "solid") +
            geom_hline(aes(yintercept = 0), colour = "#BEBEBE", linetype = "solid") +
            #point
            geom_point(aes(colour = factor(ClassNote)), size = 3, stroke = 0) +
            stat_ellipse(aes(fill = ClassNote), colour = NA, size = 0.3, level = 0.95, type = "norm",
            geom = "polygon", alpha = 0.2, show.legend = F) 
        fileName <- paste0("PC", p1Index, p2Index, "_Score_2D.pdf")
        ggsave(fileName, p, width = 5, height = 4)
    }     
}else{
    for (i in 1 : ncol(cn)) {
        row <- cn[, i]
        p1Name <- row[1]
        p2Name <- row[2]
        p1Index <- str_replace(p1Name, "PC", "")
        p2Index <- str_replace(p2Name, "PC", "")
        pc12 <- pcData \%>\%
            select(c("SampleID", p1Name, p2Name)) \%>\%
            set_colnames(c("SampleID", "p1", "p2")) 
        
        impoPc1 <- parameterData[1, str_replace(p1Name, "p", "PC")]
        impoPc1 <- round(impoPc1 * 100, 2)
        impoPc2 <- parameterData[1, str_replace(p2Name, "p", "PC")]
        impoPc2 <- round(impoPc2 * 100, 2)
        p <- ggplot(pc12, mapping = aes(x = p1, y = p2, label = SampleID)) +
            xlab(paste0("PC", p1Index, "(", impoPc1, "%)", sep = "")) +
            ylab(paste0("PC", p2Index, "(", impoPc2, "%)", sep = "")) +
            theme_bw(base_size = 8.8, base_family = "Times") +
            theme(axis.text.x = element_text(size = 9, vjust = 0.5),
            axis.text.y = element_text(size = 8.8), legend.position = 'right',
            axis.title.y = element_text(size = 11), legend.margin = margin(t = 0.3, b = 0.1, unit = 'cm'),
            legend.text = element_text(size = 6), axis.title.x = element_text(size = 11),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),

            ) +
            #0 line
            geom_vline(aes(xintercept = 0), colour = "#BEBEBE", linetype = "solid") +
            geom_hline(aes(yintercept = 0), colour = "#BEBEBE", linetype = "solid") +
            #point
            geom_point(colour="blue", size = 3, stroke = 0) +
            geom_text_repel(segment.size=0.2,size = 2, family = "Times")
        fileName <- paste0("PC", p1Index, p2Index, "_Score_2D_Label.pdf")
        ggsave(fileName, p, width = 5, height = 4)
    }
}

lineData <- parameterData \%>\%
    column_to_rownames("X") \%>\%
    t() \%>\%
    as.data.frame() \%>\%
    select(c("R2")) \%>\%
    rename(value = R2) \%>\%
    mutate(value=value*100) \%>\%
    mutate(sum = cumsum(value)) \%>\%
    rownames_to_column("index")

plotData <- lineData \%>\%
    gather("Class", "value", - index)

labelData <- plotData \%>\%
    select(- c("Class")) \%>\%
    as.data.frame() \%>\%
    unique() \%>\%
    mutate(label = paste0(value, "\%"))

p <- ggplot(plotData, aes(x = index, y = value)) +
    xlab("PC index") +
    ylab("Variance Explained (%)") +
    ggtitle("Scree plot") +
    theme_bw(base_size = 8.8, base_family = "Times") +
    theme(axis.text.x = element_text(size = 9, vjust = 0.5),
    axis.text.y = element_text(size = 8.8), legend.position = 'none',
    axis.title.y = element_text(size = 11), legend.margin = margin(t = 0.3, b = 0.1, unit = 'cm'),
    legend.text = element_text(size = 6), axis.title.x = element_text(size = 11),
    panel.grid.major.y = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 13)
    ) +
    geom_line(aes(x = index, y = value, group = 1), lineData, size = 0.5, linetype = 1, color = "blue") +
    geom_line(aes(x = index, y = sum, group = 1), lineData, size = 0.5, linetype = 1, color = "green") +
    geom_point(size = 2, pch = 21, color = "red") +
    geom_text_repel(segment.size=0.2,data = labelData, aes(x = index, y = value, label = label), color = "black", size = 2, family = "Times")
ggsave("PCA_Screeplot.pdf", p, width = 5, height = 4)


__EORSCRIPT__
close RCMD;
&run_cmd("R --restore --no-save<cmd.r");


my $html_dir = dirname($html);
getRTable("PCA_Score.txt", "$html_dir/PCA_Score.txt");
getRTable("PCA_R2.txt", "$html_dir/PCA_R2.txt");
if($group eq "TRUE"){
    &get_html_link_v3($html, "PCA Analysis Results", "$html_dir/PCA_Score.txt", "tabH", "$html_dir/PCA_R2.txt", "tabH", 
    "PCA_Screeplot.pdf", "pdf", 
    "PC12_Score_2D_Label.pdf", "pdf", 
    "PC12_Score_2D.pdf", "pdf", 
    "PC13_Score_2D_Label.pdf", "pdf", 
    "PC13_Score_2D.pdf", "pdf", 
    "PC23_Score_2D_Label.pdf", "pdf", 
    "PC23_Score_2D.pdf", "pdf", 
    );
}else{
    &get_html_link_v3($html, "PCA Analysis Results", "$html_dir/PCA_Score.txt", "tabH", "$html_dir/PCA_R2.txt", "tabH", 
    "PCA_Screeplot.pdf", "pdf", 
    "PC12_Score_2D_Label.pdf", "pdf", 
    "PC13_Score_2D_Label.pdf", "pdf", 
    "PC23_Score_2D_Label.pdf", "pdf", 
    );
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
            }else{
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


sub getRTable{
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


sub isnumeric {  
        my $val = shift;  
        looks_like_number $val;
}  


sub run_cmd {
        my $cmd = shift;
#        print "$cmd\n";

        my $ret = system( $cmd );
        if ( $ret ){
                die "Error, died with $ret";
        }

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