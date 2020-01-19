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
GetOptions (\%opts, "i1=s", "i2=s",  
"html=s",  
);



my $input = $opts{i1};
my $group_file = $opts{i2};
my $html = $opts{html};


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
deal_NA_data( $input, "__noNA_data__" );
$input = "__noNA_data__";


check_data_v3($input, $group_file, "F");
copy($input, '__matrix__');
copy($group_file, '__group__');

open RCMD, ">cmd.r";
print RCMD <<"RSCRIPT";

library(tibble)
library(magrittr)
library(dplyr)
library(stringr)
library(readr)
library(ggplot2)
library(ggrepel)
library(tidyr)
library(purrr)
library(pROC)
library(ggpubr)


data <- read.table("__matrix__", header = T, sep = "\t", check.names = FALSE,quote="", row.names = 1, com='')
group <- read.table("__group__", header = F, sep = "\t",check.names = FALSE,quote="",  row.names = 1, com='', colClasses=c("character","character"))
colnames(group) = c("Group")

group\$Group = as.character(group\$Group)
group\$Group = as.factor(group\$Group)

sampleNames = rownames(group)
data = t(data[sampleNames])
data = as.data.frame(data)
group =as.data.frame(group)
data1 = cbind(data, group)

in_data = subset(data1, select=-Group)

getRoc <- function(name){
    x <- in_data[[name]] 
    y <- data1\$Group
    rocRs <- roc(y, x, ci = T, quiet=T)
    return(list(name = name, roc = rocRs))
}
rocs <- colnames(in_data) \%>\%
    map(getRoc)

aucData <- rocs \%>\%
    map_dfr(function(rocList){
        rocRs <- rocList\$roc
        name <- rocList\$name
        auc <- round(rocRs\$auc, 3)
        ci <- rocRs\$ci
        ci1 <- ci[1]
        ci2 <- ci[3]
        pointData <- coords(rocRs, "best", transpose = FALSE) \%>\%
            round(3) \%>\%
            as.data.frame()
        return(data.frame(Metabolite = name, AUC = auc, CI1 = ci1, CI2 = ci2, Thres = pointData\$threshold,
        Specificity = pointData\$specificity, Sensitivity = pointData\$sensitivity))
    }) \%>\%
    arrange(desc(AUC))

write.csv(aucData, "Single_Met_ROC.csv", row.names = F)
write.table(aucData, file = "Single_Met_ROC.txt", quote = FALSE, sep = "\t", row.names = F)

aucData <- rocs \%>\%
    map_dfr(function(rocList){
        rocRs <- rocList\$roc
        name <- rocList\$name
        auc <- round(rocRs\$auc, 3)
        return(data.frame(name = name, auc = auc))
    }) \%>\%
    arrange(desc(auc))

aucNames <- aucData\$name

breaks <- seq(0, 1, by = 0.2)

colors <- hsv(h = seq(0, 1, length = 100) * 0.8, s = 1, v = 1)
yBreaks <- seq(0, 1, 0.2)

getP <- function(name){
    list <- rocs \%>\%
    detect(function(rocList){
        rocList\$name == name
    })
    rocRs <- list\$roc
    ci <- rocRs\$ci
    ci1 <- ci[1] \%>\%
    round(3)
    ci2 <- ci[3] \%>\%
    round(3)
    auc <- round(rocRs\$auc, 3)
    pointData <- coords(rocRs, "best", transpose = FALSE) \%>\%
        round(3) \%>\%
        as.data.frame()

    plotData <- tibble(sensitivity = rocRs\$sensitivities, specificity = rocRs\$specificities, threshold = rocRs\$thresholds) \%>\%
        mutate(specificity = 1 - specificity) \%>\%
        arrange(sensitivity)

    p <- ggplot(data = plotData, aes(x = specificity, y = sensitivity, color = threshold)) +
        theme_bw(base_size = 8.8) +
        theme(plot.title = element_text(size = 14, hjust = 0.5, face = "bold"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        plot.margin = unit(c(1, 0.5, 1, 0.5), "cm"), panel.border = element_rect(size = 0.75)
        ) +
        ggtitle(str_replace(name, "/", "/\n")) +
        labs(colour = "Cutoff") +
        geom_line() +
        geom_abline(intercept = 0, slope = 1, color = "grey", linetype = 1, size = 0.4) +
        annotate("text", x = 0.5 + 0.01 , y = 0.5, label = str_c("AUC:", auc, "(", ci1, ", ", ci2, ")"), color = "black",
        hjust = 0, size = 3.5, family = "Times") +
        geom_point(data = pointData, aes(x = 1 - specificity, y = sensitivity), size = 1, color = "red") +
        geom_text(data = pointData, aes(x = 1 - specificity, y = sensitivity, label = paste0(threshold,
        "(", specificity, ", ", sensitivity, ")")), color = "black", hjust = 0, vjust = 1, size = 3, nudge_x = 0.015,
        nudge_y = - 0.015, family = "Times") +
        scale_x_continuous("1 - Specificity", breaks = yBreaks) +
        scale_y_continuous("Sensitivity", breaks = seq(0, 1, 0.2)) +
        scale_colour_gradientn(colours = colors)
    return(p)
}

p <- aucNames \%>\%
    map(getP)
                          
pdfFileName <- "Single_Met_ROC.pdf"
pdf(pdfFileName, 6, 6)
for (i in seq(1, length(p), 1)) {
    print(p[i])
}
dev.off()

                                         
RSCRIPT
close RCMD;
my $ret = system("R --restore --no-save < cmd.r");
if ( $ret ) {
        die "Error, died with $ret";
}



my $html_dir = dirname($html);
copy("Single_Met_ROC.txt", "$html_dir/Single_Met_ROC.txt");
&get_html_link_v3($html, "ROC Analysis Results", 
    "$html_dir/Single_Met_ROC.txt", "tabH", 
    "Single_Met_ROC.pdf", "pdf", 
);


sub deal_roc_summary{
    my $output =shift;
    open OUT, ">$output";
    open IN,"names.txt";
    my @variables;
    while(<IN>){
        my $variable = trim($_);
        push @variables, $variable;
    }

    close IN;
    my @data;
    open IN,"roc_sum.txt";
    while(my $blank = <IN>){
        my @values;
        for(1..5){
            my $line = <IN>;
            my $value = (split /:/, $line)[1];
            $value = trim($value);
            push @values, $value;
         
        }
        push @data, [@values];
    }
    close IN;

    print OUT "\tNumber of cut-points\tBest cut-point\tSensitivity\tSpecificity\tAUC\n";
    for(my $i = 0; $i <= $#variables; $i++){
        print OUT "$variables[$i]\t";
        my @values = @{$data[$i]};
        my $tmp = join "\t", @values;
        print OUT "$tmp\n";
    }

    close OUT;
}


sub trim
{
        my $string = shift;
        $string =~ s/^\s+//;
        $string =~ s/\s+$//;
        return $string;
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



##check matrix and sample group data

sub check_data_v5{
	my $matrix = shift;
	my $sample_group = shift;
    my $labpos = shift;
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
        unless ( defined($cnt2{$labpos})) {
            die "Error: the group name: $labpos that defines a 'positive' event is not in the sample group file!\n";
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
            }elsif(  $type eq "txt_v2" ){
                copy($file, "$out_dir/$name");
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
            }elsif(  $type eq "txt_v2" ){
                print HTML "<li><a href=\"$basename/$name\" target=\"_parent\">$name</a></li>\n";
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