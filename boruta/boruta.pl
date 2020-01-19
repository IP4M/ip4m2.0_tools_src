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
GetOptions (\%opts, "input1=s", "input2=s",  "pValue=s", "maxRuns=s", 
 "html=s", 
 "scale=s",
);



my $input = $opts{input1};
my $group_file = $opts{input2};
my $zscal = $opts{scale};
my $pValue = $opts{pValue};
my $maxRuns = $opts{maxRuns};

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


check_data_v4($input, $group_file);
copy($input, '__matrix__');
copy($group_file, '__group__');
open RCMD, ">cmd.r";
print RCMD <<"RSCRIPT";

library(Boruta)
library(caret)
library(tibble)
library(magrittr)
library(dplyr)
library(stringr)
library(readr)
library(ggplot2)
library(ggrepel)
library(tidyr)
library(purrr)

zscal = $zscal

data <- read.table("__matrix__", header = T, sep = "\t", check.names = FALSE,quote="", row.names = 1)
group <- read.table("__group__", header = F, sep = "\t",check.names = FALSE,quote="",  row.names = 1, colClasses=c("character","character"))
colnames(group) = c("Group")

group\$Group = as.character(group\$Group)
group\$Group = as.factor(group\$Group)

sampleNames = rownames(group)
data = t(data[sampleNames])
data = as.data.frame(data)
if( zscal){
    data = scale(data)
}
group =as.data.frame(group)
data1 = cbind(data, group)

## rename variable
name_pairs = data.frame( oldname = colnames(data1), newname = make.names(colnames(data1)))
colnames(data1) = make.names(colnames(data1))
write.table(name_pairs, file = "pairs.txt", quote = FALSE, sep = "\t")

borutaRs <- Boruta(Group~.,data=data1, maxRuns=$maxRuns, pValue=$pValue )
decData <- borutaRs\$finalDecision \%>\%
    as.data.frame() \%>\%
    rownames_to_column("Metabolite") \%>\%
    set_colnames(c("Metabolite", "decision"))

statData <- attStats(borutaRs) \%>\%
    rownames_to_column("Metabolite")
outData <- statData \%>\%
    arrange(desc(medianImp)) \%>\%
    column_to_rownames("Metabolite")

write.csv(outData, "Decision_Info.csv")
write.table(outData, "Decision_Info.txt", quote = FALSE, sep = "\t" )


sortData <- borutaRs\$ImpHistory \%>\%
    t() \%>\%
    as.data.frame() \%>\%
    rownames_to_column("Metabolite") \%>\%
    rowwise() \%>\%
    do({
        result <- as.data.frame(.)
        data <- result[- 1] \%>\%
            unlist() \%>\%
            discard(~ is.infinite(.x))
        median <- median(data, na.rm = T)
        result\$median <- median
        result
    }) \%>\%
    arrange(median) \%>\%
    as.data.frame()

sortNames <- sortData\$Metabolite

decLelvels <- c("Shadow", "Rejected", "Tentative", "Confirmed")
fillColors <- setNames(c("#0000FF", "#FF0000", "#FFFF00", "#00FF00"), decLelvels)

imp <- borutaRs\$ImpHistory \%>\%
    t() \%>\%
    as.data.frame() \%>\%
    rownames_to_column("Metabolite") \%>\%
    gather("Sample", "Value", - Metabolite) \%>\%
    left_join(decData, by = c("Metabolite")) \%>\%
    mutate_at(vars("decision"), function(x){
        ifelse(is.na(x), "Shadow", as.character(x))
    }) \%>\%
    mutate(Metabolite = factor(Metabolite, levels = sortNames)) \%>\%
    mutate(decision = factor(decision, levels = decLelvels)) \%>\%
    as.data.frame()

p <- ggplot(imp, mapping = aes(x = Metabolite, y = Value, fill = decision)) +
    ylab("Importance") +
    xlab("") +
    theme_bw(base_size = 8.8, base_family = "Times") +
    theme(axis.text.x = element_text(angle = 45, size = 9, hjust = 1, vjust = 1), legend.position = 'top',
    legend.text = element_text(size = 9), legend.title = element_text(size = 11),
    axis.text.y = element_text(size = 8.8), axis.title.y = element_text(size = 11),
    axis.title.x = element_text(size = 12), panel.border = element_rect(size = 0.75)
    ) +
    geom_boxplot() +
    scale_fill_manual("finalDecision:", values = fillColors)

rowNum <- length(unique(imp\$Metabolite))
width <- max(2 + (rowNum - 10) * 0.25, 2) + 2
ggsave("Decision_Boxplot.pdf", p, width = width, height = 6, limitsize = FALSE)


RSCRIPT
close RCMD;
my $ret = system("R --restore --no-save < cmd.r");
if ( $ret ) {
        die "Error, died with $ret";
}


my $html_dir = dirname($html);
old_new_name_convert("Decision_Info.txt", "$html_dir/Decision_Info.txt");
&get_html_link_v3($html, "Boruta Results", 
    "$html_dir/Decision_Info.txt", "tabH", 
    "Decision_Boxplot.pdf", "pdf",  
);



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