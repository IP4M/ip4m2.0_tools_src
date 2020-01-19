#!/usr/bin/perl
use strict;
use warnings;
use File::Copy;
use Getopt::Long;
use FindBin qw($Bin);
my $r_path = "$Bin/../../bin/R/bin";
$ENV{'PATH'} = "$r_path;$ENV{'PATH'}";

my %opts;
GetOptions (\%opts,"w=s",  "h=s", "i=s", "s=s", "min=s", "pdf=s", "outline=s", "html=s");
my $usage = <<"USAGE";
       Program : $0
       Version : 1.0
       Discription: 
       Usage :perl $0 [options]
                   -i  input fpkm matrix
                   -s  sample group file
                   -w  width
                   -h  height
                   -min  counts >= min fpkm
                   -pdf  output pdf
                   -outline  default F
USAGE

die $usage if ( !$opts{i} );

my $matrix_file = $opts{i};
my $sample_group_file = $opts{s};
my $width = defined($opts{w}) ? $opts{w} : 10;
my $height = defined($opts{h}) ? $opts{h} : 7;

 $width =  $width / 2.54;
 $height =  $height / 2.54;

my $min = defined($opts{min}) ? $opts{min} : 1;
my $outline = defined($opts{outline}) ? $opts{outline} : "F";
my $pdf = $opts{pdf};

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
deal_NA_data( $matrix_file, "__noNA_data__" );
$matrix_file = "__noNA_data__";


&check_data_v4( $matrix_file, $sample_group_file );

copy($matrix_file, '__matrix__');
if ( defined($sample_group_file) ){
	copy($sample_group_file, 'sample_group');
}

open RCMD, ">cmd.r";
print RCMD "library(\"ggplot2\")";

unless ( defined($sample_group_file) ){
print RCMD <<"RSCRIPT";

pdf("out.pdf", width=$width,height=$height)
data = read.table("__matrix__", header=T, com='', sep="\t", row.names=1, check.names=F, quote="")
    head(data)
data = as.matrix(data)
    head(data)
boxplot_data = data
boxplot_data[boxplot_data<$min] = NA
boxplot_data = log2(boxplot_data+1)

#boxplot(boxplot_data[,1:length(boxplot_data[1,])], outline=$outline, main="fpkm boxplot ( analize fpkm >= $min )", ylab="log2(fpkm+1)" )
require(reshape2)
plot_data = melt(boxplot_data)
#plot_data
head(plot_data)
ggplot(plot_data, aes(x=Var2, y=value)) + geom_violin() + ylab("log2(v+1)") + xlab("") + ggtitle("Violinplot ( analize v >= $min )")+ theme(panel.background = element_rect(fill = NA), panel.grid.major = element_line(colour = NA ), panel.grid.minor = element_line(colour = NA), panel.border = element_rect( colour = "gray30", fill = NA, size=1.5 ) )


RSCRIPT
close RCMD;

} else {

print RCMD <<"RSCRIPT";

#pdf("out.pdf", width=$width,height=$height)
data = read.table("__matrix__", header=T, com='', sep="\t", row.names=1, check.names=F, quote="")
data = as.matrix(data)
samples_data = read.table("sample_group", header=F, check.names=F, quote="")
samples_data = data.frame( samples_data[,2], samples_data[,1] )
colnames(samples_data) = c("V1", "V2")
sample_types = as.character(unique(samples_data[,1]))
data = data[, colnames(data) %in% samples_data[,2], drop=F ]
nsamples = length(sample_types)
sample_colors = rainbow(nsamples)
names(sample_colors) = sample_types
sample_type_list = list()
for (i in 1:nsamples) {
    samples_want = samples_data[samples_data[,1]==sample_types[i], 2]
    sample_type_list[[sample_types[i]]] = as.vector(samples_want)
}
sample_factoring = colnames(data)
for (i in 1:nsamples) {
    sample_type = sample_types[i]
    replicates_want = sample_type_list[[sample_type]]
    sample_factoring[ colnames(data) %in% replicates_want ] = sample_type
}
# reorder according to sample type.
tmp_sample_reordering = order(sample_factoring)
data = data[,tmp_sample_reordering,drop=F]
sample_factoring = sample_factoring[tmp_sample_reordering]
# set up barplot colors:
sample_cols = rainbow(nsamples)
barplot_cols = c()
for (i in 1:nsamples) {
    barplot_cols[ sample_factoring %in% sample_types[i] ] = sample_cols[i]
}
names(barplot_cols) = colnames(data)
boxplot_data = data

boxplot_data[boxplot_data<$min] = NA
boxplot_data = log2(boxplot_data+1)

#boxplot( boxplot_data[,1:length(boxplot_data[1,])], outline=$outline, main="FPKM Boxplot ( analize fpkm >= $min )", ylab="log2(fpkm+1)", border=barplot_cols )
#legend("topright", pch=15,legend=names(sample_colors), col=sample_colors, cex=1)

require(reshape2)
melt_data = melt(boxplot_data)
write.table(melt_data, "new_data", quote=F, sep="\t")

save(list=ls(all=TRUE), file="all.RData")

RSCRIPT

close RCMD;

my $ret = system("R --restore --no-save < cmd.r");
if ( $ret ) {
        die "Error, died with $ret";
}


my %hs;
open IN,"sample_group";
while(<IN>){
	chomp;
	my @s = split /\t/, $_;
	$hs{$s[0]} = $s[1];
}
close IN;

open IN,"new_data";
open OUT,">new_data2";
my $header = <IN>;
chomp $header;
print OUT "\t$header\tGroup\n";
while (<IN>) {
	chomp;
	my @s = split /\t/, $_;
	print OUT "$_\t$hs{$s[2]}\n";
}
close IN;
close OUT;


open RCMD,">cmd.r";
print RCMD <<"RSCRIPT";

pdf("out.pdf", width=$width,height=$height)
load("all.RData")
library("ggplot2")
melt_data = read.table("new_data2", header=T, com='', sep="\t", row.names=1, check.names=F, quote="")

melt_data\$Group = as.character(melt_data\$Group)
melt_data\$Group = as.factor(melt_data\$Group)

#head(plot_data)
ggplot(melt_data, aes(x=Var2, y=value, colour =Group) ) + geom_violin() + ylab("log2(v+1)") + xlab("") + ggtitle("Violinplot ( analize v >= $min )")+ theme(panel.background = element_rect(fill = NA), panel.grid.major = element_line(colour = NA ), panel.grid.minor = element_line(colour = NA), panel.border = element_rect( colour = "gray30", fill = NA, size=1.5 ) ) + scale_colour_manual(values = sample_colors )
dev.off()

RSCRIPT
close RCMD;

}


my $ret = system("R --restore --no-save < cmd.r");
if ( $ret ) {
        die "Error, died with $ret";
}

move('out.pdf', $pdf);
unlink 'cmd.r';
unlink '__matrix__';

if ( defined($sample_group_file) ){
	unlink 'sample_group';
	unlink 'new_data';
	unlink 'new_data2';
	unlink 'all.RData';
}


get_html_link_v2($opts{html}, "Box plot Results", $pdf, "pdf");

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