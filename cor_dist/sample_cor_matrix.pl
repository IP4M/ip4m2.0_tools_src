#!/usr/bin/perl
use strict;
use warnings;
use File::Copy;
use Getopt::Long;
use File::Copy;
use FindBin qw($Bin);
my $r_path = "$Bin/../../bin/R/bin";
$ENV{'PATH'} = "$r_path;$ENV{'PATH'}";

my %opts;
GetOptions (\%opts,"w=s",  "h=s", "i=s", "s=s", "min=s", "pdf=s", "log=s","method=s", "color=s", "clust=s", "data=s", "note=s", "rcex=s", "ccex=s", "html=s", "output3=s", "output4=s", "corrected=s",
);
my $usage = <<"USAGE";
       Program : $0
       Version : 1.0
       Discription: 
       Usage :perl $0 [options]
                   -i  input fpkm matrix
                   -s  sample group file
                   -w  width
                   -h  height
                   -min  min row sum
                   -pdf  output pdf
                   -data touput data
                   -log  log2 transfer
                   -method   cor method
                   -color   headtmap color 
                   -clust   clust method
                   -note    note number
                   -rcex    row font size 
                   -ccex    col font size
USAGE

die $usage if ( !$opts{i} );

my $matrix_file = $opts{i};
my $sample_group_file = $opts{s};
my $width = defined($opts{w}) ? $opts{w} : 6;
my $height = defined($opts{h}) ? $opts{h} : 6;

my $min = defined($opts{min}) ? $opts{min} : 0;
my $log = defined($opts{log}) ? $opts{log} : "F";
my $note = defined($opts{note}) ? $opts{note} : "F";
my $method = defined($opts{method}) ? $opts{method} : "pearson";
my $color = defined($opts{color}) ? $opts{color} : "purple,black,yellow";
my $clust = defined($opts{clust}) ? $opts{clust} : "complete";
my $ccex = defined($opts{ccex}) ? $opts{ccex} : 1;
my $rcex = defined($opts{rcex}) ? $opts{rcex} : 1;

my $pdf = $opts{pdf};
my $output = $opts{data};

&check_data_v4( $matrix_file, $sample_group_file );

my @colors = split(/,/, $color);
foreach my $color (@colors) {
	$color = "\'$color\'";
}
my $HEATMAP_COLORS = join(",", @colors);
my $r_script = "$Bin/sample_cor_matrix.R";

copy($matrix_file, '__matrix__');
copy($r_script, '__tmp.R');
if ( defined($sample_group_file) ){
        copy($sample_group_file, 'sample_group');
}
if ( $note eq "T" ){
	$note = "cellnote=round(sample_cor,3),";
} else {
	$note = "";
}

my $dend = 'both';
if( $clust eq  "none"){
    $dend = "none";
    $clust = "complete";
}


open RCMD, ">cmd.r";

unless ( defined($sample_group_file) ){
print RCMD <<"RSCRIPT";
source('__tmp.R');
myheatcol = colorpanel(75, $HEATMAP_COLORS)
pdf("out.pdf", width=$width,height=$height)
data = read.table("__matrix__", header=T, com='', sep="\t", row.names=1, check.names=F,quote="")
data = as.matrix(data)
dend = "$dend"

RSCRIPT
if ( $min >= 0 ){
	print RCMD "data = data[rowSums(data)>=$min,]\n";
}
if ( $log eq "T" ){
	print RCMD "data = log2(data+1)\n";
}
print RCMD <<"RSCRIPT";
library(Hmisc)
rcorr = rcorr(data, type="$method" )
#sample_cor = cor(data, method="$method", use='pairwise.complete.obs')
sample_cor = rcorr\$r
p_value = rcorr\$P
q_value = p.adjust(as.vector(p_value), method = "$opts{corrected}")
q_value = matrix(q_value, nrow=nrow(p_value), dimnames=dimnames(p_value))
write.table(p_value, "__p_matrix__", quote=F, sep="\t")
write.table(q_value, "__q_matrix__", quote=F, sep="\t")
write.table(sample_cor, "__cor_matrix__", quote=F, sep="\t")
sample_dist = as.dist(1-sample_cor)
hc_samples = hclust(sample_dist, method="$clust")
heatmap.3(sample_cor, $note dendrogram=dend, Rowv=as.dendrogram(hc_samples), Colv=as.dendrogram(hc_samples), col=myheatcol , scale='none', symm=TRUE, key=TRUE,density.info='none', trace='none', symkey=FALSE, symbreaks=F, margins=c(10,10), cexCol=$ccex, cexRow=$rcex )


dev.off()

RSCRIPT
} else {
print RCMD <<"RSCRIPT";
source('__tmp.R');
myheatcol = colorpanel(75, $HEATMAP_COLORS)
pdf("out.pdf", width=$width,height=$height)
data = read.table("__matrix__", header=T, com='', sep="\t", row.names=1, check.names=F, quote="")
data = as.matrix(data)
dend = "$dend"

samples_data = read.table("sample_group", header=F, check.names=F, quote="", com='', colClasses=c("character","character"))
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


RSCRIPT
if ( $min >= 0 ){
        print RCMD "data = data[rowSums(data)>=$min,]\n";
}
if ( $log eq "T" ){
        print RCMD "data = log2(data+1)\n";
}
print RCMD <<"RSCRIPT";
sample_factoring = colnames(data)
for (i in 1:nsamples) {
    sample_type = sample_types[i]
    replicates_want = sample_type_list[[sample_type]]
    sample_factoring[ colnames(data) %in% replicates_want ] = sample_type
}
sampleAnnotations = matrix(ncol=ncol(data),nrow=nsamples)
for (i in 1:nsamples) {
  sampleAnnotations[i,] = colnames(data) %in% sample_type_list[[sample_types[i]]]
}
sampleAnnotations = apply(sampleAnnotations, 1:2, function(x) as.logical(x))
sampleAnnotations = sample_matrix_to_color_assignments(sampleAnnotations, col=sample_colors)
rownames(sampleAnnotations) = as.vector(sample_types)
colnames(sampleAnnotations) = colnames(data)

library(Hmisc)
rcorr = rcorr(data, type="$method" )
#sample_cor = cor(data, method="$method", use='pairwise.complete.obs')
sample_cor = rcorr\$r
p_value = rcorr\$P
q_value = p.adjust(as.vector(p_value), method = "$opts{corrected}")
q_value = matrix(q_value, nrow=nrow(p_value), dimnames=dimnames(p_value))
write.table(p_value, "__p_matrix__", quote=F, sep="\t")
write.table(q_value, "__q_matrix__", quote=F, sep="\t")
write.table(sample_cor, "__cor_matrix__", quote=F, sep="\t")
sample_dist = as.dist(1-sample_cor)
hc_samples = hclust(sample_dist, method="$clust")
heatmap.3(sample_cor, $note dendrogram=dend, Rowv=as.dendrogram(hc_samples), Colv=as.dendrogram(hc_samples), col=myheatcol , scale='none', symm=TRUE, key=TRUE,density.info='none', trace='none', symkey=FALSE, symbreaks=F, margins=c(10,10), cexCol=$ccex, cexRow=$rcex, ColSideColors=sampleAnnotations, RowSideColors=t(sampleAnnotations))

dev.off()

RSCRIPT

}


close RCMD;

my $cmd = "R --vanilla -q < cmd.r";
my $ret = system($cmd);
if ($ret) {
    die "Error, cmd: $cmd died with ret $ret";
}

open IN,"__cor_matrix__";
open OUT,">$output";
my $h = <IN>;
print OUT "\t$h";
while (<IN>){
        print OUT "$_";
}
close IN;
close OUT;

move('out.pdf', $pdf);
unlink '__tmp.R';
unlink 'cmd.r';
unlink '__matrix__';
unlink '__cor_matrix__';
if ( defined($sample_group_file) ){
        unlink 'sample_group';
}

getRTable("__p_matrix__", $opts{output3});
getRTable("__q_matrix__", $opts{output4});


get_html_link_v2($opts{html}, "Correlation analysis Results", "$pdf", "pdf", $output , "relation", $opts{output3}, "relation", $opts{output4}, "relation",
);

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