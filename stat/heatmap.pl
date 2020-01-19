#!/usr/bin/perl
use strict;
use warnings;
use File::Copy;
use Getopt::Long;
use File::Copy;
use FindBin qw($Bin);
use File::Basename;
my $r_path = "$Bin/../../bin/R/bin";
$ENV{'PATH'} = "$r_path;$ENV{'PATH'}";

my %opts;
GetOptions (\%opts,"w=s",  "h=s", "i=s", "s=s", "pdf=s", "log=s", "color=s", "rcex=s", "ccex=s", "data=s",
"gene_dist=s", "sample_dist=s", "gene_clust=s", "sample_clust=s", "gene_cor=s", "sample_cor=s", "html=s", "dendrogram=s",
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

                   -pdf  output pdf
                   -data matrix for heatmap

                   -log  log2 transfer/log2 center/center/none [none,log2,log2_center,center]
                   -color   headtmap color 
                   -rcex    row font size 
                   -ccex    col font size

                   -gene_dist
                   -sample_dist
                   -gene_clust
                   -sample_clust
                   -gene_cor
                   -sample_cor

USAGE

die $usage if ( !$opts{i} );

my $matrix_file = $opts{i};
my $sample_group_file = $opts{s};
my $pdf = $opts{pdf};
my $output = $opts{data};

my $width = defined($opts{w}) ? $opts{w} : 7;
my $height = defined($opts{h}) ? $opts{h} : 10;

 $width =  $width / 2.54;
 $height =  $height / 2.54;

my $log = defined($opts{log}) ? $opts{log} : "log2";
my $color = defined($opts{color}) ? $opts{color} : "purple,black,yellow";
my $ccex = defined($opts{ccex}) ? $opts{ccex} : 1;
my $rcex = defined($opts{rcex}) ? $opts{rcex} : 1;

my $gene_dist = defined($opts{gene_dist}) ? $opts{gene_dist} : "euclidean";
my $sample_dist = defined($opts{sample_dist}) ? $opts{sample_dist} : "euclidean";
my $gene_clust = defined($opts{gene_clust}) ? $opts{gene_clust} : "complete";
my $sample_clust = defined($opts{sample_clust}) ? $opts{sample_clust} : "complete";
my $gene_cor = defined($opts{gene_cor}) ? $opts{gene_cor} : "pearson";
my $sample_cor = defined($opts{sample_cor}) ? $opts{sample_cor} : "pearson";
my $dendrogram = defined($opts{dendrogram}) ? $opts{dendrogram} : "both";

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

my @colors = split(/,/, $color);
foreach my $color (@colors) {
	$color = "\'$color\'";
}
my $HEATMAP_COLORS = join(",", @colors);

my $r_script = "$Bin/heatmap.R";
copy($matrix_file, '__matrix__');
copy($r_script, '__tmp.R');
if ( defined($sample_group_file) ){
	copy($sample_group_file, 'sample_group');
}


open RCMD, ">cmd.r";
print RCMD <<"RSCRIPT";
source('__tmp.R');
myheatcol = colorpanel(75, $HEATMAP_COLORS)
pdf("out.pdf", width=$width,height=$height)
data = read.table("__matrix__", header=T, com='', sep="\t", row.names=1, check.names=F, quote="")
data = as.matrix(data)
RSCRIPT

if ( defined($sample_group_file) ){
	print RCMD <<"RSCRIPT";
samples_data = read.table("sample_group", header=F, check.names=F, quote="", com='', colClasses=c("character","character") )
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
}


if ( $log eq "log2"){
	print RCMD "data = log2(data+1)\n";
} elsif ( $log eq "log2_center"){
	print RCMD "data = log2(data+1)\n";
	print RCMD "data = t(scale(t(data), scale=F))\n";
} elsif( $log eq "center" ){
	print RCMD "data = t(scale(t(data), scale=F))\n";
}
print RCMD "write.table(data, \"__heat_matrix__\", quote=F, sep=\"\t\")\n";

if ( defined( $sample_group_file ) ){
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
RSCRIPT
}


if ($gene_dist =~ /gene_cor/) {
	print RCMD "data = as.matrix(data)\n";
	print RCMD "gene_cor = cor(t(data), method=\"$gene_cor\", use='pairwise.complete.obs')\n";
	print RCMD "gene_dist = as.dist(1-gene_cor)\n";
} else {
	print RCMD "gene_dist = dist(data, method=\"$gene_dist\")\n";
}
print RCMD "if (nrow(data) <= 1) { message('Too few genes to generate heatmap'); quit(status=0); }\n";
print RCMD "hc_genes = hclust(gene_dist, method=\"$gene_clust\")\n";
if ($sample_dist =~ /sample_cor/) {
	print RCMD "sample_cor = cor(data, method=\"$sample_cor\", use='pairwise.complete.obs')\n";
	print RCMD "sample_dist = as.dist(1-sample_cor)\n";
} else {
	print RCMD "sample_dist = dist(t(data), method=\"$sample_dist\")\n";
}
print RCMD "hc_samples = hclust(sample_dist, method=\"$sample_clust\")\n";

if ( defined($sample_group_file) ){
	print RCMD <<"RSCRIPT";
heatmap.3( data, dendrogram="$dendrogram", Rowv=as.dendrogram(hc_genes), Colv=as.dendrogram(hc_samples), col=myheatcol, scale="none", density.info="none", trace="none", key=TRUE, keysize=0.8, cexCol=$ccex, cexRow=$rcex, margins=c(10,10), ColSideColors=sampleAnnotations)
RSCRIPT
} else {
	print RCMD <<"RSCRIPT";
heatmap.3( data, dendrogram="$dendrogram", Rowv=as.dendrogram(hc_genes), Colv=as.dendrogram(hc_samples), col=myheatcol, scale="none", density.info="none", trace="none", key=TRUE, keysize=0.8, cexCol=$ccex, cexRow=$rcex, margins=c(10,10) )
RSCRIPT
}



print RCMD "dev.off()\n";
close RCMD;
my $cmd = "R --vanilla -q < cmd.r";
my $ret = system($cmd);
if ($ret) {
    die "Error, cmd: $cmd died with ret $ret";
}

open IN,"__heat_matrix__";
open OUT,">$output";
my $h = <IN>;
print OUT "\t$h";
while (<IN>){
        print OUT "$_";
}
close IN;
close OUT;

move('out.pdf', $pdf);
unlink '__heat_matrix__';
unlink '__tmp.R';
unlink 'cmd.r';
unlink '__matrix__';
if ( defined($sample_group_file) ){
	unlink 'sample_group';
}


get_html_link_v2($opts{html}, "Plot heatmap with tree", "$pdf", "pdf", $output , "tabH",
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