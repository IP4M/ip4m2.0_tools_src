#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use File::Copy;
use File::Path;
use File::Basename;
use FindBin qw($Bin);
my $r_path = "$Bin/../../bin/R/bin";
$ENV{'PATH'} = "$r_path;$ENV{'PATH'}";

my $usage = <<__EOUSAGE__;

  perl expressPattern.pl [options]

  Required:
  --matrix      fpkm matrix

  Optional:

  --gene_dist   <string>      
  --gene_clust  <string>       
  --gene_cor    <string>
  --K           <int>          
  --Ktree       <int>               
  --Ptree       <float>             
  --width       <float>
  --height      <float>
  --output                    

__EOUSAGE__

my $matrix_file;
my $log = "log2";
my $gene_dist = "euclidean";
my $gene_clust = "complete";
my $gene_cor = "pearson";
my $Kmeans;
my $Ktree;
my $Ptree;
my $width = 7;
my $height = 7;

my $output;
my $output_directory;

&GetOptions (  
	'matrix=s' => \$matrix_file,
	'log=s' => \$log,
	'gene_dist=s' => \$gene_dist,
	'gene_clust=s' => \$gene_clust,
	'gene_cor=s' => \$gene_cor,
	'K=s' => \$Kmeans,
	'Ktree=s' => \$Ktree,
	'Ptree=s' => \$Ptree,
	'width=f' => \$width,
	'height=f' => \$height,
	'output=s' => \$output, 
    'out_dir=s' => \$output_directory,
);

 $width =  $width / 2.54;
 $height =  $height / 2.54;


if ( ! $matrix_file ) {
	die $usage;
}

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


&check_data( $matrix_file );

copy($matrix_file, '__matrix__');

open RCMD, ">cmd.r";
print RCMD <<"RSCRIPT";

data = read.table("__matrix__", header=T, com='', sep="\t", row.names=1, check.names=F, quote="")
data = as.matrix(data)
RSCRIPT

my $ylab = "v";
if ( $log eq "log2"){
        print RCMD "data = log2(data+1)\n";
	$ylab = "log2(v+1)";
} elsif ( $log eq "log2_center"){
        print RCMD "data = log2(data+1)\n";
        print RCMD "data = t(scale(t(data), scale=F))\n";
	$ylab = "centered log2(v+1)";
} elsif( $log eq "center" ){
        print RCMD "data = t(scale(t(data), scale=F))\n";
	$ylab = "centered v";
}

if ($gene_dist =~ /gene_cor/) {
        print RCMD "data = as.matrix(data)\n";
        print RCMD "gene_cor = cor(t(data), method=\"$gene_cor\", use='pairwise.complete.obs')\n";
        print RCMD "gene_dist = as.dist(1-gene_cor)\n";
} else {
        print RCMD "gene_dist = dist(data, method=\"$gene_dist\")\n";
}
print RCMD "if (nrow(data) <= 1) { message('Too few genes to analize'); quit(status=0); }\n";
print RCMD "hc_genes = hclust(gene_dist, method=\"$gene_clust\")\n";

if ($Kmeans) {
	print RCMD "kmeans_clustering <- kmeans(data, centers=$Kmeans, iter.max=100, nstart=5)\n";
	print RCMD "gene_partition_assignments = kmeans_clustering\$cluster\n";    
} elsif ( $Ktree ) {
	print RCMD "gene_partition_assignments <- cutree(as.hclust(hc_genes), k=$Ktree)\n";
} elsif ( $Ptree ) {
	print RCMD "gene_partition_assignments <- cutree(as.hclust(hc_genes), h=$Ptree/100*max(hc_genes\$height))\n";
}

print RCMD <<"RSCRIPT";
max_cluster_count = max(gene_partition_assignments)
gene_names = rownames(data)
num_cols = length(data[1,])
for (i in 1:max_cluster_count) {
	partition_i = (gene_partition_assignments == i)
	partition_data = data[partition_i,,drop=F]
	if (sum(partition_i) == 1) {
		dim(partition_data) = c(1,num_cols)
		colnames(partition_data) = colnames(data)
		rownames(partition_data) = gene_names[partition_i]
	}
	outfile = paste( "subcluster_", i, ".matrix.txt", sep='' )
	pdffile = paste( "subcluster_", i, ".pdf", sep='' )
	write.table(partition_data, file=outfile, quote=F, sep="\\t")
	pdf(pdffile, width=$width,height=$height)
	par(mar=c(14,6,6,4))
	ymin = min(partition_data); ymax = max(partition_data);
	plot(as.numeric(partition_data[1,]), type='l', ylim=c(ymin,ymax), col='lightgray', xaxt='n', xlab='', ylab='$ylab')
	axis(side=1, at=1:length(partition_data[1,]), labels=colnames(partition_data), las=2)
	if( length(partition_data[,1]) > 1 ){
		for(r in 2:length(partition_data[,1])) {
			points(as.numeric(partition_data[r,]), type='l', col='lightgray')
		}
	}
	points(as.numeric(colMeans(partition_data)), type='o', col='blue')
	dev.off()
}
RSCRIPT

close RCMD;
my $cmd = "R --vanilla -q < cmd.r";
my $ret = system($cmd);
if ($ret) {
	die "Error, cmd: $cmd died with ret $ret";
}
unlink '__matrix__';
unlink 'cmd.r';



my $dir = "./";
opendir(DIR, $dir) || die "Can't open directory $dir";
my @file_list;
foreach(readdir(DIR)){
        unless ( $_=~/^\./){
                if ( $_=~/^subcluster/ ){
                        push @file_list, $_;
                }
        }
}
close DIR;

@file_list = sort { $a cmp $b } @file_list;
my @out_file_list;


for (my $i = 0; $i <= $#file_list; $i++){
        my $input_file = "$dir/$file_list[$i]";
        my $output_file = "$output_directory/$file_list[$i]";
	if ( $file_list[$i]=~/txt$/ ){
		open IN,"$input_file";
		open OUT,">$output_file";
		my $h = <IN>;
		print OUT "\t$h";
		while (<IN>){
			print OUT "$_";
		}
		close IN;
		close OUT;
		unlink $input_file;
        push @out_file_list, ("$output_file", "tabH");
	} else {
		move($input_file, $output_file);
        push @out_file_list, ("$output_file", "pdf");
	}
        
}

get_html_link_v2($output, "Subcluster analysis Results", @out_file_list);

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