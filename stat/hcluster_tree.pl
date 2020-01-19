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
GetOptions (\%opts,"w=s",  "h=s", "i=s", "pdf=s", "log=s", "data=s",
 "sample_dist=s", "sample_clust=s", "sample_cor=s", "html=s",
);
my $usage = <<"USAGE";
       Program : $0
       Version : 1.0
       Discription: 
       Usage :perl $0 [options]
                   -i  input fpkm matrix

                   -w  width
                   -h  height

                   -pdf  output pdf
                   -data output nwk tree file

                   -log  log2 transfer/log2 center/center/none [none,log2,log2_center,center]

                   -sample_dist
                   -sample_clust
                   -sample_cor

USAGE



die $usage if ( !$opts{i} );

my $matrix_file = $opts{i};
my $pdf = $opts{pdf};
my $output = $opts{data};

my $width = defined($opts{w}) ? $opts{w} : 8;
my $height = defined($opts{h}) ? $opts{h} : 6;

 $width =  $width / 2.54;
 $height =  $height / 2.54;

my $log = defined($opts{log}) ? $opts{log} : "log2";

my $sample_dist = defined($opts{sample_dist}) ? $opts{sample_dist} : "euclidean";
my $sample_clust = defined($opts{sample_clust}) ? $opts{sample_clust} : "complete";
my $sample_cor = defined($opts{sample_cor}) ? $opts{sample_cor} : "pearson";

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



&check_data_v4( $matrix_file );

copy($matrix_file, '__matrix__');

open RCMD, ">cmd.r";
print RCMD <<"RSCRIPT";
library('ape')
pdf("out.pdf", width=$width,height=$height)
par(mar=c(3,2,2,12))
data = read.table("__matrix__", header=T, com='', sep="\t", row.names=1, check.names=F, quote="")
data = as.matrix(data)
RSCRIPT

if ( $log eq "log2"){
        print RCMD "data = log2(data+1)\n";
} elsif ( $log eq "log2_center"){
        print RCMD "data = log2(data+1)\n";
        print RCMD "data = t(scale(t(data), scale=F))\n";
} elsif( $log eq "center" ){
        print RCMD "data = t(scale(t(data), scale=F))\n";
}

if ($sample_dist =~ /sample_cor/) {
        print RCMD "sample_cor = cor(data, method=\"$sample_cor\", use='pairwise.complete.obs')\n";
        print RCMD "sample_dist = as.dist(1-sample_cor)\n";
} else {
        print RCMD "sample_dist = dist(t(data), method=\"$sample_dist\")\n";
}
print RCMD "hc_samples = hclust(sample_dist, method=\"$sample_clust\")\n";


print RCMD "
tree <-as.dendrogram(hc_samples)
plot(tree,type=\"rectangle\",horiz=TRUE)
tr <-as.phylo.hclust(hc_samples)
write.tree(tr,\"__nwk__\")
";


print RCMD "dev.off()\n";
close RCMD;
my $cmd = "R --vanilla -q < cmd.r";
my $ret = system($cmd);
if ($ret) {
    die "Error, cmd: $cmd died with ret $ret";
}

move('out.pdf', $pdf);
move('__nwk__', $output);
unlink 'cmd.r';
unlink '__matrix__';


get_html_link_v2($opts{html}, "Plot hcluster tree Results", $output, "txt", $pdf, "pdf",);



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