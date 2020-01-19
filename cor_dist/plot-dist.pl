#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use File::Copy;
use File::Basename;
use FindBin qw($Bin);
my $r_path = "$Bin/../../bin/R/bin";
$ENV{'PATH'} = "$r_path;$ENV{'PATH'}";

my %opts;
GetOptions (\%opts,"i=s","g=s","w=s","h=s","clc=f","rlc=f","col=s","lw=s","lh=s","mcr=s","cn=s", "o=s", "log=s", "html=s");
my $usage = <<"USAGE";
	Usage:perl $0 [options] 
		base options:
                -i	input dist file
                -o      output pdf file
                #-g	group file (a one column file telling what samples will be put in the figure)
			[
			Sample1
			Sample2
			Sample3
			...] 
		-w	the width of the figure,default:7
		-h	the height of the figure,default:8
                -clc	collab cex defalt:1.2
		-rlc	rowlab cex defalt:1.2
		-col	two or more colors to ramp. eg: black-red  or green-red-yellow and so on. colors Seprated by "-". default : darkblue-darkgreen-yellow-darkred
		#-lw	split plot in width,defalt (four numbers) :0.1:0.2:4:1
                #-lh	split plot in heigth,defalt (three numbers) :0.3:5.5:1.2
		#-mcr	the margins for column and row names,default: 5:7 
		-cn	add cellnote or not [yes/no],default:yes
                -log 
USAGE
#die $usage if (!($opts{i}&&$opts{g}));
die $usage if (!($opts{i}));

$opts{g}=defined$opts{g}?$opts{g}:"all";
$opts{w}=defined$opts{w}?$opts{w}:7;
$opts{h}=defined$opts{h}?$opts{h}:8;
$opts{clc}=defined$opts{clc}?$opts{clc}:1.2;
$opts{rlc}=defined$opts{rlc}?$opts{rlc}:1.2;
$opts{col}=defined$opts{col}?$opts{col}:"darkblue-darkgreen-yellow-darkred";
$opts{lw}=defined$opts{lw}?$opts{lw}:"0.1:0.2:4:1";
$opts{lh}=defined$opts{lh}?$opts{lh}:"0.3:5.5:1.2";
$opts{mcr}=defined$opts{mcr}?$opts{mcr}:"10:30";
$opts{cn}=defined$opts{cn}?$opts{cn}:"yes";
$opts{log}=defined$opts{log}?$opts{log}:"none";

$opts{w} = $opts{w} / 2.54;
$opts{h} = $opts{h} / 2.54;

my $lib_path = "$Bin/heatmap.d.r";
copy( $lib_path, "./__libs.r__" );

copy( $opts{i}, "__input_dm__");
$opts{i} = "__input_dm__";

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
deal_NA_data( $opts{i}, "__noNA_data__" );
$opts{i} = "__noNA_data__";


check_data($opts{i});

open RCMD, ">cmd.r";

print RCMD "
library(\"gplots\")
library(\"gtools\")
tab<-read.delim(file=\"$opts{i}\",header=TRUE,sep=\"\\t\", row.names=1)
tab<-as.matrix(tab) 
hmc<-tab

if (\"$opts{g}\"!=\"all\"){
group <- as.character(read.table(\"$opts{g}\",header=FALSE)[,1])
hmc<-tab[group,group]
}
";

if ( $opts{log} eq "log2"){
        print RCMD "hmc = log2(hmc+1)\n";
}

print RCMD "
#heatmap_color
kbn<-200
cramp <-\"$opts{col}\"
colramp <-unlist(strsplit(cramp,\"-\"))
ramp <-colorRamp(as.vector(colramp))
heatcol <- paste(rgb(ramp(seq(0,1,length=kbn)),max=255),\"FF\",sep=\"\")

##

pdf(\"__out.pdf__\",width=as.numeric(\"$opts{w}\"),height=as.numeric(\"$opts{h}\"))
#layout(matrix(1:3, 3, 1), heights =c(1,1,1))
#par(mfrow=c(3,1))
#par(mar=c(14,6,6,14))

l_w <-as.numeric(unlist(strsplit(\"$opts{lw}\",\":\")))
l_h <-as.numeric(unlist(strsplit(\"$opts{lh}\",\":\")))
mcr <-as.numeric(unlist(strsplit(\"$opts{mcr}\",\":\")))

source(\"./__libs.r__\")

if (\"$opts{cn}\"==\"yes\"){
heatmap.d(hmc,cellnote=round(hmc,3),margins = mcr ,cexCol=\"$opts{clc}\",cexRow=\"$opts{rlc}\",Rowv=FALSE, symm=TRUE,trace=\"none\",density.info=\"none\",dendrogram=\"none\",col=heatcol,keysize=1.5,notecol=\"white\",lmat=rbind( c(0,0, 3,0), c(2,0,1,1), c(0,4,4,0) ), lhei=l_h,lwid=l_w)
}

if (\"$opts{cn}\"==\"no\"){
heatmap.d(hmc,margins = mcr ,cexCol=\"$opts{clc}\",cexRow=\"$opts{rlc}\",Rowv=FALSE, symm=TRUE,trace=\"none\",density.info=\"none\",dendrogram=\"none\",col=heatcol,keysize=1.5,notecol=\"white\",lmat=rbind( c(0,0, 3,0), c(2,0,1,1), c(0,4,4,0) ), lhei=l_h,lwid=l_w)
}
dev.off()



";

my $ret = system("R --restore --no-save < cmd.r");
if ( $ret ) {
        die "Error, died with $ret";
}

move( '__out.pdf__', $opts{o});

unlink 'cmd.r';
unlink '__libs.r__';
unlink '__input_dm__';
unlink 'Rplots.pdf' if -e "Rplots.pdf";


get_html_link_v2($opts{html}, "Distance matrix heatmap plot Results", $opts{o} , "pdf",
);


#####
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