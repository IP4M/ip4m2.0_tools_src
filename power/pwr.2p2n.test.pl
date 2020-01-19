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
GetOptions (\%opts,"h=s",  "n1=s", "n2=s", "sig=s", "power=s", "alternative=s", "output=s", "p1=s", "p2=s", "html=s",
);

open RCMD, ">cmd.r";
if( defined($opts{h}) ){
print RCMD <<"RSCRIPT";

library(pwr)
res= pwr.2p2n.test(h=$opts{h},n1=$opts{n1},n2=$opts{n2},sig.level=$opts{sig},power=$opts{power},alternative="$opts{alternative}")
sink("results")
res
sink(file=NULL)

RSCRIPT
}else{
print RCMD <<"RSCRIPT";

library(pwr)
h=2*asin(sqrt($opts{p1})) - 2*asin(sqrt($opts{p2}))
res= pwr.2p2n.test(h=h,n1=$opts{n1},n2=$opts{n2},sig.level=$opts{sig},power=$opts{power},alternative="$opts{alternative}")
sink("results")
res
sink(file=NULL)

RSCRIPT
}
close RCMD;
&run_cmd("R --restore --no-save<cmd.r");





sub run_cmd {
        my $cmd = shift;
#        print "$cmd\n";

        my $ret = system( $cmd );
        if ( $ret ){
                die "Error, died with $ret";
        }

}



copy('results', $opts{output});
get_html_link_v2($opts{html}, "Power analysis Results", $opts{output}, "txt" );



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
                copy($file, "$out_dir/$name.html");
            }
    }

    my $basename = basename $out_dir;
	for (my $i = 0; $i <= $#par; $i+=2){
            my $file = $par[$i];
            my $name = basename $file;
	        print HTML "<li><a href=\"$basename/$name.html\" target=\"_parent\">$name</a></li>\n";
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
