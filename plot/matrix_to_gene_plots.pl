#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Getopt::Long qw(:config no_ignore_case bundling pass_through);
use POSIX qw (floor ceil);
use FindBin qw($Bin);
use File::Copy;
my $r_path = "$Bin/../../bin/R/bin";
$ENV{'PATH'} = "$r_path;$ENV{'PATH'}";

my $usage = <<__EOUSAGE__;


########################################################################
#
# --matrix <string>              matrix file
#
# --pdf <string>                 output pdf file
#
########################################################################
#
# Dimensions:
#
#  --width_per_plot <float>      default: 2.5
#  --height_per_plot <float>     default: 2.5
#
# Layout:
#
#  --plots_per_row <int>         default: 2
#  --plots_per_col <int>         default: 3
#
# Misc:
#
#  --barplot
#  --log2
#
########################################################################

__EOUSAGE__

    ;



my $width_per_plot = 2.5;
my $height_per_plot = 2.5;

my $matrix;

my $plots_per_row = 2;
my $plots_per_col = 3;
my $help_flag;

my $pdf_out;

my $barplot_flag = 0;
my $log2_flag = 0;
my $html;

&GetOptions( 'h' => \$help_flag,
             
             'matrix=s' => \$matrix,
             'width_per_plot=f' => \$width_per_plot,
             'height_per_plot=f' => \$height_per_plot,
             'plots_per_row=i' => \$plots_per_row,
             'plots_per_col=i' => \$plots_per_col,
             'barplot' => \$barplot_flag,
             'log2_flag' => \$log2_flag,
             'pdf=s' => \$pdf_out,
             'html=s' => \$html,
    );





if ($help_flag) {
    die $usage;
}

unless ($matrix) {
    die $usage;
}

main: {

    
    my $R_script = "__tmp_plot_clusters.R";

    deal_NA_data( $matrix, "__noNA_data__" );
    $matrix = "__noNA_data__";

    check_data_v4($matrix);
    copy( $matrix, "__matrix__" );
    $matrix = "__matrix__";
    my $pdf_out_ = $pdf_out;
    $pdf_out = "__out.pdf__";

    open (my $ofh, ">$R_script") or die "Error, cannot write to $R_script";
   
    #print $ofh "postscript(file=\"my_cluster_plots.eps\", horizontal=FALSE, width=$width, height=$height, paper=\"special\")\n";
    print $ofh "pdf(file=\"$pdf_out\")\n";
    print $ofh "par(mfrow=c($plots_per_col, $plots_per_row))\n";
    #print $ofh "# png(file=\"my_cluster_plots.png\");\n"; 
    print $ofh "all_data = read.table(file=\"$matrix\", header=T, com=\'\', row.names=1, quote=\"\", check.names=F, sep=\"\\t\")\n";
    print $ofh "gene_names = rownames(all_data)\n";

    if ($log2_flag) {
        print $ofh "all_data = log2(all_data+1);\n";
    }
    
    print $ofh "for (i in 1:length(all_data[,1])) {\n";
    print $ofh "    data = all_data[i,]\n";
    print $ofh "    ymin = min(data); ymax = max(data);\n";
    if ($barplot_flag) {
        print $ofh "    barplot(as.numeric(data), cex.names=0.75, names.arg=colnames(data), las=2, main=gene_names[i])\n";
    }
    else {
        print $ofh "    plot(as.numeric(data), type='l', ylim=c(ymin,ymax), main=gene_names[i], col='blue', xaxt='n', xlab='', ylab='')\n";
        print $ofh "    axis(side=1, at=1:length(data), labels=colnames(all_data), las=2)\n";
    }
    print $ofh "}\n";
    print $ofh "dev.off()\n";
    
    close $ofh;
    

    &process_cmd("R --vanilla -q < $R_script");
    
    move( $pdf_out, $pdf_out_);
    unlink $R_script;

    &get_html_link_v2($html, "Metabolites plot Results", $pdf_out_ , "pdf");


    exit(0);
}


####
sub process_cmd {
    my ($cmd) = @_;

    print STDERR "CMD: $cmd\n";
    
    my $ret = system($cmd);
    if ($ret) {
        die "Error, cmd: $cmd died with ret $ret";
    }

    return;
}


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
