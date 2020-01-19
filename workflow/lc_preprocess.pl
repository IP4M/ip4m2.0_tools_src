#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use File::Copy;
use FindBin qw($Bin);
use Data::Dumper;
use File::Basename;
use Scalar::Util qw(looks_like_number);
use File::Copy::Recursive qw(dircopy);

my %opts;
GetOptions (\%opts, "input=s",    "html=s", 
"db=s",
"ions=s", 
"mz_cut=s", 
);

my $perl = $^X;


my $runlc_bin = "$Bin/../runlc/metams_runLC.pl";
my $lc_anno_bin = "$Bin/../lc_anno/lib_search_ppm.pl";
my $outlier_bin = "$Bin/../table_pretreatment/dataFilter.pl";
my $zf_bin = "$Bin/../table_pretreatment/dataReplace.pl";
my $norm_bin = "$Bin/../table_pretreatment/peakAreaNormal.pl";


my $dir = dirname $opts{html};
my $num = 1;
my @output_files;


mkdir "01_lc_ana";
my $cmd = "\"$perl\" \"$runlc_bin\" -i \"$opts{input}\" -html \"01_lc_ana/outputs.html\" -output1 \"01_lc_ana/lcms_raw_pkTable.txt\" -output2 \"01_lc_ana/lcms_isotopes_removed_pkTable.txt\" -nSlaves \"2\"   -settings  \"user\"         -peakpicking \"matchedFilter\" -step \"0.1\" -fwhm \"20\" -max \"50\" -snthresh \"4\" -minfrac \"0.3\" -minsamp \"3\" -mzwid \"0.1\" -bws \"130,10\" -missingratio \"0.2\" -extraratio \"0.1\" -retcor \"linear\" -family \"symmetric\"            -fillPeaks \"TRUE\"  -camera_perfwhm \"0.6\"  -camera_cor_eic_th \"0.7\"   -camera_ppm \"5\"  -intensity \"into\"   -polarity \"positive\"";
run_cmd($cmd);
push @output_files,  (
    "01_lc_ana/lcms_raw_pkTable.txt", "tabH", "$num",
    "01_lc_ana/lcms_isotopes_removed_pkTable.txt", "tabH", "$num",
    "01_lc_ana/raw_tics.pdf", "pdf", "$num",
    "01_lc_ana/raw_bpcs.pdf", "pdf", "$num",
    "01_lc_ana/rtcorrected_tics.pdf", "pdf", "$num",
    "01_lc_ana/rtcorrected_bpcs.pdf", "pdf", "$num",
    "01_lc_ana/rt_deviation_plot.pdf", "pdf", "$num",
    "01_lc_ana/EICs", "dir", "$num",
);
$num += 1;


mkdir "02_lc_ana";
my $cmd = "\"$perl\" \"$lc_anno_bin\" -matrix \"01_lc_ana/lcms_isotopes_removed_pkTable.txt\" -html \"02_lc_ana/outputs.html\" -output \"02_lc_ana/identified_pkTable.txt\" -output2 \"02_lc_ana/identified_uniq_pkTable.txt\" -log \"02_lc_ana/detailed_information.txt\" -db \"$opts{db}\" -rt_cut \"0.2\" -ions \"$opts{ions}\" -mz_cut \"$opts{mz_cut}\"";
run_cmd($cmd);
push @output_files,  (
    "02_lc_ana/identified_pkTable.txt", "tabH", "$num",
    "02_lc_ana/identified_uniq_pkTable.txt", "tabH", "$num",
    "02_lc_ana/detailed_information.txt", "tabH", "$num",
);
$num += 1;


mkdir "03_lc_ana";
$cmd = "\"$perl\" \"$outlier_bin\" -html \"03_lc_ana/outputs.html\" -i \"02_lc_ana/identified_uniq_pkTable.txt\"  -o \"03_lc_ana/outliers_processed_pkTable.txt\" -c \"1.5\"";
run_cmd($cmd);
push @output_files,  (
    "03_lc_ana/outliers_processed_pkTable.txt", "tabH", "$num",
);
$num += 1;


mkdir "04_lc_ana";
$cmd = "\"$perl\" \"$zf_bin\" -html \"04_lc_ana/outputs.html\" -i \"03_lc_ana/outliers_processed_pkTable.txt\" -o \"04_lc_ana/zero_filled_pkTable.txt\" -m \"user\" -val \"0.00001\"";
run_cmd($cmd);
push @output_files,  (
    "04_lc_ana/zero_filled_pkTable.txt", "tabH", "$num",
);
$num += 1;



mkdir "05_lc_ana";
my $cmd = "\"$perl\" \"$norm_bin\" -html \"05_lc_ana/outputs.html\" -i \"04_lc_ana/zero_filled_pkTable.txt\" -c \"1000\" -o \"05_lc_ana/total_area_1000norm_pkTable.txt\"";
run_cmd($cmd);
push @output_files,  (
    "05_lc_ana/total_area_1000norm_pkTable.txt", "tabH", "$num",
);
$num += 1;




get_html_link_workflow($opts{html}, "LC-MS Data Preprocess Workflow Results", @output_files);





sub run_cmd {
        my $cmd = shift;
        print "$cmd\n";
        my $ret = system( $cmd );
        if ( $ret ){
                die "Error, died with $ret";
        }
}










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




sub get_html_link_workflow{
	my @param = @_;
	my $html = shift @param;
	my $title = shift @param;
	my @par = @param;

    my $out_dir = "$html.files";
    mkdir $out_dir;
    my $html_dir = dirname($html);

	open HTML, ">$html";
	print HTML "<html><head><title>$title</title></head><body><h3>Output Files:</h3><p><ul>\n";

    for(my $i = 0; $i <= $#par; $i+=3){
            my $file = $par[$i];
            my $type = $par[$i+1];
            my $num = $par[$i+2];
            my $name = basename $file;
            $name = "$num.$name";
            if( $type eq "tabH" ){
                copy($file, "$html_dir/$name");
                get_html_table($file, "$out_dir/$name.html", "y");
            }elsif(  $type eq "tabN" ){
                copy($file, "$html_dir/$name");
                get_html_table($file, "$out_dir/$name.html", "n");
            }elsif(  $type eq "txt" ){
                copy($file, "$html_dir/$name");
                txt2html($file, "$out_dir/$name.html");
            }elsif(  $type eq "pdf"){
                copy($file, "$html_dir/$name");
            }elsif(  $type eq "dir" ){
                dircopy(  $file , "$html_dir/$name") or die $!;
            }elsif(  $type eq "relation" ){
                copy($file, "$html_dir/$name");
                get_relation_table($file, "$out_dir/$name.html", "y");
            }
            else{
                copy($file, "$out_dir/$name");
            }
    }

    my $basename = basename $out_dir;
	for (my $i = 0; $i <= $#par; $i+=3){
            my $file = $par[$i];
            my $type = $par[$i+1];
            my $num = $par[$i+2];
            my $name = basename $file;
            $name = "$num.$name";
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

