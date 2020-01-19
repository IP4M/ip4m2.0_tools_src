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
GetOptions (\%opts, "input=s", 
"pp_method=s", 
"cutoff=s", 
"db=s", 
"html=s",
);

my $perl = $^X;


my $rungc_bin = "$Bin/../rungc/metams_runGC_v2.pl";
my $erah_bin = "$Bin/../erah/erah.pl";
my $gc_anno_bin = "$Bin/../gc_anno/lib_search.pl";
my $outlier_bin = "$Bin/../table_pretreatment/dataFilter.pl";
my $zf_bin = "$Bin/../table_pretreatment/dataReplace.pl";
my $norm_bin = "$Bin/../table_pretreatment/peakAreaNormal.pl";


my $dir = dirname $opts{html};
my $num = 1;
my @output_files;


mkdir "01_gc_ana";
if( $opts{pp_method} eq "metaMS"){
    my $cmd = "\"$perl\" \"$rungc_bin\"  -i  \"$opts{input}\"  -peaktable2 \"01_gc_ana/gcms_raw_pkTable.txt\" -peakspectra \"01_gc_ana/gcms_mass_spectra.msp\" -normpeakspectra \"01_gc_ana/gcms_mass_spectra_999norm.msp\" -html \"01_gc_ana/outputs.html\" -nSlaves \"2\" -param \"default\"";
    run_cmd($cmd);
    push @output_files,  (
        "01_gc_ana/gcms_raw_pkTable.txt", "tabH", "$num", 
        "01_gc_ana/gcms_mass_spectra.msp", "txt", "$num", 
        "01_gc_ana/gcms_mass_spectra_999norm.msp", "txt", "$num", 
        "01_gc_ana/TICs.pdf", "pdf",  "$num", 
        "01_gc_ana/BPCs.pdf", "pdf",  "$num", 
        "01_gc_ana/EICs", "dir",  "$num", 
    );
}else{
    my $cmd = "\"$perl\" \"$erah_bin\"  -i  \"$opts{input}\"  -analysis_time \"0\" -min_peak_width \"1\" -min_peak_height \"2500\" -noise_threshold \"500\" -avoid_processing_mz \"c(73:75,147:149)\" -min_spectra_cor \"0.9\" -max_time_dist \"3\" -mz_range \"c(70:600)\" -min_samples \"1\" -blocks_size \"NULL\" -peaktable \"01_gc_ana/gcms_raw_pkTable.txt\" -peakspectra \"01_gc_ana/gcms_mass_spectra.msp\" -html \"01_gc_ana/outputs.html\" -normpeakspectra   \"01_gc_ana/gcms_mass_spectra_999norm.msp\"";
    run_cmd($cmd);
    push @output_files,  (
        "01_gc_ana/gcms_raw_pkTable.txt", "tabH", "$num", 
        "01_gc_ana/gcms_mass_spectra.msp", "txt", "$num", 
        "01_gc_ana/gcms_mass_spectra_999norm.msp", "txt", "$num", 
        "01_gc_ana/TICs.pdf", "pdf",  "$num", 
        "01_gc_ana/BPCs.pdf", "pdf",  "$num", 
        "01_gc_ana/EICs", "dir",  "$num", 
    );
}
$num += 1;


mkdir "02_gc_ana";
my $cmd = "\"$perl\" \"$gc_anno_bin\" -matrix \"01_gc_ana/gcms_raw_pkTable.txt\" -query \"01_gc_ana/gcms_mass_spectra_999norm.msp\"  -output2 \"02_gc_ana/identified_uniq_pkTable.txt\" -method \"dot\" -cutoff \"$opts{cutoff}\" -output  \"02_gc_ana/identified_pkTable.txt\" -html \"02_gc_ana/outputs.html\" -mzRes \"0\" -rt_cut \"0.2\" -log \"02_gc_ana/detailed_information.txt\" -db \"$opts{db}\" -user \"no\"";
run_cmd($cmd);
push @output_files,  (
    "02_gc_ana/identified_pkTable.txt", "tabH", "$num",
    "02_gc_ana/identified_uniq_pkTable.txt", "tabH", "$num",
    "02_gc_ana/detailed_information.txt", "tabH", "$num",
);
$num += 1;


mkdir "03_gc_ana";
$cmd = "\"$perl\" \"$outlier_bin\" -html \"03_gc_ana/outputs.html\" -i \"02_gc_ana/identified_uniq_pkTable.txt\"  -o \"03_gc_ana/outliers_processed_pkTable.txt\" -c \"1.5\"";
run_cmd($cmd);
push @output_files,  (
    "03_gc_ana/outliers_processed_pkTable.txt", "tabH", "$num",
);
$num += 1;


mkdir "04_gc_ana";
$cmd = "\"$perl\" \"$zf_bin\" -html \"04_gc_ana/outputs.html\" -i \"03_gc_ana/outliers_processed_pkTable.txt\" -o \"04_gc_ana/zero_filled_pkTable.txt\" -m \"user\" -val \"0.00001\"";
run_cmd($cmd);
push @output_files,  (
    "04_gc_ana/zero_filled_pkTable.txt", "tabH", "$num",
);
$num += 1;



mkdir "05_gc_ana";
my $cmd = "\"$perl\" \"$norm_bin\" -html \"05_gc_ana/outputs.html\" -i \"04_gc_ana/zero_filled_pkTable.txt\" -c \"1000\" -o \"05_gc_ana/total_area_1000norm_pkTable.txt\"";
run_cmd($cmd);
push @output_files,  (
    "05_gc_ana/total_area_1000norm_pkTable.txt", "tabH", "$num",
);
$num += 1;






get_html_link_workflow($opts{html}, "GC-MS Data Preprocess Workflow Results", @output_files);








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