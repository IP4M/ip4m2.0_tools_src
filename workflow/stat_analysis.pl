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
GetOptions (\%opts, "matrix=s", "samples_file=s",  "html=s",
);

my $perl = $^X;


my $basic_sum_bin = "$Bin/../table_op/summary.pl";
my $ttest_bin = "$Bin/../diff/t_test.pl";
my $wtest_bin = "$Bin/../diff/wilcoxon_test.pl";
my $annova_bin = "$Bin/../diff/annova.pl";
my $kwtest_bin = "$Bin/../diff/kw_test.pl";
my $pca_bin = "$Bin/../stat/plot-pca.pl";
my $oplsda_bin = "$Bin/../diff/oplsda.pl";
my $plsda_bin = "$Bin/../diff/plsda.pl";
my $svm_bin = "$Bin/../diff/svm.pl";
my $rf_bin = "$Bin/../diff/randomforeast.pl";
my $biosigner_bin = "$Bin/../diff/biosigner.pl";
my $boruta_bin = "$Bin/../boruta/boruta.pl";


my $matrix = $opts{matrix};
my $samples_file = $opts{samples_file};
my $html = $opts{html};

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
deal_NA_data( $matrix, "__noNA_data_wfw__" );
$matrix = "__noNA_data_wfw__";


check_data_v3($matrix, $samples_file, "F");

my $dir = dirname $opts{html};
my $num = 1;
my @output_files;

mkdir "01_ana";
my $cmd = "\"$perl\" \"$basic_sum_bin\"  -i  \"$matrix\"  -o \"01_ana/pkTable_summary.txt\"  -html \"01_ana/outputs.html\" ";
run_cmd($cmd);
push @output_files,  (
    "01_ana/pkTable_summary.txt", "tabH", "$num",
);
$num += 1;

mkdir "02_ana";
$cmd = "\"$perl\" \"$ttest_bin\"  -i  \"$matrix\"  -o \"02_ana/t_test_results.txt\"  -g \"$samples_file\" -o2 \"02_ana/t_test_significant_results.txt\" -html \"02_ana/outputs.html\" -corrected \"fdr\" -paired F -var_equal F -thrN 0.05 ";
run_cmd($cmd);
push @output_files,  (
    "02_ana/t_test_results.txt", "tabH", "$num",
    "02_ana/t_test_significant_results.txt", "tabH", "$num",
);
$num += 1;

mkdir "03_ana";
$cmd = "\"$perl\" \"$wtest_bin\"  -i  \"$matrix\"  -o \"03_ana/wilcox_test_results.txt\"  -g \"$samples_file\" -o2 \"03_ana/wilcox_test_significant_results.txt\" -html \"03_ana/outputs.html\" -corrected \"fdr\" -paired F -tf T -thrN 0.05 ";
run_cmd($cmd);
push @output_files,  (
    "03_ana/wilcox_test_results.txt", "tabH", "$num",
    "03_ana/wilcox_test_significant_results.txt", "tabH", "$num",
);
$num += 1;

mkdir "04_ana";
$cmd = "\"$perl\" \"$annova_bin\"  -i  \"$matrix\"  -o \"04_ana/aov_results.txt\"  -g \"$samples_file\" -o2 \"04_ana/aov_significant_results.txt\" -html \"04_ana/outputs.html\" -corrected \"fdr\"  -thrN 0.05 ";
run_cmd($cmd);
push @output_files,  (
    "04_ana/aov_results.txt", "tabH", "$num",
    "04_ana/aov_significant_results.txt", "tabH", "$num",
);
$num += 1;

mkdir "05_ana";
$cmd = "\"$perl\" \"$kwtest_bin\"  -i  \"$matrix\"  -o \"05_ana/kw_test_results.txt\"  -g \"$samples_file\" -o2 \"05_ana/kw_test_significant_results.txt\" -html \"05_ana/outputs.html\" -corrected \"fdr\"  -thrN 0.05 ";
run_cmd($cmd);
push @output_files,  (
    "05_ana/kw_test_results.txt", "tabH", "$num",
    "05_ana/kw_test_significant_results.txt", "tabH", "$num",
);
$num += 1;

mkdir "06_ana";
$cmd = "\"$perl\" \"$pca_bin\"  -i  \"$matrix\"  -html \"06_ana/outputs.html\"  -g \"$samples_file\" -zscal FALSE";
run_cmd($cmd);
push @output_files,  (
    "06_ana/PCA_Score.txt", "tabH", "$num",
    "06_ana/PCA_R2.txt", "tabH", "$num",
    "06_ana/PCA_Screeplot.pdf", "pdf", "$num",
    "06_ana/PC12_Score_2D_Label.pdf", "pdf", "$num",
    "06_ana/PC12_Score_2D.pdf", "pdf", "$num",
    "06_ana/PC13_Score_2D_Label.pdf", "pdf", "$num",
    "06_ana/PC13_Score_2D.pdf", "pdf", "$num",
    "06_ana/PC23_Score_2D_Label.pdf", "pdf", "$num",
    "06_ana/PC23_Score_2D.pdf", "pdf", "$num",
);
$num += 1;

mkdir "plsda";
$cmd = "\"$perl\" \"$plsda_bin\"  -i  \"$matrix\"  -html \"plsda/outputs.html\"  -g \"$samples_file\"  -zscal FALSE";
run_cmd($cmd);
push @output_files,  (
    "plsda/PLSDA_Score.txt", "tabH", "$num",
    "plsda/PLSDA_R2X_R2Y_Q2.txt", "tabH", "$num",
    "plsda/PLSDA_Score_2D_Label.pdf", "pdf", "$num",
    "plsda/PLSDA_Score_2D.pdf", "pdf", "$num",
);
$num += 1;

mkdir "07_ana";
$cmd = "\"$perl\" \"$oplsda_bin\"  -i  \"$matrix\"  -html \"07_ana/outputs.html\"  -g \"$samples_file\" -vip 1 -zscal FALSE";
run_cmd($cmd);
push @output_files,  (
    "07_ana/OPLSDA_Score.txt", "tabH", "$num",
    "07_ana/OPLSDA_VIP.txt", "tabH", "$num",
    "07_ana/OPLSDA_VIP_Sig.txt", "tabH", "$num",
    "07_ana/OPLSDA_Permutation.txt", "tabH", "$num",
    "07_ana/Fitted_Curve_Parameter.txt", "tabH", "$num",
    "07_ana/OPLSDA_R2X_R2Y_Q2.txt", "tabH", "$num",
    "07_ana/OPLSDA_VPlot.pdf", "pdf", "$num",
    "07_ana/OPLSDA_Score_2D_Label.pdf", "pdf", "$num",
    "07_ana/OPLSDA_Score_2D.pdf", "pdf", "$num",
    "07_ana/OPLSDA_Permutation.pdf", "pdf", "$num",
    "07_ana/OPLSDA_R2X_R2Y_Q2.pdf", "pdf", "$num",
);
$num += 1;

mkdir "08_ana";
$cmd = "\"$perl\" \"$svm_bin\"  -input1  \"$matrix\"  -input2  \"$samples_file\"   -html \"08_ana/outputs.html\"  -scale FALSE  -kernel  linear -cost 1 -tolerance 0.001 -epsilon 0.1";
run_cmd($cmd);
push @output_files,  (
    "08_ana/Support_Vectors.txt", "tabH", "$num",
    "08_ana/SVM_Prediction.txt", "tabH", "$num",
    "08_ana/SVM_Prediction_Summary.txt", "relation", "$num",
    "08_ana/SVM_Imp_Rank.txt", "tabH", "$num",
    "08_ana/SVM_Imp.pdf", "pdf", "$num",
    "08_ana/SVM_Top10_Imp.pdf", "pdf", "$num",
);
$num += 1;


mkdir "09_ana";
$cmd = "\"$perl\" \"$rf_bin\"  -input1  \"$matrix\"  -input2  \"$samples_file\"   -html \"09_ana/outputs.html\" -scale FALSE -ntree 500 -replace TRUE ";
run_cmd($cmd);
push @output_files,  (
    "09_ana/RF_Prediction.txt", "tabH", "$num",
    "09_ana/RF_Prediction_Summary.txt", "relation", "$num",
    "09_ana/RF_Imp_Rank.txt", "tabH", "$num",
    "09_ana/RF_Imp.pdf", "pdf", "$num",
    "09_ana/RF_Top10_Imp.pdf", "pdf", "$num",
);
$num += 1;

mkdir "10_ana";
$cmd = "\"$perl\" \"$boruta_bin\"  -input1  \"$matrix\"  -input2  \"$samples_file\"   -html \"10_ana/outputs.html\" -scale FALSE -pValue 0.01 -maxRuns 100  ";
run_cmd($cmd);
copy("10_ana/Decision_Info.txt","10_ana/Boruta_Decision_Info.txt");
copy("10_ana/Decision_Boxplot.pdf","10_ana/Boruta_Decision_Boxplot.pdf");
push @output_files,  (
    "10_ana/Boruta_Decision_Info.txt", "tabH", "$num",
    "10_ana/Boruta_Decision_Boxplot.pdf", "pdf",  "$num",
);
$num += 1;

mkdir "11_ana";
$cmd = "\"$perl\" \"$biosigner_bin\"  -input1  \"$matrix\"  -input2  \"$samples_file\"  -html \"11_ana/outputs.html\" -methodVc all  -bootI 50 -pvalN 0.05 -permI 1  -tierC S -output1 \"11_ana/biosigner_variable_results.txt\" -output2 \"11_ana/biosigner_variable_significant_results.txt\"  -output4 \"11_ana/biosigner_figure-tier.pdf\"  -output5 \"11_ana/biosigner_figure-boxplot.pdf\" ";
run_cmd($cmd);
push @output_files,  (
    "11_ana/biosigner_variable_results.txt", "tabH", "$num",
    "11_ana/biosigner_variable_significant_results.txt", "tabH", "$num",
    "11_ana/biosigner_figure-tier.pdf", "pdf", "$num",
    "11_ana/biosigner_figure-boxplot.pdf", "pdf", "$num",
);
$num += 1;










get_html_link_workflow($opts{html}, "Statistical Analysis Results", @output_files);


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






##check matrix and sample group data, group number must be 2
sub check_data_v3{
	my $matrix = shift;
	my $sample_group = shift;
    my $paired = shift;
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
        if ( @gps != 2) {
            die "Error: groups number must be 2 in the sample group file!\n";
        }

        if( $paired eq "T"){
            if ( $cnt2{$gps[0]} != $cnt2{$gps[1]} ) {
                die "Error: for paired test, samples number must be equal in the two groups!\n";
            }
        }
	}
}


use Scalar::Util qw(looks_like_number);
sub isnumeric {  
        my $val = shift;  
        looks_like_number $val;
}  





