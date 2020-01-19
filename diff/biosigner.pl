#!/usr/bin/perl -w
use strict;
use warnings;

use Getopt::Long;
use File::Copy;
use FindBin qw($Bin);
my $r_path = "$Bin/../../bin/R/bin";
$ENV{'PATH'} = "$r_path;$ENV{'PATH'}";

my %opts;
GetOptions (\%opts, "input1=s", "input2=s", "methodVc=s", "bootI=s", "pvalN=s", "permI=s", "tierC=s", 
"parCexN=s", "output1=s", "output2=s",  "output4=s", "html=s", "output5=s",
);

my $input1 = $opts{input1};
my $input2 = $opts{input2};
my $methodVc = $opts{methodVc};
my $bootI = $opts{bootI};
my $pvalN = $opts{pvalN};
my $permI = $opts{permI};
my $tierC = $opts{tierC};

my $output1 = $opts{output1};
my $output2 = $opts{output2};

my $output4 = $opts{output4};
my $output5 = $opts{output5};
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
deal_NA_data( $input1, "__noNA_data__" );
$input1 = "__noNA_data__";



check_data_v3($input1, $input2, "F");
copy($input1, '__matrix__');
copy($input2, '__group__');
open RCMD, ">cmd.r";
print RCMD <<"RSCRIPT";
library(biosigner)

xMN <- read.table(quote="","__matrix__", check.names = FALSE, header = TRUE, row.names = 1, sep = "\t", com='')
tmpDf = xMN
samDF <- read.table("__group__", header = F, sep = "\t",check.names = FALSE,quote="",  row.names = 1,com = '', colClasses=c("character","character"))
colnames(samDF) = c("Group")

samDF\$Group = as.character(samDF\$Group)
samDF\$Group = as.factor(samDF\$Group)

std <- function(x) sd(x) / sqrt(length(x))
uniq.group <- as.character(unique(samDF\$Group))
group1 <- rownames(subset(samDF, Group == uniq.group[1]))
group2 <- rownames(subset(samDF, Group == uniq.group[2]))
group1 <- as.character(group1)
group2 <- as.character(group2)
data = tmpDf
mean1 = apply(data[, group1], 1, mean)
okk = as.data.frame(mean1)
rownames(okk) = rownames(data)
okk\$var1 = apply(data[, group1], 1, var)
okk\$std1 = apply(data[, group1], 1, std)
okk\$mean2 = apply(data[, group2], 1, mean)
okk\$var2 = apply(data[, group2], 1, var)
okk\$std2 = apply(data[, group2], 1, std)
okk\$logFC = log2(okk\$mean2 / okk\$mean1)
okk\$logFC = log2(okk\$mean2 / okk\$mean1)

names(okk)[names(okk) == "mean1"] = paste("mean(", uniq.group[1], ")", sep = "")
names(okk)[names(okk) == "var1"] = paste("variance(", uniq.group[1], ")", sep = "")
names(okk)[names(okk) == "std1"] = paste("stderr(", uniq.group[1], ")", sep = "")
names(okk)[names(okk) == "mean2"] = paste("mean(", uniq.group[2], ")", sep = "")
names(okk)[names(okk) == "var2"] = paste("variance(", uniq.group[2], ")", sep = "")
names(okk)[names(okk) == "std2"] = paste("stderr(", uniq.group[2], ")", sep = "")
names(okk)[names(okk) == "logFC"] = paste("log2FC(", uniq.group[2], "/", uniq.group[1], ")", sep = "")

xMN = xMN[rownames(samDF)]
xMN = t(xMN)

respC = "Group"
tierMaxC <- "$tierC"
pvalN <- $pvalN
methodVc = "$methodVc"
bootI = $bootI
permI = $permI


respVc <- samDF[, respC]
respFc <- factor(respVc)
set.seed(123)
bsnLs <- biosign(x = xMN, y = respFc, methodVc = methodVc, bootI = bootI,  pvalN = pvalN, permI=permI, printL = FALSE, plotL = FALSE)
set.seed(NULL)


tierMC <- bsnLs\@tierMC

plot(bsnLs, tierMaxC = tierMaxC, fig.pdfC = "tier.pdf")
plot(bsnLs, tierMaxC = tierMaxC, typeC = "boxplot", fig.pdfC = "boxplot.pdf")


allDF = cbind( okk[rownames(tierMC) ,], as.data.frame(tierMC))
write.table(allDF,"out.all",quote=FALSE ,sep="\\t")

tierFullVc <- c("S", LETTERS[1 : 5])
tierVc <- tierFullVc[1 : which(tierFullVc == tierMaxC)]
    tierDF <- data.frame(tier = sapply(rownames(tmpDf),
    function(varC) {
        varTirVc <- tierMC[varC,]
        names(varTirVc) = colnames(tierMC)
        varTirVc <- names(varTirVc)[varTirVc %in% tierVc]
        paste(varTirVc, collapse = "|")
    }),
    stringsAsFactors = FALSE)
    colnames(tierDF) <- paste(
    colnames(tierDF),
    paste(tierVc, collapse = "+"),
    sep = "_")



outDf = cbind(okk, tierDF )
outDf = outDf[ rownames(tierMC),]
lastName = colnames(tierDF)
outDf = outDf[outDf[[lastName]] != "",]
write.table(outDf,"out.filter",quote=FALSE ,sep="\\t")

RSCRIPT



close RCMD;
my $ret = system("R --restore --no-save < cmd.r");
if ( $ret ) {
        die "Error, died with $ret";
}

copy('tier.pdf', $output4);
copy('boxplot.pdf', $output5);
getRTable( 'out.all', $output1 );
getRTable( 'out.filter', $output2 );

get_html_link_v2($html, "Biosigner analysis Results",  "$output1", "tabH","$output2", "tabH",  "$output4", "pdf", "$output5", "pdf",  );





# getRTable(in, out)
sub getRTable{
    my $in = shift;
    my $out = shift;
    my $name = shift;
    open IN,"$in";
    open OUT,">$out";
    my $h = <IN>;
    if( defined($name)){
        print OUT "$name\t$h";
    }else{
        print OUT "\t$h";
    }
    while (<IN>){
          print OUT "$_";
    }
    close IN;
    close OUT;
}


sub getRTable_v2{
    my $in = shift;
    my $out = shift;
    open IN,"$in";
    open OUT,">$out";
    my $h = <IN>;
    $h=~s/holm/p_holm/;
    $h=~s/bonferroni/p_bonferroni/;
    print OUT "\t$h";
    while (<IN>){
          print OUT "$_";
    }
    close IN;
    close OUT;
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