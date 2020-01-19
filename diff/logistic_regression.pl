#!/usr/bin/perl -w
use strict;
use warnings;

use Getopt::Long;
use Scalar::Util qw(looks_like_number);
use File::Copy;
use File::Copy::Recursive qw(dircopy);
use File::Basename;
use FindBin qw($Bin);
my $r_path = "$Bin/../../bin/R/bin";
$ENV{'PATH'} = "$r_path;$ENV{'PATH'}";

my %opts;
GetOptions (\%opts,"i=s", "g=s", "html=s", 
);


my $input = $opts{i};
my $group_file = $opts{g};
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
deal_NA_data( $input, "__noNA_data__" );
$input = "__noNA_data__";


check_data_v3($input, $group_file, "F");
copy($input, '__matrix__');
copy($group_file, '__group__');
open RCMD, ">cmd.r";
print RCMD <<"RSCRIPT";

library(caret)
library(car)
library(tibble)
library(magrittr)
library(dplyr)
library(stringr)
library(readr)
library(ggplot2)
library(ggrepel)
library(tidyr)
library(PRROC)
library(pROC)

data = read.table(quote="", "__matrix__", check.names = FALSE, header = TRUE, row.names = 1, sep = "\t", comment.char = "")
group = read.table(quote="", "__group__", check.names = FALSE, header = F, row.names = 1,  com = '', sep = "\t", colClasses=c("character","character"))

colnames(group) = c("Group")

group\$Group = as.character(group\$Group)
group\$Group = as.factor(group\$Group)

sampleNames = rownames(group)
data = t(data[sampleNames])
data = as.data.frame(data)
group =as.data.frame(group)
data1 = cbind(data, group)


glm = glm(Group ~ ., data1, family=binomial(link="logit"))
varImp = varImp(glm, scale = TRUE)


varImpDf <- varImp \%>\%
    rownames_to_column("Metabolite") \%>\%
    rename(VarImp = Overall) \%>\%
    arrange(desc(VarImp)) \%>\%
    mutate(Metabolite=str_replace_all(Metabolite,"`",""))


write.csv(varImpDf, "LR_VarImp.csv", row.names = F)
write.table(varImpDf, file = "LR_VarImp.txt", quote = FALSE, sep = "\t", row.names = F)

uniq.group <- as.character(unique(data1\$Group))
Yhat = fitted(glm)
thresh = 0.5
YhatFac = cut(Yhat, breaks=c(-Inf, thresh, Inf), labels=c( uniq.group[1] ,  uniq.group[2]))
Yfac = as.factor(data1\$Group)
sum = length(YhatFac)
count = 0

for( i in 1:length(YhatFac)){
    if( as.character(YhatFac[i]) == as.character(Yfac[i]) ){
        count = count + 1 
    }
}
if( count / sum < 0.5 ){
    YhatFac = cut(Yhat, breaks=c(-Inf, thresh, Inf), labels=c( uniq.group[2] ,  uniq.group[1]))
}

pre_summary = table(YhatFac, Yfac)
write.table(pre_summary, file = "pre_sum.txt", quote = FALSE, sep = "\t")

cTab = table(Yfac, YhatFac)
odd_ratio = (cTab[1, 1] / cTab[1, 2]) / (cTab[2, 1] / cTab[2, 2])
write.table(odd_ratio, file = "or.txt", quote = FALSE, sep = "\t", row.names=F, col.names=F)

Probability=mutate_all(as_tibble(Yhat), 
    function(x){
        ifelse(x > 0.5, x, 1 - x)
    })
names(Probability) = "Probability"
samDF <- cbind.data.frame(group, LR_prediction = YhatFac, LR_Fitted_Values = Yhat, Probability)
write.table(samDF, file = "samDF.txt", quote = FALSE, sep = "\t")

uniqGroup <- unique(samDF\$Group)
group1Data <- samDF \%>\%
    filter(Group == uniqGroup[1])
group2Data <- samDF \%>\%
    filter(Group == uniqGroup[2])

x <- samDF\$LR_Fitted_Values
y <- samDF\$Group

otherRocRs <- roc(y, x, ci = T)
rocRs <- roc.curve(group2Data\$LR_Fitted_Values, group1Data\$LR_Fitted_Values, curve = T)
auc <- rocRs\$auc
if (auc < 0.5) {
    rocRs <- roc.curve(group1Data\$LR_Fitted_Values, group2Data\$LR_Fitted_Values, curve = T)
    auc <- rocRs\$auc
}
rocDf <- rocRs\$curve \%>\%
  as.data.frame() \%>\%
  set_colnames(c("FPR", "Sensitivity", "Cutoff")) \%>\%
  select("Cutoff", everything())

write.csv(rocDf, "ROC_Curve_Data.csv", row.names = F)
write.table(rocDf, file = "ROC_Curve_Data.txt", quote = FALSE, sep = "\t", row.names = F)

plotData <- rocRs\$curve \%>\%
  as.data.frame()

colors <- hsv(h = seq(0, 1, length = 100) * 0.8, s = 1, v = 1)
yBreaks <- seq(0, 1, 0.2)

pointData <- coords(otherRocRs, "best", transpose = FALSE) \%>\%
  round(3) \%>\%
  as.data.frame()
ci <- otherRocRs\$ci
ci1 <- ci[1] \%>\%
  round(3)
ci2 <- ci[3] \%>\%
  round(3)

p <- ggplot(data = plotData, aes(x = V1, y = V2, color = V3)) +
  theme_bw(base_size = 8.8) +
  theme(plot.title = element_text(size = 14, hjust = 0.5, face = "bold"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        plot.margin = unit(c(1, 0.5, 1, 0.5), "cm"), panel.border = element_rect(size = 0.75)
  ) +
  geom_line() +
  geom_abline(intercept = 0, slope = 1, color = "grey", linetype = 1, size = 0.4) +
  annotate("text", x = 0.5 + 0.01, y = 0.5, label = str_c("AUC:", format(auc, digits = 3), "(", ci1, ", ", ci2, ")"), color = "black",
           hjust = 0, size = 3.5, family = "Times") +
  geom_point(data = pointData, aes(x = 1 - specificity, y = sensitivity), size = 1, color = "red") +
  xlab("FPR") +
  geom_text(data = pointData, aes(x = 1 - specificity, y = sensitivity, label = paste0(threshold,
                                                                                       "(", specificity, ", ", sensitivity, ")")), color = "black", hjust = 0, vjust = 1, size = 3, nudge_x = 0.015,
            nudge_y = -0.015, family = "Times") +
  ylab("Sensitivity") +
  ggtitle(str_c("AUC=", format(auc, digits = 3))) +
  labs(colour = "Cutoff") +
  scale_x_continuous("1 - Specificity", breaks = yBreaks) +
  scale_y_continuous("Sensitivity", breaks = seq(0, 1, 0.2)) +
  scale_colour_gradientn(colours = colors, breaks = yBreaks)
ggsave("ROC_Curve.pdf", p, width = 6, height = 6)


prRs <- pr.curve(group2Data\$LR_Fitted_Values, group1Data\$LR_Fitted_Values, curve = T)
auc <- prRs\$auc.integral
if (auc < 0.5) {
    prRs <- pr.curve(group1Data\$LR_Fitted_Values, group2Data\$LR_Fitted_Values, curve = T)
    auc <- prRs\$auc.integral
}

prDf <- prRs\$curve \%>\%
    as.data.frame() \%>\%
    set_colnames(c("Recall", "Precision", "Cutoff")) \%>\%
    select("Cutoff", everything())
write.csv(prDf, "PR_Curve_Data.csv", row.names = F)
write.table(prDf, file = "PR_Curve_Data.txt", quote = FALSE, sep = "\t", row.names = F)

colors <- hsv(h = seq(0, 1, length = 100) * 0.8, s = 1, v = 1)
yBreaks <- seq(0, 1, 0.2)

plotData <- prRs\$curve \%>\%
    as.data.frame() \%>\%
    arrange(desc(V2), V1)

p <- ggplot(data = plotData, aes(x = V1, y = V2, color = V3)) +
    theme_bw(base_size = 8.8) +
    theme(plot.title = element_text(size = 14, hjust = 0.5, face = "bold"), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
    plot.margin = unit(c(1, 0.5, 1, 0.5), "cm"), panel.border = element_rect(size = 0.75)
    ) +
    geom_line() +
    xlab("FPR") +
    ylab("Sensitivity") +
    ggtitle(str_c("AUC=", format(auc, digits = 3))) +
    labs(colour = "Cutoff") +
    scale_x_continuous("Recall", breaks = yBreaks) +
    scale_y_continuous("Precision", breaks = seq(0, 1, 0.2), limit = c(0, 1)) +
    scale_colour_gradientn(colours = colors, breaks = yBreaks)

ggsave("PR_Curve.pdf", p, width = 6, height = 6)

RSCRIPT
close RCMD;
my $ret = system("R --restore --no-save < cmd.r");
if ( $ret ) {
        die "Error, died with $ret";
}


my $html_dir = dirname($html);
getRTable('pre_sum.txt', "$html_dir/LR_Prediction_Summary.txt",  "Prediction");
getRTable('samDF.txt', "$html_dir/LR_Prediction.txt" );
copy("LR_VarImp.txt", "$html_dir/LR_VarImp.txt" );
copy("PR_Curve_Data.txt", "$html_dir/PR_Curve_Data.txt" );
copy("ROC_Curve_Data.txt", "$html_dir/ROC_Curve_Data.txt" );

&get_html_link_4lr($html, "Logistic Regression Result", 
    "$html_dir/LR_VarImp.txt", "tabH", 
    "$html_dir/LR_Prediction.txt", "tabH", 
    "$html_dir/LR_Prediction_Summary.txt", "or_table", 
    
    "$html_dir/ROC_Curve_Data.txt", "tabH", 
    "ROC_Curve.pdf", "pdf", 
    
    "$html_dir/PR_Curve_Data.txt", "tabH",  
    "PR_Curve.pdf", "pdf", 
);


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
            elsif(  $type eq "or_table" ){
                get_html_table_or($file, "$out_dir/$name.html", "y");
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
            }elsif(  $type eq "txt_v2" ){
                print HTML "<li><a href=\"$basename/$name\" target=\"_parent\">$name</a></li>\n";
            }elsif(  $type eq "relation" ){
                print HTML "<li><a href=\"$basename/$name.html\" target=\"_parent\">$name</a></li>\n";
            }
            elsif(  $type eq "or_table" ){
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


sub get_html_table_or{
	my $in = shift;
	my $out = shift;
	my $isHeader = shift;
    open IN,"or.txt";
    my $or = <IN>;
    $or = trim($or);
    close IN;


	open IN,"$in";
	open OUT,">$out";
    print OUT "<!DOCTYPE html>\n";
    print OUT "<html>\n";
	print OUT "<table border=\"1\">\n";
    
    print OUT "<tr align=\"center\"><th colspan=\"3\">Odd Ratio: $or</th></tr>";
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



sub trim
{
        my $string = shift;
        $string =~ s/^\s+//;
        $string =~ s/\s+$//;
        return $string;
}



sub get_html_link_4lr{
	my @param = @_;
	my $html = shift @param;
	my $title = shift @param;
	my @par = @param;

    my $out_dir = "$html.files";
    mkdir $out_dir;
    my $html_dir = dirname($html);

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
            }elsif(  $type eq "pdf"){
                copy($file, "$html_dir/$name");
            }elsif(  $type eq "dir" ){
                dircopy(  $file , "$html_dir/$name") or die $!;
            }elsif(  $type eq "or_table" ){
                get_html_table_or($file, "$out_dir/$name.html", "y");
            }
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
            }elsif( $type eq "pdf" ){
                print HTML "<li><a href=\"$name\" target=\"_parent\">$name</a></li>\n";
            }elsif(  $type eq "dir" ){
                print HTML "<li><a href=\"$name\" target=\"_parent\">$name</a></li>\n";
            }elsif(  $type eq "or_table" ){
                print HTML "<li><a href=\"$basename/$name.html\" target=\"_parent\">$name</a></li>\n";
            }
            else{
                print HTML "<li><a href=\"$basename/$name\" target=\"_parent\">$name</a></li>\n";
            }
	        
	}
	print HTML "</ul></p>\n";


}