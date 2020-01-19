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
GetOptions (\%opts, "i=s", "g=s", "vip=s",  "html=s",  "zscal=s", 
);

my $matrix = $opts{i};
my $zscal = $opts{zscal};
my $html = $opts{html};
my $group_file = $opts{g};
my $vip = $opts{vip};


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
deal_NA_data( $matrix, "__noNA_data__" );
$matrix = "__noNA_data__";


check_data_v3($matrix, $group_file, "F");

copy($matrix, '__matrix__');
copy($group_file, '__group__');
get_csv_group_file($group_file, '__group2__');

open RCMD, ">cmd.r";
print RCMD <<"RSCRIPT";

library(ropls)
library(tibble)
library(magrittr)
library(dplyr)
library(stringr)
library(readr)
library(ggplot2)
library(ggrepel)
library(tidyr)
zscal = $zscal
xMN <- read.table(quote="", "__matrix__", check.names = FALSE, header = TRUE, row.names = 1, sep = "\t", comment.char = "")
tmpDf = xMN
samDF <- read.table(quote="", "__group__", check.names = FALSE, header = F, row.names = 1,  com = '', sep = "\t")
colnames(samDF) = c("Group")

samDF\$Group = as.character(samDF\$Group)
samDF\$Group = as.factor(samDF\$Group)

sampleNames = rownames(samDF)
sampleSize = length(sampleNames)
xMN = xMN[rownames(samDF)]
xMN = t(xMN)
yMCN <- matrix(samDF[, "Group"], ncol = 1, dimnames = list(rownames(xMN), "Group"))

if( zscal){
    xMN = scale(xMN)
}
crossValI <- min(nrow(xMN), 7)
plsdaRs <- opls(x = xMN, y = yMCN, orthoI = 1, plotL = F, crossvalI = crossValI, permI = 999)
modelDf <- plsdaRs\@modelDF \%>\%
  rownames_to_column("poName") \%>\%
  mutate(poName = poName \%>\%
    str_replace("p1", "P1") \%>\%
    str_replace("o1", "O1")) \%>\%
  column_to_rownames("poName")

scoreMN <- plsdaRs\@scoreMN \%>\%
  as.data.frame() \%>\%
  rownames_to_column("SampleID")
orthoScoreMN <- plsdaRs\@orthoScoreMN \%>\%
  as.data.frame() \%>\%
  rownames_to_column("SampleID")
plotData <- scoreMN \%>\%
  inner_join(orthoScoreMN, by = c("SampleID")) \%>\%
  rename_at(vars(-c("SampleID")), function(x) {
    x \%>\%
      str_replace("p1", "P1") \%>\%
      str_replace("o1", "O1")
  }) \%>\%
  column_to_rownames("SampleID")

write.csv(plotData, "OPLSDA_Score.csv")
write.table(plotData, "OPLSDA_Score.txt", quote = FALSE, sep = "\t")
write.csv(modelDf, "OPLSDA_R2X_R2Y_Q2.csv", quote = FALSE)
write.table(modelDf, "OPLSDA_R2X_R2Y_Q2.txt", quote = FALSE, sep = "\t")


method <- "pearson"
pData <- plsdaRs\@suppLs[["xModelMN"]] \%>\%
    as.data.frame() \%>\%
    summarise_all(function(x){cor.test(x, plsdaRs\@scoreMN, method = method)\$p.value}) \%>\%
    unlist()
vipVn = getVipVn(plsdaRs)
corData <- cor(plsdaRs\@suppLs[["xModelMN"]], plsdaRs\@scoreMN, method = method) \%>\%
    as.data.frame() \%>\%
    .\$p1
vipData <- data.frame(VIP = vipVn, Corr.Coeffs. = corData, Corr.P = pData) %>%
    rownames_to_column("Metabolite")
allData <- vipData \%>\%
     mutate(FDR = p.adjust(Corr.P, method = "fdr")) \%>\%
     arrange(desc(VIP))
write.csv(allData, "OPLSDA_VIP.csv", row.names = FALSE)
write.table(allData, "OPLSDA_VIP.txt", quote = FALSE, sep = "\t" ,row.names = FALSE )

sigallData =  subset(allData, VIP >= $vip)
write.table(sigallData, "OPLSDA_VIP_Sig.txt", quote = FALSE, sep = "\t" ,row.names = FALSE )

perMN <- plsdaRs\@suppLs\$permMN \%>\%
    as.data.frame() \%>\%
    set_colnames(c("R2X", "R2Y", "Q2", "RMSEE", "pre", "ort", "Similarity")) \%>\%
    select(- c("pre", "ort")) \%>\%
    arrange(Similarity) \%>\%
    filter(Similarity > 0)

write.csv(perMN, "OPLSDA_Permutation.csv")
write.table(perMN, "OPLSDA_Permutation.txt", quote = FALSE, sep = "\t")

fixData <- perMN \%>\%
    filter(Similarity == 1) \%>\%
    head(1)

offsetPlotData <- perMN \%>\%
    mutate(Similarity = Similarity - 1, Q2 = Q2 - fixData\$Q2, R2Y = R2Y - fixData\$R2Y)

lmRs <- lm(Q2 ~ Similarity + 0, offsetPlotData)
q2A <- coef(lmRs)["Similarity"]
q2B <- fixData\$Q2 - q2A * 1

lmRs <- lm(R2Y ~ Similarity + 0, offsetPlotData)
r2A <- coef(lmRs)["Similarity"]
r2B <- fixData\$R2Y - r2A * 1

tbf <- tibble(` ` = c("Q2", "R2"), `vertical intercept` = c(q2B, r2B), slope = c(q2A, r2A))
write_csv(tbf, "Fitted_Curve_Parameter.csv")
write.table(tbf, "Fitted_Curve_Parameter.txt", quote = FALSE, sep = "\t", row.names = FALSE )

sampleInfo <- read.csv("__group2__", header = T, stringsAsFactors = F) \%>\%
    select(c("SampleID", "ClassNote")) \%>\%
    mutate(ClassNote = as.character(ClassNote))
plotData <- read_csv("OPLSDA_Score.csv") \%>\%
    rename(SampleID = X1) \%>\%
    inner_join(sampleInfo, by = c("SampleID")) \%>\%
    mutate(ClassNote = factor(ClassNote, levels = unique(ClassNote)))
modelDf <- read.csv("OPLSDA_R2X_R2Y_Q2.csv", header = T, stringsAsFactors = F, comment.char = "", row.names = 1)
p <- ggplot(plotData, mapping = aes(x = P1, y = O1, color = ClassNote, label = SampleID, fill = ClassNote)) +
  xlab(paste0("P1 (", modelDf["P1", "R2X"] * 100, "%)")) +
  ylab(paste0("O1 (", modelDf["O1", "R2X"] * 100, "%)")) +
  theme_bw(base_size = 8.8, base_family = "Times") +
  theme(axis.text.x = element_text(size = 9, vjust = 0.5),
        axis.text.y = element_text(size = 8.8), legend.position = 'right',
        axis.title.y = element_text(size = 12), legend.margin = margin(t = 0.3, b = 0.1, unit = 'cm'),
        legend.text = element_text(size = 9), axis.title.x = element_text(size = 12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  ) +
  #0 line
  geom_vline(aes(xintercept = 0), colour = "#BEBEBE", linetype = "solid") +
  geom_hline(aes(yintercept = 0), colour = "#BEBEBE", linetype = "solid") +
  #point
  geom_point(aes(colour = factor(ClassNote)), size = 3, stroke = 0) +
  stat_ellipse(aes(fill = ClassNote), colour = NA, size = 0.3, level = 0.95, type = "norm",
               geom = "polygon", alpha = 0.2, show.legend = F) +
  geom_text_repel(segment.size = 0.2, size = 2, family = "Times") 

ggsave("OPLSDA_Score_2D_Label.pdf", p, width = 5, height = 4)

p <- ggplot(plotData, mapping = aes(x = P1, y = O1, color = ClassNote)) +
  xlab(paste0("P1 (", modelDf["P1", "R2X"] * 100, "%)")) +
  ylab(paste0("O1 (", modelDf["O1", "R2X"] * 100, "%)")) +
  theme_bw(base_size = 8.8, base_family = "Times") +
  theme(axis.text.x = element_text(size = 9, vjust = 0.5),
        axis.text.y = element_text(size = 8.8), legend.position = 'right',
        axis.title.y = element_text(size = 12), legend.margin = margin(t = 0.3, b = 0.1, unit = 'cm'),
        legend.text = element_text(size = 9), axis.title.x = element_text(size = 12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  ) +
  #0 line
  geom_vline(aes(xintercept = 0), colour = "#BEBEBE", linetype = "solid") +
  geom_hline(aes(yintercept = 0), colour = "#BEBEBE", linetype = "solid") +
  #point
  geom_point(aes(colour = factor(ClassNote)), size = 4, stroke = 0) +
  stat_ellipse(aes(fill = ClassNote), colour = NA, size = 0.3, level = 0.95, type = "norm",
               geom = "polygon", alpha = 0.2, show.legend = F)

ggsave("OPLSDA_Score_2D.pdf", p, width = 5, height = 4)

p.cor <- function(n, p, method="pearson"){
    df <- n - 2
    t <- qt(p / 2, df, lower.tail = FALSE)
    r <- switch(method,
    "pearson" = sqrt(t ^ 2 / (df + t ^ 2))
    )
    return(r)
}

vipNum = $vip

method <- "pearson"
n <- nrow(sampleInfo)
r5 <- p.cor(n, 0.05, method)
r5s <- c(- r5, r5)

plotData <- read.csv("OPLSDA_VIP.csv", header = T) \%>\%
    rowwise() \%>\%
    do({
        result <- as.data.frame(.)
        p <- result["Corr.P"]
        vip <- result["VIP"]
        col <- if (p <= 0.05 && vip >= vipNum) {
            "height"
        }else if (p > 0.05 && vip >= vipNum) {
            "median"
        }else "low"
        result\$col <- col
        result
    }) \%>\%
    mutate(col = factor(col, levels = c("height", "median", "low"))) \%>\%
    rename(vip = VIP, cor = Corr.Coeffs.) \%>\%
    as.data.frame()
leftPlotData <- plotData \%>\%
    filter(cor < 0) \%>\%
    filter(col != "low")

rightPlotData <- plotData \%>\%
    filter(cor > 0) \%>\%
    filter(col != "low")

yMax <- max(plotData\$vip)
r5Col <- "#006400"
breaks <- seq(- 1, 1, 0.2)
fontSize <- 1.5

pointCols <- tibble(sig = c("height", "median", "low"), col = c("#006400", "#0000FF", "grey")) \%>\%
  deframe()

p <- ggplot(plotData, mapping = aes(x = cor, y = vip)) +
    theme_bw(base_size = 8.8, base_family = "Times") +
    theme(axis.text.x = element_text(size = 9, vjust = 0.5),
    axis.text.y = element_text(size = 8.8), legend.position = 'none',
    axis.title.y = element_text(size = 11), legend.margin = margin(t = 0.3, b = 0.1, unit = 'cm'),
    legend.text = element_text(size = 6), axis.title.x = element_text(size = 11), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5, size = 13)
    ) +
    ylab("VIP") +
    geom_hline(aes(yintercept = vipNum), colour = "#A9A9A9", linetype = 2, size = 0.375) +
    geom_vline(xintercept = r5s, colour = r5Col, linetype = 2, size = 0.4) +
    annotate("text", x = r5 + 0.01 , y = yMax, label = "P=0.05", color = r5Col, hjust = 1 ,
    size = 2, family = "Times", angle = 90, vjust = 1,) +
    annotate("text", x = - r5 - 0.01 , y = yMax, label = "P=0.05", color = r5Col, hjust = 1 ,
    size = 2, family = "Times", angle = 90, vjust = 0) +
    scale_colour_manual("", values = pointCols) +
    scale_x_continuous("Corr.Coeffs.", breaks = breaks, limits = c(- 1.2, 1.2)) +
    geom_text(aes(x = cor, y = vip, colour = col), label = "+", size = 4, family = "Times",
    hjust = 0.5, vjust = 0.3, fontface = "bold") +
    geom_text(data = leftPlotData, aes(x = cor, y = vip, label = Metabolite, color = col), size = fontSize, family = "Times",
    nudge_x = - 0.015, hjust = 1) +
# geom_text_repel(segment.size=0.2,data = rightPlotData, aes(x = cor, y = vip, label = Metabolite,color=col), size = 2, family = "Times")
    geom_text(data = rightPlotData, aes(x = cor, y = vip, label = Metabolite, color = col), size = fontSize, family = "Times",
    hjust = 0, nudge_x = 0.015,)

ggsave("OPLSDA_VPlot.pdf", p, width = 6, height = 4)

modelDf <- read.csv("OPLSDA_R2X_R2Y_Q2.csv", header = T, stringsAsFactors = F, comment.char = "", row.names = 1)
plotData <- modelDf \%>\%
    rownames_to_column("poName") \%>\%
    filter(poName \%in\% c("P1", "O1")) \%>\%
    select(c("poName", "R2X", "R2Y", "Q2")) \%>\%
    gather("group", "value", - poName) \%>\%
    mutate(group = factor(group, levels = c("R2X", "R2Y", "Q2"))) \%>\%
    mutate(poName = factor(poName, levels = c("P1", "O1")))

p <- ggplot(plotData, mapping = aes(x = poName, y = value, fill = group, label = value)) +
    xlab("") +
    ylab("") +
    theme_bw(base_size = 8.8, base_family = "Times") +
    theme(axis.text.x = element_text(size = 9, vjust = 0.5),
    axis.text.y = element_text(size = 8.8), legend.position = 'right',
    axis.title.y = element_text(size = 12), legend.margin = margin(t = 0.3, b = 0.1, unit = 'cm'),
    legend.text = element_text(size = 9), axis.title.x = element_text(size = 12),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    ) +
    geom_bar(stat = "identity", position = position_dodge()) +
    geom_text(size = 3.5, family = "Times", vjust = - 0.3, position = position_dodge(0.9)) + 
    scale_colour_manual("", values =  c("#ADD8E6", "#FFE4E1", "#E6E6FA")) +
    scale_fill_manual("", values = c("#ADD8E6", "#FFE4E1", "#E6E6FA"))
ggsave("OPLSDA_R2X_R2Y_Q2.pdf", p, width = 7, height = 5)

plotData <- read.csv("OPLSDA_Permutation.csv", header = T)
fixData <- plotData \%>\%
    filter(Similarity == 1) \%>\%
    head(1)

parameterData <- read_csv("Fitted_Curve_Parameter.csv")
q2Row <- parameterData \%>\%
    filter(X1 == "Q2")
q2A <- q2Row[1, "slope"] \%>\%
    as.numeric()
q2B <- q2Row[1, "vertical intercept"] \%>\%
    as.numeric()

r2Row <- parameterData \%>\%
    filter(X1 == "R2")
r2A <- r2Row[1, "slope"] \%>\%
    as.numeric()
r2B <- r2Row[1, "vertical intercept"] \%>\%
    as.numeric()

lineX <- plotData\$Similarity
lineY <- q2A * lineX + q2B
plotData\$y <- lineY

lineY <- r2A * lineX + r2B
plotData\$y1 <- lineY

plotData1 <- plotData \%>\%
    select(c("X", "Similarity", "Q2")) \%>\%
    rename(value = Q2) \%>\%
    mutate(class = "Q2", fill = "#0000FF", pch = 22)

plotData2 <- plotData \%>\%
    select(c("X", "Similarity", "R2Y")) \%>\%
    rename(value = R2Y) \%>\%
    mutate(class = "R2", fill = "#4CAF50", pch = 21)
plotDf <- rbind(plotData1, plotData2) \%>\%
    mutate(pch = as.numeric(pch)) \%>\%
    mutate(class = factor(class, levels = c("R2", "Q2"))) \%>\%
    as.data.frame()

colDf <- plotDf \%>\%
    select(c("class", "fill")) \%>\%
    deframe()
pchDf <- plotDf %>%
    select(c("class", "pch")) \%>\%
    deframe()

labels <- c(substitute(paste(R ^ 2, "Y=", R2Y), list(R2Y = fixData\$R2Y)),
substitute(paste(Q ^ 2, "Y=", Q2), list(Q2 = fixData\$Q2)))

p <- ggplot(plotDf, mapping = aes(x = Similarity, y = value)) +
    xlab("Similarity") +
    ylab("") +
    ylim(c(- 1, 1)) +
    xlim(c(0, 1.1)) +
    theme_bw(base_size = 8.8, base_family = "Times") +
    theme(axis.text.x = element_text(size = 9, vjust = 0.5),
    axis.text.y = element_text(size = 8.8), legend.position = 'right',
    axis.title.y = element_text(size = 12), legend.margin = margin(t = 0.3, b = 0.1, unit = 'cm'),
    legend.text = element_text(size = 9, margin = margin(l = 0, unit = "cm"), hjust = 0), axis.title.x = element_text(size = 12), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), legend.spacing.x = unit(0, 'cm')
    ) +
#0 line
    geom_vline(aes(xintercept = 0), colour = "#BEBEBE", linetype = "solid") +
    geom_hline(aes(yintercept = 0), colour = "#BEBEBE", linetype = "solid") +
#point
    geom_line(data = plotData, aes(x = Similarity, y = y, group = 1), linetype = 2, size = 1.2) +
    geom_line(data = plotData, aes(x = Similarity, y = y1, group = 1), linetype = 2, size = 1.2) +
    geom_point(aes(fill = class, pch = class), size = 3, color = "#000000") +
    scale_fill_manual("", values = colDf, labels = labels) +
    scale_shape_manual("", values = pchDf, labels = labels)

ggsave("OPLSDA_Permutation.pdf", p, width = 6, height = 5)


RSCRIPT
close RCMD;
my $ret = system("R --restore --no-save < cmd.r");
if ( $ret ) {
        die "Error, died with $ret";
}

my $html_dir = dirname($html);
getRTable_oplsda("OPLSDA_Score.txt", "$html_dir/OPLSDA_Score.txt");
copy("OPLSDA_VIP.txt", "$html_dir/OPLSDA_VIP.txt");
getRTable_oplsda("OPLSDA_Permutation.txt", "$html_dir/OPLSDA_Permutation.txt");
copy("Fitted_Curve_Parameter.txt", "$html_dir/Fitted_Curve_Parameter.txt");
getRTable_oplsda("OPLSDA_R2X_R2Y_Q2.txt", "$html_dir/OPLSDA_R2X_R2Y_Q2.txt");
copy("OPLSDA_VIP_Sig.txt", "$html_dir/OPLSDA_VIP_Sig.txt");

&get_html_link_v3($html, "OPLS-DA Analysis Results", 
    "$html_dir/OPLSDA_Score.txt", "tabH", 
    "$html_dir/OPLSDA_VIP.txt", "tabH", 
    "$html_dir/OPLSDA_VIP_Sig.txt", "tabH", 
    "$html_dir/OPLSDA_Permutation.txt", "tabH", 
    "$html_dir/Fitted_Curve_Parameter.txt", "tabH", 
    "$html_dir/OPLSDA_R2X_R2Y_Q2.txt", "tabH", 
    "OPLSDA_VPlot.pdf", "pdf", 
    "OPLSDA_Score_2D_Label.pdf", "pdf", 
    "OPLSDA_Score_2D.pdf", "pdf", 
    "OPLSDA_Permutation.pdf", "pdf", 
    "OPLSDA_R2X_R2Y_Q2.pdf", "pdf", 
);




sub getRTable_oplsda{
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



# getRTable(in, out)
sub getRTable_v10{
    my $in = shift;
    my $out = shift;

    open IN,"$in";
    open OUT,">$out";
    my $h = <IN>;
    chomp $h;
    my @s = split /\t/, $h;
    pop @s;
    $h = join "\t", @s;
    print OUT "\t$h\n";
    while (<IN>){
          chomp;
          my @s = split /\t/, $_;
          pop @s;
          my $line = join "\t", @s;
          print OUT "$line\n";
    }
    close IN;
    close OUT;
}


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



sub get_csv_group_file{
    my $in = shift;
    my $out = shift;
    open IN,"$in";
    open OUT,">$out";
    print OUT "SampleID,ClassNote\n";
    while (<IN>){
          chomp;
          my @s = split /\t/, $_;
          my $l = join ",", @s;
          print OUT "$l\n";
    }
    close IN;
    close OUT;
}



sub get_html_link_v3{
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
            }
            else{
                print HTML "<li><a href=\"$basename/$name\" target=\"_parent\">$name</a></li>\n";
            }
	        
	}
	print HTML "</ul></p>\n";


}