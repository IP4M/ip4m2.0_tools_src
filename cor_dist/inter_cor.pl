#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Scalar::Util qw(looks_like_number);
use File::Copy;
use File::Copy::Recursive qw(dircopy);
use File::Basename;
use FindBin qw($Bin);
use FindBin qw($Bin);
my $r_path = "$Bin/../../bin/R/bin";
$ENV{'PATH'} = "$r_path;$ENV{'PATH'}";

my %opts;
GetOptions (\%opts, "i1=s",  "method=s", "html=s", "corrected=s",  "scale=s", "i2=s",
"pcut=s", "mtcuf=s", "rcuf=s",
);

my $input1 = $opts{i1};
my $input2 = $opts{i2};
my $zscal = $opts{scale};
my $html = $opts{html};
my $method = $opts{method};
my $corrected = $opts{corrected};
my $pcut = $opts{pcut};
my $mtcuf = $opts{mtcuf};
my $rcuf = $opts{rcuf};

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
deal_NA_data( $input1, "__noNA_data1__" );
$input1 = "__noNA_data1__";

deal_NA_data( $input2, "__noNA_data2__" );
$input2 = "__noNA_data2__";


&check_data_v4( $input1  );
&check_data_v4( $input2  );
&check_data_sets( $input1, $input2 );

copy($input1, '__matrix1__');
copy($input2, '__matrix2__');

open RCMD, ">cmd.r";
print RCMD <<"RSCRIPT";

corTest = function(da1, da2, method){
    name1 <- colnames(da1)
    name2 <- colnames(da2)
    tr_da1 <- as.matrix(da1)
    tr_da2 <- as.matrix(da2)
    num1 <- ncol(da1)
    num2 <- ncol(da2)
    pvalue <- vector()
    value <- vector()
    id1 <- vector()
    id2 <- vector()
    rec <- 1;
    for (i in 1 : num1) {
        for (j in 1 : num2) {
            a1 <- tr_da1[, i]
            a2 <- tr_da2[, j]
            id1[rec] <- name1[i]
            id2[rec] <- name2[j]
            a1 <- as.numeric(a1)
            a2 <- as.numeric(a2)
            
            tryCatch(corr <- cor.test(a1, a2, method = method), error = function(e) {
                value[rec] <- 0
                pvalue[rec] <- 1              
            },finally = {
                esti <- as.matrix(corr\$estimate)[1]
                value[rec] <- esti
                pvalue[rec] <- corr\$p.value
            })

            rec <- rec + 1
        }
    }
    list(node1 = id1, node2 = id2, cor = value, p = pvalue)
}


library(tibble)
library(magrittr)
library(dplyr)
library(stringr)
library(readr)
library(tidyr)
library(purrr)
library(circlize)
zscal = $zscal

data <- read.table("__matrix1__", header = T, sep = "\t", check.names = FALSE,quote="", row.names = 1)
data = t(data)
if( zscal){
    data = scale(data)
}

extraData <- read.table("__matrix2__", header = T, sep = "\t", check.names = FALSE,quote="", row.names = 1)
extraData = t(extraData)
if( zscal){
    extraData = scale(extraData)
}

dataSampleIds = rownames(data)
extraSampleIds = rownames(extraData)
interSampleIds <- intersect(dataSampleIds, extraSampleIds)

data = data[interSampleIds,]
extraData = extraData[interSampleIds,]

listRs <- corTest(data, extraData, "$method")

allData <- data.frame(Node1 = listRs\$node1, Node2 = listRs\$node2, r = listRs\$cor, P = listRs\$p,
    stringsAsFactors = F) \%>\%
        mutate_at(vars(c("r")), function(x){
            ifelse(is.na(x), 0, x)
        }) \%>\%
        mutate_at(vars(c("P")), function(x){
            ifelse(is.na(x), 1, x)
        }) \%>\%
        mutate(FDR = p.adjust(P, method = "$corrected"))

corData <- allData \%>\%
    select(c("Node1", "Node2", "r")) \%>\%
    spread(Node1, "r") \%>\%
    rename(` ` = Node2)
    
write_csv(corData,  "Inter_r_Matrix.csv")
write.table(corData, file = "Inter_r_Matrix.txt", quote = FALSE, sep = "\t", row.names=F)

pData <- allData \%>\%
    select(c("Node1", "Node2", "P")) \%>\%
    spread(Node1, "P") \%>\%
    rename(` ` = Node2)

write_csv(pData, "Inter_P_Matrix.csv")
write.table(pData, file = "Inter_P_Matrix.txt", quote = FALSE, sep = "\t", row.names=F)

fdrData <- allData %>%
    select(c("Node1", "Node2", "FDR")) \%>\%
    spread(Node1, "FDR") \%>\%
    rename(` ` = Node2)

write_csv(fdrData,  "Inter_Corrected_P_Matrix.csv")
write.table(fdrData, file = "Inter_Corrected_P_Matrix.txt", quote = FALSE, sep = "\t", row.names=F)

edgeData <- allData \%>\%
        filter(P < $pcut & FDR < $mtcuf & abs(r) > $rcuf) \%>\%
        filter(Node1 != Node2) \%>\%
        mutate(distName = {
            Node1 \%>\%
            map2_chr(Node2, function(x, y){
                vec <- c(x, y) \%>\%
                sort()
                str_c(vec, collapse = ";")
            })
        }) \%>\%
        distinct(distName, .keep_all = T) \%>\%
        select(- c("distName")) \%>\%
        arrange(desc(abs(r)))

write_csv(edgeData, "Network_Edges_for_Cytoscape.csv")
write.table(edgeData, file = "Network_Edges_for_Cytoscape.txt", quote = FALSE, sep = "\t", row.names=F)

nodes <- unique(c(edgeData\$Node1, edgeData\$Node2))
infoData <- tibble(Node = nodes) \%>\%
        mutate(Size = {
            Node \%>\%
            map_int(function(x){
                edgeData \%>\%
                    filter(Node1 == x | Node2 == x) \%>\%
                    nrow()
            })
        }) 

write_csv(infoData, "Network_Nodes_for_Cytoscape.csv")
write.table(infoData, file = "Network_Nodes_for_Cytoscape.txt", quote = FALSE, sep = "\t", row.names=F)

library(ComplexHeatmap)



corData <- read_csv(  "Inter_r_Matrix.csv" )
plotData <- corData \%>\%
    column_to_rownames("X1") \%>\%
    as.matrix()
pData <- read_csv( "Inter_P_Matrix.csv" ) \%>\%
    column_to_rownames("X1") \%>\%
    as.matrix()

colors <- rev(colorRampPalette(RColorBrewer::brewer.pal(10, "RdBu"))(256))

rowNum <- nrow(plotData)
colNum <- ncol(plotData)
eachHeight <- 1 / rowNum

row_names = rownames(plotData)
col_names = colnames(plotData)
max_row_chr_num = max(nchar(row_names))
max_col_chr_num = max(nchar(col_names))

height <- max(2 + (rowNum - 10) * 0.1, 2) + 1 / 4 * max_col_chr_num
height
width <- max(2 + (colNum - 10) * 0.1, 2) + 1 / 4 * max_row_chr_num
width

fileName <- "Inter_Correlation_Heatmap.pdf"
pdf(fileName, width = width, height = height)

  htList <- Heatmap(plotData, col = colors, show_column_names = T, cluster_rows = T, cluster_columns = T,
                    #name = "", 
                    row_names_gp = gpar(fontsize = 6, fontfamily = "Times"), show_row_names = T, show_heatmap_legend = F,
                    column_names_gp = gpar(fontsize = 6, fontfamily = "Times"), cell_fun = function(j, i, x, y, width, height, fill) {
      value <- pData[i, j]
      str <- if (value < 0.01) {
        "+"
      }else if (value < 0.05) {
        "*"
      }else ""
      top <- if (str == "*") {
        y - unit(eachHeight * 0.12, "npc")
      }else y
      grid.text(str, x, top, gp = gpar(fontsize = 9), hjust = 0.5, vjust = 0.5, just = "center")
    })

  rLgd <- Legend(col_fun = colorRamp2(seq(-1, 1, length.out = 256), colors), at = seq(-1, 1, 0.5), title = "",
                 title_gp = gpar(fontsize = 12, fontfamily = "Times"), grid_height = unit(7.5, "mm"),
                 labels_gp = gpar(fontsize = 12, fontfamily = "Times"))

  pLegend <- Legend(pch = c("+", "*"), type = "points", labels = c("p<0.01", "0.01<=p<0.05"), title = "",
                    labels_gp = gpar(fontsize = 9, fontfamily = "Times"), title_gp = gpar(fontsize = 12, fontfamily = "Times"),
                    ncol = 1, grid_height = unit(5, "mm"),
  )
  pd <- packLegend(list = list(rLgd, pLegend), row_gap = unit(0.5, "cm"))

  draw(htList, heatmap_legend_side = "right", annotation_legend_side = "right",
       padding = unit(c(0.5, 0.5, 0.5, 0.5), "cm"), annotation_legend_list = pd
  )
  dev.off()







RSCRIPT
close RCMD;
my $ret = system("R --restore --no-save < cmd.r");
if ( $ret ) {
        die "Error, died with $ret";
}




my $html_dir = dirname($html);
copy("Inter_r_Matrix.txt", "$html_dir/Inter_r_Matrix.txt");
copy("Inter_P_Matrix.txt", "$html_dir/Inter_P_Matrix.txt");
copy("Inter_Corrected_P_Matrix.txt", "$html_dir/Inter_Corrected_P_Matrix.txt");
copy("Network_Edges_for_Cytoscape.txt", "$html_dir/Network_Edges_for_Cytoscape.txt");
copy("Network_Nodes_for_Cytoscape.txt", "$html_dir/Network_Nodes_for_Cytoscape.txt");

&get_html_link_v3($html, "Inter Features Correlation Analysis Results", 
    "$html_dir/Inter_r_Matrix.txt", "relation", 
    "Inter_Correlation_Heatmap.pdf", "pdf", 
    "$html_dir/Inter_P_Matrix.txt", "relation", 
    "$html_dir/Inter_Corrected_P_Matrix.txt", "relation", 
    "$html_dir/Network_Edges_for_Cytoscape.txt", "tabH",     
    "$html_dir/Network_Nodes_for_Cytoscape.txt", "tabH", 


);




sub check_data_sets {
    my $in1 = shift;
    my $in2 = shift;
    my $f1 = basename $in1;
    my $f2 = basename $in2;
    open IN1, "$in1";
    open IN2, "$in2";
    my $h1 = <IN1>;
    my $h2 = <IN2>;
    my @h1 = split /\t/, $h1;
    my @h2 = split /\t/, $h2;
    shift @h1;
    shift @h2;
    my %h1;
    foreach(@h1){
        $h1{$_} = 1;
    }
    my $cnt = 0;
    foreach(@h2){
        if(defined($h1{$_})){
            $cnt ++;
        }
    }
    if ( $cnt < 2) {
        die "Error: The intersection samples of file \"$f1\" and \"$f2\" < 2 \n";
    }

    my %feature1;
    while(<IN1>){
        chomp;
        my @s = split /\t/, $_;
        $feature1{$s[0]} = 1;
    }
    while(<IN2>){
        chomp;
        my @s = split /\t/, $_;
        if ( defined($feature1{$s[0]}) ) {
            die "Error: The feature \"$s[0]\" is replicated in \"$f1\" and \"$f2\" \n";
        }
    }

    close IN1;
    close IN2;
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