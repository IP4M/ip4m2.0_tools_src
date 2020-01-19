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
GetOptions (\%opts, "i1=s",  "method=s", "html=s", "corrected=s",  "scale=s", "i2=s", "i3=s",
"pcut=s", "mtcuf=s", "rcuf=s",
);

my $input1 = $opts{i1};
my $input2 = $opts{i2};
my $input3 = $opts{i3};
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


deal_NA_data( $input3, "__noNA_data3__" );
$input3 = "__noNA_data3__";

&check_data_v4( $input1  );
&check_data_v4( $input2  );
&check_data_v4( $input3  );
&check_data_sets_v2( $input1, $input2, $input3 );

copy($input1, '__matrix1__');
copy($input2, '__matrix2__');
copy($input3, '__matrix3__');

open RCMD, ">cmd.r";
print RCMD <<"RSCRIPT";

zd_pcor <- function (left = left , right =right  , confounding = confounding , unique_id_colname = "ID" , mannual_n_for_fdr = F, filtered_pairs_saved_in = "global_filtered_pairs"  ,
method = "spearman"
){
    if (exists(filtered_pairs_saved_in) == F) {
        global_filtered_pairs <<- data.frame()

    }

    left = left \%>\% as.data.frame() \%>\% rownames_to_column("ID")
    right = right \%>\% as.data.frame() \%>\% rownames_to_column("ID")
    confounding = confounding \%>\% as.data.frame() \%>\% rownames_to_column("ID")

    part_1 <- left[, unique_id_colname] \%>\% as.character()
    part_2 <- right[, unique_id_colname] \%>\% as.character()
    part_3 <- confounding[, unique_id_colname] \%>\% as.character ()

    pooled_ID <- c(part_1, part_2, part_3) \%>\% unique()

    neworder <- pooled_ID \%>\% order()
    pooled_ID <- pooled_ID[neworder] \%>\% data.frame(. , stringsAsFactors = F)
    colnames(pooled_ID) <- unique_id_colname

    left <- merge (pooled_ID , left , by = c(unique_id_colname) , all.x = T)
    right <- merge (pooled_ID , right , by = c(unique_id_colname) , all.x = T)
    confounding <- merge (pooled_ID , confounding , by = c(unique_id_colname) , all.x = T)

    left <- left[, - 1]
    right <- right[, - 1]
    confounding <- confounding[, - 1]


    for (i in 1 : ncol (left)) {
        if (i == 1) { pooled_r_and_p <- data.frame()
            pooled_error <- data.frame ()
        }
        for (j in 1 : ncol (right)) {
            pair_1 <- colnames (left)[i]
            pair_2 <- colnames (right)[j]

            temp_frame <- cbind (left[, i], right[, j], confounding)
            temp_frame <- na.omit(temp_frame)
            temp_frame <- apply (temp_frame , 2, function (each_col){
                result <- each_col \%>\% as.character () \%>\% as.numeric()
            })
            IS_error <<- F
            if (T) {
                tryCatch(each_line <- ppcor::pcor.test(x = temp_frame[, 1], y = temp_frame[, 2],
                z = temp_frame[, 3 : ncol(temp_frame)], method = c(method)),
                error = function(e) {
                    IS_error <<- T


                    each_filtered_pair <- data.frame (pair_1 = pair_1, pair_2 = pair_2)
                    global_filtered_pairs <<- rbind (global_filtered_pairs, each_filtered_pair)
                })
            }

            if (IS_error == F) {
                each_line <- data.frame (pair_1 = pair_1 , pair_2 = pair_2 , each_line , stringsAsFactors = F)
                colnames(each_line) <- c("pair_1", "pair_2", "r", "p", "statistic", "n", "gp", "Method")
                pooled_r_and_p <- rbind (pooled_r_and_p , each_line)
            }
        }
    }

    IS_filtered <- is.na (pooled_r_and_p[, "p"]) == T
    temp_filtered_data <- pooled_r_and_p[IS_filtered  , c("pair_1", "pair_2")]
    global_filtered_pairs <<- rbind (global_filtered_pairs , temp_filtered_data)
    pooled_r_and_p <- pooled_r_and_p[! IS_filtered,]

    if (mannual_n_for_fdr > 0) {n <- mannual_n_for_fdr} else {
        n <- nrow(pooled_r_and_p)}
    return (pooled_r_and_p)
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

confounderData <- read.table("__matrix3__", header = T, sep = "\t", check.names = FALSE,quote="", row.names = 1)
confounderData = t(confounderData)
if( zscal){
    confounderData = scale(confounderData)
}

dataSampleIds = rownames(data)
extraSampleIds = rownames(extraData)
confounderIds = rownames(confounderData)
interIds <- Reduce(intersect, list(dataSampleIds, extraSampleIds, confounderIds))

inData = data[interIds,]
extraData = extraData[interIds,]
confounderData =  confounderData[interIds,]


pcorData <- zd_pcor(left = inData, right = extraData, confounding = confounderData, method="$method")

allData <- tibble(Node1 = pcorData\$pair_1, Node2 = pcorData\$pair_2, r = pcorData\$r, P = pcorData\$p) \%>\%
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
            mutate_at(vars(- c("Node2")), function(x){
                ifelse(is.na(x), 0, x)
            }) \%>\%
            rename(` ` = Node2)

    
write_csv(corData,  "Inter_r_Matrix.csv")
write.table(corData, file = "Inter_r_Matrix.txt", quote = FALSE, sep = "\t", row.names=F)

pData <- allData \%>\%
            select(c("Node1", "Node2", "P")) \%>\%
            spread(Node1, "P") \%>\%
            mutate_at(vars(- c("Node2")), function(x){
                ifelse(is.na(x), 1, x)
            }) \%>\%
            rename(` ` = Node2)

write_csv(pData, "Inter_P_Matrix.csv")
write.table(pData, file = "Inter_P_Matrix.txt", quote = FALSE, sep = "\t", row.names=F)

fdrData <- allData \%>\%
            select(c("Node1", "Node2", "FDR")) \%>\%
            spread(Node1, "FDR") \%>\%
            mutate_at(vars(- c("Node2")), function(x){
                ifelse(is.na(x), 1, x)
            }) \%>\%
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




sub check_data_sets_v2 {
    my @files = @_;
    my @fs;
    my %hs1;
    my %hs2;
    foreach(@files){
        push @fs, basename $_;
        open IN, "$_";
        my $h = <IN>;
        my @h = split /\t/, $h;
        shift @h;
        foreach( @h ){
            $hs1{$_} += 1;
        }
        while(<IN>){
            chomp;
            my @s = split /\t/, $_;
            $hs2{$s[0]} += 1;
        }
        close IN;
    }
    
    my $cnt = 0;
    foreach(keys %hs1){
        if ( $hs1{$_} == 3) {
            $cnt ++;
        }
    }
    if ( $cnt < 2) {
        die "Error: The intersection samples of files: @fs  < 2 \n";
    }

    foreach(keys %hs2){
        if ( $hs2{$_} != 1) {
            die "Error: The feature \"$_\" is replicated in files: @fs \n";
        }
    }
}


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