#!/usr/bin/perl -w
use strict;
use warnings;

use Getopt::Long;
use File::Copy;
use FindBin qw($Bin);
my $r_path = "$Bin/../../bin/R/bin";
$ENV{'PATH'} = "$r_path;$ENV{'PATH'}";

my %opts;
GetOptions (\%opts,"i=s","o=s","g=s","html=s", "corrected=s",   "o2=s",  "thrN=s",
);


my $input = $opts{i};
my $group_file = $opts{g};
my $output = $opts{o};
my $output_sign = $opts{o2};
my $html = $opts{html};

my $multi_test_method = $opts{corrected};
my $thrN = $opts{thrN};


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




check_data_v4($input, $group_file);
copy($input, '__matrix__');
copy($group_file, '__group__');
open RCMD, ">cmd.r";
print RCMD <<"RSCRIPT";
std <- function(x) sd(x) / sqrt(length(x))
data = read.table("__matrix__", header=T, quote="",com='', sep="\t", row.names=1, check.names=F)
group <- read.table(quote="", com='', file="__group__", header = F, sep = "\t",check.names=FALSE, colClasses=c("character","character"))
uniq.group <- as.character(unique(group\$V2))
group1_samples <- subset(group, V2 == uniq.group[1])\$V1
group1_samples <- as.character(group1_samples)

mean = apply(data[, group1_samples], 1, mean)
var = apply(data[, group1_samples], 1, var)
stderr = apply(data[, group1_samples], 1, std)
okk = data.frame(mean, var, stderr)
rownames(okk)=rownames(data)
names(okk)[names(okk) == "mean"] = paste("mean(", uniq.group[1], ")", sep = "")
names(okk)[names(okk) == "var"] = paste("variance(", uniq.group[1], ")", sep = "")
names(okk)[names(okk) == "stderr"] = paste("stderr(", uniq.group[1], ")", sep = "")


for( i in 2:length(uniq.group) ){
    subgroup_samples <- subset(group, V2 == uniq.group[i])\$V1
    subgroup_samples <- as.character(subgroup_samples)
    mean = apply(data[, subgroup_samples], 1, mean)
    var = apply(data[, subgroup_samples], 1, var)
    stderr = apply(data[, subgroup_samples], 1, std)
    okk = cbind(okk, mean)
    okk = cbind(okk, var)
    okk = cbind(okk, stderr)
    names(okk)[names(okk) == "mean"] = paste("mean(", uniq.group[i], ")", sep = "")
    names(okk)[names(okk) == "var"] = paste("variance(", uniq.group[i], ")", sep = "")
    names(okk)[names(okk) == "stderr"] = paste("stderr(", uniq.group[i], ")", sep = "")
}

for (i in 1 : nrow(data)) {
    varVn = vector()
    facFcVn = vector()
    for( j in 1:length(uniq.group) ){
        subgroup_samples <- subset(group, V2 == uniq.group[j])\$V1
        subgroup_samples <- as.character(subgroup_samples)
        varVn = c( varVn, as.numeric(data[i, subgroup_samples]) )
        facFcVn = c(facFcVn, rep(j, length(subgroup_samples)))
    }
    okk\$p[i] = tryCatch(
    {
        row.aov = aov(varVn ~ facFcVn)
        summary(row.aov)[[1]][1, "Pr(>F)"]
    }, error = function(e){
        1
    }
    )
}
okk\$p_cor = p.adjust(okk\$p, method = "$multi_test_method")
okk = okk[order(okk[, "p"]),]
names(okk)[names(okk) == "p_cor"] = "$multi_test_method"
write.table(okk,"out.all",quote=FALSE ,sep="\\t")
okk = subset(okk, $multi_test_method <= $thrN)
write.table(okk,"out.filter",quote=FALSE ,sep="\\t")

RSCRIPT





close RCMD;
my $ret = system("R --restore --no-save < cmd.r");
if ( $ret ) {
        die "Error, died with $ret";
}
getRTable_v2('out.all', $output);
getRTable_v2('out.filter', $output_sign);
get_html_link_v2($html, "Analysis of variance Results", "$output", "tabH", "$output_sign", "tabH");

unlink 'cmd.r';
unlink '__matrix__';
unlink '__group__';
unlink 'out.all';
unlink 'out.filter';








# getRTable(in, out)
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
