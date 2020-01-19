#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use File::Copy;
use FindBin qw($Bin);
use Scalar::Util qw(looks_like_number);
my $r_path = "$Bin/../../bin/R/bin";
$ENV{'PATH'} = "$r_path;$ENV{'PATH'}";

my %opts;
GetOptions (\%opts,"i=s","o=s","c=s","html=s");

my $html = $opts{html};
my $input = $opts{i};
my $output = $opts{o};
my $coef = $opts{c};


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


&check_data( $input );
copy($input, '__matrix__');
open RCMD, ">cmd.r";
print RCMD <<"RSCRIPT";

tlog = read.table("__matrix__", header=T, quote="",com='', sep="\t", row.names=1, check.names=F)
ndt01<-data.frame()
for(i in 1: dim(tlog)[1]){
    value<-as.numeric(tlog[i,])
    QL <- quantile(value, probs = 0.25,na.rm=T)
    QU <- quantile(value, probs = 0.75,na.rm=T)
    QU_QL <- QU-QL
    coef <- $coef
    # QL;QU;QU_QL
    test01 <- value
    out_imp01 <- max(test01[which(test01 <= QU + coef*QU_QL)])
    test01[which(test01 > QU + coef*QU_QL)] <- out_imp01
    dt01<-as.data.frame(t(test01))
    ndt01<-rbind(ndt01,dt01)
}
ndt01<-t1<-as.matrix(ndt01)
rownames(ndt01) = rownames(tlog)
colnames(ndt01) = colnames(tlog)
write.table(ndt01,"out.table",quote=FALSE ,sep=\"\\t\")



RSCRIPT



close RCMD;
my $ret = system("R --restore --no-save < cmd.r");
if ( $ret ) {
        die "Error, died with $ret";
}
getRTable('out.table', $output);
unlink 'cmd.r';
unlink '__matrix__';
unlink 'out.table';


#get_html_table($output, $html, "y");
get_html_link_v2($html, "Outlier Processing Results", "$output", "tabH");

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

sub findParam{
	my $opt = shift;
	for(my $i = 0; $i <= $#ARGV; $i++){
		if ($ARGV[$i] eq $opt) {
			return $ARGV[$i+1];
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

sub isnumeric {  
        my $val = shift;  
        if(looks_like_number $val  || $val eq "NA"){
            return 1
        }else{
            return  0   
        }
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



