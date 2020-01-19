#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use File::Copy;
use FindBin qw($Bin);
use Data::Dumper;
use File::Basename;

my %opts;
GetOptions (\%opts,"matrix=s",  "output=s", "html=s", "log=s", "db=s", "ions=s", "mz_cut=s", "output2=s", "rt_cut=s",
);

my $matrix = $opts{matrix};
my $log = $opts{log};
my $output = $opts{output};
my $output2 = $opts{output2};
my $html = $opts{html};

my $db = $opts{db};
my $ions = $opts{ions};
my $mz_cut = $opts{mz_cut};
my $rt_cut = $opts{rt_cut};

open IN,"$Bin/adducts.txt";
my %mass2add;
while(<IN>){
    my @s = split /\s+/, $_;
    $mass2add{$s[0]}= $s[1];
}

close IN;

$ions = trim($ions);
my @adducts = split /,/, $ions;


## initilize db 

my $database;
open IN,"$db";
while(<IN>){
    chomp;
    $_ = trim($_);
    my @s = split /\t/, $_;

    if( @s == 2 ){
        my $int = int($s[1]);
        push @{ $database->{$int}}, [$s[0], $s[1], "NA"];
    } elsif ( @s == 3 ) {
        my $int = int($s[1]);
        push @{ $database->{$int}}, [$s[0], $s[1], $s[2]];
    }
}
close IN;

## search 

open IN,"$matrix";
open OUT, ">$output";
open LOG, ">$log";
my $h = <IN>;
chomp $h;
my @s = split /\t/, $h;
my $line = join "\t", @s[0,3..$#s];
print OUT "$line\n";
print LOG "Query\tQuery_mz\tQuery_rt\tBest_hit\tBest_hit_adducts_type\tBest_hit_mz\tBest_hit_delta_mz\tBest_hit_rt\tBest_hit_delta_rt\tOther_hits\tOther_hits_adducts_type\tOther_hits_mz\tOther_hits_delta_mz\tOther_hits_rt\tOther_hits_delta_rt\n";
while(<IN>){
    chomp;
    my @s = split /\t/, $_;
    my $q_rt = $s[2];
    my @q_mzs; my @q_names;
    #push @q_mzs, $s[1];
    foreach ( @adducts ) {
        my $v = $s[1] + $_;
        push @q_mzs, $v;
        push @q_names, $mass2add{$_};
    }

    my @cand;
    for(my $i = 0; $i <= $#q_mzs; $i++){
    #foreach my $q_mz ( @q_mzs){
        my $q_mz = $q_mzs[$i];
        my $q_type = $q_names[$i];
        my $q_int = int($q_mz);
        my @d_cands;
        my @v;
        push @v, ($q_int, $q_int - 1, $q_int + 1, $q_int + 2);
        foreach(@v){
            if( defined($database->{$_}) ){
            push @d_cands, @{ $database->{$_} };
            }
        }

        for(my $i = 0; $i <= $#d_cands; $i++){
            my @pairs = @{$d_cands[$i]};
            my $diff = abs($pairs[1]-$q_mz);
            my $diff2 = $pairs[1]-$q_mz;
            if( $diff <= $mz_cut ){
                if( $rt_cut >= 0 && $pairs[2] ne "NA" && abs($pairs[2]-$q_rt) <= $rt_cut ){
                    my $delta_rt = $pairs[2]-$q_rt;
                    push @cand, [ $pairs[0], $diff, $pairs[1], $pairs[2], $q_type, $diff2 , $delta_rt];
                }elsif( $rt_cut < 0 || $pairs[2] eq "NA" ){
                    push @cand, [ $pairs[0], $diff, $pairs[1], $pairs[2], $q_type, $diff2, "NA"];
                }
            }
        }
    }
    if( @cand != 0){
        my @sort = sort { $a->[1] <=> $b->[1] } @cand;
        my $q = $s[0];
        $s[0] = $sort[0] -> [0];
        my $l = join "\t", @s[0,3..$#s];
        print OUT "$l\n";
        print LOG "$q\t$s[1]\t$q_rt";
        if( @sort > 5){
            @sort = @sort[0..4];
        }
        foreach  (@sort) {
            my $hit = $_ -> [0];
            my $mz = $_ -> [2];
            my $rt = $_ -> [3];
            my $type = $_ -> [4];
            my $delta_mz = $_ -> [5];
            my $delta_rt = $_ -> [6];
            print LOG "\t$hit\t$type\t$mz\t$delta_mz\t$rt\t$delta_rt";
        }
        print LOG "\n";
    }else{
        my $l = join "\t", @s[0,3..$#s];
        print OUT "$l\n";
        print LOG "$s[0]\t$s[1]\t$s[2]\tNo Hit\n";
    }
}

close IN;
close OUT;
close LOG;



get_uniq_table($output, $output2);
get_html_link_v2($html, "LC-MS Library Search Results", $output, "tabH",  $output2, "tabH", $log, "tabH");

sub trim
{
        my $string = shift;
        $string =~ s/^\s+//;
        $string =~ s/\s+$//;
        return $string;
}




sub get_uniq_table{
    my $in =shift;
    my $out = shift;
    open IN,"$in";
    open OUT,">$out";
    my $h = <IN>;
    print OUT "$h";

    my %hs;
    my %sum;
    while(<IN>){
        chomp;
        my @s = split /\t/,$_;
        my $sum = 0;
        foreach  (@s[1..$#s]) {
            if ( $_ eq "NA") {
                next;
            }else{
                $sum += $_;
            }
        }
        if( defined($hs{$s[0]}) ){
            if( $sum >= $sum{$s[0]}){
                $hs{$s[0]} = $_;
                $sum{$s[0]} = $sum;    
            }
        
        }else{
            $hs{$s[0]} = $_;
            $sum{$s[0]} = $sum;
        }
    }

    foreach  (values %hs) {
        print OUT "$_\n";
    }

    close IN;
    close OUT;
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


use Scalar::Util qw(looks_like_number);
sub isnumeric {  
        my $val = shift;  
        looks_like_number $val;
}  