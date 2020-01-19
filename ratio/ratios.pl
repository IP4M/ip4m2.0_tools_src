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
GetOptions (\%opts, "i=s",    "html=s",  "lib=s", 
);

my $input = $opts{i};
my $html = $opts{html};
my $lib = $opts{lib};
my $db = "$Bin/../pathways/hmdb_info.txt";


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


check_data($input);

sub get_greek_hs {
    my %hs;
    open IN,"$Bin/../pathways/replace_info.txt";
    while(<IN>){
        chomp;
        my @s = split /\t/, $_;
        $hs{$s[0]} = $s[1];
    }
    close IN;
    return %hs;
}

my %gks = get_greek_hs();
my @greeks = keys %gks;

sub process_greek {
    my $string = shift;
    foreach (@greeks) {
        $string=~s/$_/$gks{$_}/g;
    }
    return $string;
}

sub process_name {
    my $string = shift;
    $string =~ s/^\s+//;
    $string =~ s/\s+$//;
    $string =~ s/\s+/ /g;
    $string =~ s/_/ /g;
    $string =~ s/\?//g;
    $string = process_greek($string);
    $string = lc($string);
    return $string;
}


my %database;
open IN, "$db";
<IN>;
while(<IN>){
    chomp;
    my @s = split /\t/, $_;
    for(my $i = 0; $i <= 11; $i++){
        unless(defined($s[$i])){
            $s[$i] = "";
        }
    }   
    my $hit = $s[8];
    if(  $hit=~/^C/){
        foreach( @s[1, 2, 3] ){
            $_ = lc($_);
            $database{$_} = $hit;
        }
        my @ss = split / \|\| /, $s[10];
        foreach( @ss ){
            $_ = lc($_);
            $database{$_} = $hit;
        }
    }
}
close IN;


my %cnum;
my %max;
my @aa;
open IN,"$input";
my $header = <IN>;
chomp $header;
while (<IN>) {
    chomp;
    my @s = split /\t/, $_;
    my $name = process_name($s[0]);
    my $sum = sum_row( @s[1..$#s] );

    if(defined($database{$name})){
        my $hit = $database{$name};
        unless(defined($max{$hit})){
            @{$cnum{$hit}} = @s;
            $max{$hit} = $sum;
        }else{
            if ( $sum > $max{$hit}) {
                @{$cnum{$hit}} = @s;
                $max{$hit} = $sum;
            }
        }
    }
    push @aa, $_;
}
close IN;


foreach(keys %cnum){
    my @dats = @{$cnum{$_}};
    print "$_\n";
    print "@dats\n";
}

my %reacs;
my $reac_file = "$Bin/reactions.txt";
open IN,"$reac_file";
<IN>;
while(<IN>){
    chomp;
    my @s = split /\t/, $_;
    if( $s[0] eq $lib ){
        if(defined($cnum{$s[3]}) && defined($cnum{$s[4]})){
            push @{$reacs{$s[3]}{$s[4]}}, $_;
        }
    }
}
close IN;

my $html_dir = dirname($html);
my $info_file = "$html_dir/reactions_info.txt";
open OUT,">$info_file";
print OUT "Pathway_Name\tReaction_Name\tSubstrate\tProduct\n";
foreach my $sub ( keys %reacs) {
    my %hs = %{$reacs{$sub}};
    foreach my $pro ( keys %hs) {
        my @reacs = @{$hs{$pro}};
        foreach(@reacs){
            my @s = split /\t/, $_;
            my $name1 = $cnum{$s[3]} -> [0];
            my $name2 = $cnum{$s[4]} -> [0];
            print OUT "$s[1]\t$s[2]\t$s[3]($name1)\t$s[4]($name2)\n";
        }
    }
}
close OUT;

my $ratio_file = "$html_dir/Ratio.txt";
open OUT,">$ratio_file";
print OUT "$header\n";
my @new_datas;
foreach my $sub ( keys %reacs) {
    my %hs = %{$reacs{$sub}};
    foreach my $pro ( keys %hs) {
        my @sub_data = @{$cnum{$sub}};
        my @pro_data = @{$cnum{$pro}};
        my @new_data;
        $new_data[0] = "$pro_data[0] / $sub_data[0]";
        for(my $i = 1; $i <= $#sub_data; $i++){
            if( $sub_data[$i] != 0 ){
                $new_data[$i] = $pro_data[$i] / $sub_data[$i];
            }else{
                $new_data[$i] = "0";
            }
        }
        my $line = join "\t", @new_data;
        push @new_datas, $line;
        print OUT "$line\n";
    }
}
close OUT;

my $all_file = "$html_dir/AllMet_and_Ratio.txt";
open OUT,">$all_file";
print OUT "$header\n";
foreach( @aa ){
    print OUT "$_\n";
}
foreach( @new_datas ){
    print OUT "$_\n";
}



close OUT;

get_html_link_v2($html, "Variable Expansion Results",  
"$ratio_file", "tabH", 
"$all_file", "tabH", 
"$info_file", "tabH", 
);







sub sum_row{
    my @data = @_;
    my $sum;
    foreach(@data){
        $sum += $_;
    }
    return $sum;
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
                get_html_table_v2($file, "$out_dir/$name.html", "y");
            }elsif(  $type eq "tabN" ){
                get_html_table_v2($file, "$out_dir/$name.html", "n");
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



sub get_html_table_v2{
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
