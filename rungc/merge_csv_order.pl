#!/usr/bin/perl
use strict;
use warnings;
use File::Copy;
use File::Path;
use File::Basename;
use FindBin qw($Bin);
use Getopt::Long;
#use Text::CSV;
use Try::Tiny;


use Scalar::Util qw(looks_like_number);
sub isnumeric {  
        my $val = shift;  
        looks_like_number $val;
}  

my %opts;
GetOptions (\%opts, "i=s", "peaktable=s", "peaktable2=s",  "html=s", 
);

my $usage = <<"USAGE";
       Program : $0
       Version : 1.0
       Discription: 
       Usage :perl $0 [options]
                   -i  input leco csv files
                   -o  output html file
USAGE

die $usage if ( !$opts{i} );



use Text::CSV;
my $csv = Text::CSV->new({ binary => 1 });

sub get_term{
    my $file = shift;
    open IN,"$file";
    my $term = <IN>;
    chomp $term;
    close IN;
    return $term;
}



my $input = get_term($opts{i});
my $details = $opts{peaktable};
my $html = $opts{html};
my $peaktable = $opts{peaktable2};

my @files = split /__SEP__/, $input;
my @names;
foreach ( @files ) {
    my $basename = basename $_;
    my @s = split /\./, $basename;
    pop @s;
    my $name = join ".", @s;
    push @names, $name;
}

my %info;
for(my $i = 0; $i <= $#files; $i++){
    my $file = $files[$i];
    my $name = $names[$i];
    open IN,"$file";
    my $h = <IN>;

    chomp $h;
    my @s;
    if( $csv->parse($h) ){
        @s = $csv->fields();
    }else{
        my $log =  $csv->error_diag();
        die "Error: \"$name\", it is not a LECO csv file in the inputs.\n$log\n$h\n";
    }

    my $compound_index = get_index( "Name", \@s, $name  );
    my $rt_index = get_index( "R.T.", \@s, $name  );
    my $mass_index;
    try{
        $mass_index = get_index( "Quant Masses", \@s, $name  );
    }catch{
        $mass_index = get_index( "UniqueMass", \@s, $name  );     
    };
    my $area_index = get_index( "Area", \@s, $name  );

    while(<IN>){
        chomp;
        my @s;
        if( $csv->parse($_) ){
            @s = $csv->fields();
        }else{
            my $log =  $csv->error_diag();
            die "Error: \"$name\", it is not a LECO csv file in the inputs.\n$log\n$_\n";
        }
        my $compound = $s[$compound_index];
        my $rt = $s[$rt_index];
        if($rt=~/^\s*$/){
            $rt = 0;
        }else{
            my @s = split /,/, $rt;
            $rt = trim($s[0]);
        }
        my $mass = $s[$mass_index];
        if($mass=~/^\s*$/){
            $mass = "-";
        }
        my $area = $s[$area_index];
        if($area=~/^\s*$/){
            $area = 0;
        }
        push @{$info{$mass}{$name}}, [$compound, $rt, $area];
    }
    close IN;
}

my @mass = keys %info;
my $flag = 1;
foreach(@mass){
    unless( isnumeric($_)){
        $flag = 0;
    }
}
if($flag == 1){
    @mass = sort{ $a <=> $b } @mass;
}else{
    @mass = sort{ $a cmp $b } @mass;
}

open OUT,">details";
my @lines;
foreach( @names ){
    push @lines, "$_.name\t$_.rt\t$_.area";
}
my $l = join "\t", @lines;
print OUT "Mass\t$l\n";
foreach my $mass ( @mass ) {
    my %data = %{$info{$mass}};

    my $num = @{$data{$names[0]}};
    foreach my $name ( @names ) {
        my @array = @{$data{$name}};
        if($num != @array){
            die "Error: for mass $mass, row numbers are not same among samples\n";
        }
        @{$data{$name}} = sort{ $a->[1] <=> $b->[1] } @array;
    }
    for(my $i = 0 ; $i < $num; $i++){
        print OUT "$mass";
        foreach my $name (@names) {
            my $cmp = $data{$name} -> [$i] -> [0];
            my $rt = $data{$name} -> [$i] -> [1];
            my $area = $data{$name} -> [$i] -> [2];
            print OUT "\t$cmp\t$rt\t$area";
        }
        print OUT "\n";       
    }
}
close OUT;


open IN,"details";
open OUT,">$details";
my $h = <IN>;
chomp $h;
my @s = split /\t/, $h;
shift @s;
$h = join "\t", @s;
print OUT "ID\tMass\trt.mean\t$h\n";
my @data;
while (<IN>) {
    chomp;
    my @s = split /\t/, $_;
    my $mass = shift @s;
    my $sum_rt = 0;
    my $counts = 0;
    for(my $i = 0; $i <= $#s; $i+=3){
        if( $s[$i+1] ne "-"){
            $sum_rt +=  $s[$i+1];
            $counts += 1;
        }
    }
    my $mean_rt;
    if($counts == 0){
        $mean_rt = "-";
    } else {
       $mean_rt = $sum_rt / $counts;
    }
    push @data, [$mass, $mean_rt, @s];
}

if($flag == 1){
    @data = sort {  $a->[0] <=> $b->[0] or  $a->[1] <=> $b->[1] } @data;
}else{
    @data = sort {  $a->[0] cmp $b->[0] or  $a->[1] <=> $b->[1] } @data;
}

for(my $i = 0; $i<= $#data; $i++){
    my $id = $i + 1;
    my @aa = @{$data[$i]};
    my $tmp = join "\t", @aa;
    print OUT "$id\t$tmp\n";
}


close IN;
close OUT;



open IN,"$details";
open OUT,">$peaktable";
$h = <IN>;
chomp $h;
@s = split /\t/, $h;
print OUT "name";
for(my $i = 3; $i <= $#s; $i+=3){
    my $name = $s[$i];
    $name=~s/\.name$//;
    print OUT "\t$name";
}
print OUT "\n";
while (<IN>) {
    chomp;
    my @s = split /\t/, $_;

    my %freq;
    my %sum;
    my %aver;
    for(my $i = 3; $i <= $#s; $i+=3){
        if( $s[$i+1] ne "-"){
            $freq{$s[$i]} += 1;
            $sum{$s[$i]} += $s[$i+2];
        }
    }
    my @names = keys %freq;
    foreach(@names){
        $aver{$_} = $sum{$_} / $freq{$_};
    }

    @names = sort { $freq{$b} <=> $freq{$a} or $aver{$b}  <=>  $aver{$a}} @names;  
    my $name = $names[0];
    unless(defined($name)){
        $name = $s[3];
    }

    print OUT "$name";
    for(my $i = 5; $i <= $#s; $i+=3){
        print OUT "\t$s[$i]";
    }
    print OUT "\n";
}
close IN;
close OUT;

get_html_link_v2($html, "LECO CSV Merge Results", $peaktable, "tabH",  $details, "tabH",);
unlink 'details';

sub get_cls{
    my $relation = shift;
    my $mark = shift;
    my $i = shift;
    my @cls;
    unless (defined( $mark->{$i} )) {
        push @cls, $i;
        $mark->{$i} = 1;
        my @nexts = keys %{$relation->{$i}};
        foreach (@nexts) {
            push @cls, &get_cls($relation, $mark, $_);
        }
    }
    return @cls;
}


## name, array
sub get_index{
    my $string = shift;
    my $array = shift;
    my $file = shift;
    my $index = -1;
    my @aa = @{$array};
    for(my $i = 0; $i <= $#aa; $i++){
        if($aa[$i] =~/$string/){
            $index = $i;
        }
    }
    if( $index == -1){
        die "Error: \"$file\", it is not a LECO csv file in the inputs\n";
    }
    return $index;
}

#source dest
sub mklink{
    my $source = shift;
    my $destination = shift;
    my $cmd = "mklink \"$destination\" \"$source\"";
    &run_cmd( $cmd );
}


sub run_cmd {
        my $cmd = shift;
#        print "$cmd\n";

        my $ret = system( $cmd );
        if ( $ret ){
                die "Error, died with $ret";
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



sub trim
{
        my $string = shift;
        $string =~ s/^\s+//;
        $string =~ s/\s+$//;
        return $string;
}