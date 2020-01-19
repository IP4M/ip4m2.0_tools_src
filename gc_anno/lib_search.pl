#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use File::Copy;
use FindBin qw($Bin);
use Data::Dumper;
use File::Basename;

my %opts;
GetOptions (\%opts,"matrix=s", "query=s", "db=s", "method=s", "cutoff=s", "output=s", "html=s", "mzRes=s", "rt_cut=s", "log=s","output2=s", "user=s",
);

my $matrix = $opts{matrix};
my $query = $opts{query};
my $db = $opts{db};
my $log = $opts{log};
my $output = $opts{output};
my $output2 = $opts{output2};
my $html = $opts{html};

my $method = $opts{method};
my $cutoff = $opts{cutoff};
my $mzRes = $opts{mzRes};
my $rt_cut = $opts{rt_cut};
my $user = $opts{user};

my %id2name;
my %id2db;
my $id = 1;


my ($query_data, $query_rt )= readMSP($query, $mzRes);

my ($db_data, $db_rt);
if( $user eq "yes"){
    ($db_data, $db_rt) = readMSP2($db, $mzRes, "user");
}else{
    my @files = split /,/, $db;
    my $cfg_file = "$Bin/gcdb_cfg.txt";
    my %file2db = get_cfg($cfg_file);
    foreach( @files){
        my $file = "$Bin/$_";
        my $name = $file2db{$_};
        my ($data, $rt) = readMSP2($file, $mzRes, $file2db{$_});
        merge_hash1( $data);
        merge_hash2( $rt);
    }
}


sub merge_hash1{
    my $hs = shift;
    my @key = keys %{$hs};
    foreach(@key){
        $db_data -> {$_} = $hs -> {$_};
    }
}

sub merge_hash2{
    my $hs = shift;
    my @key = keys %{$hs};
    foreach(@key){
        $db_rt -> {$_} = $hs -> {$_};
    }
}


sub get_cfg{
    my $in = shift;
    my %hs;
    open IN,"$in";
    while(<IN>){
        chomp;
        my @s = split /\t/, $_;
        $hs{$s[0]} = $s[1];
    }
    close IN;
    return %hs;
}


sub get_query_ms{
    my $query = shift;
    my %hs = %{$query_data->{$query}};
    my @mz = keys %hs;
    @mz = sort { $a <=> $b} @mz;
    my $tmp = " ";
    foreach(@mz) {
        $tmp .= "$_ $hs{$_}; ";
    }
    return $tmp;
}

sub get_db_ms{
    my $query = shift;
    my %hs = %{$db_data->{$query}};
    my @mz = keys %hs;
    @mz = sort { $a <=> $b} @mz;
    my $tmp = " ";
    foreach(@mz) {
        $tmp .= "$_ $hs{$_}; ";
    }
    return $tmp;
}




my $query_sos = get_sum_of_squares($query_data);
my $db_sos = get_sum_of_squares($db_data);




open IN,"$matrix";
open OUT,">$output";
open LOG,">$log";
my $h = <IN>;
print OUT "$h";
print LOG "Query\tDB\tmatchFactor\trt_diff\tquery_rt\tquery_ms\tdb_rt\tdb_ms\tdb_name\n";
while(<IN>){
    chomp;
    my @ar = split /\t/, $_;
    my $info_line = join "\t", @ar[1..$#ar];
    unless( defined($query_data->{$ar[0]}) ){
        print OUT "$ar[0]\t$info_line\n";
        print LOG "$ar[0]\tNo Mass Spectrum\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n";
    }else{
        my $query_name = $ar[0];
        my @dbs = keys %{$db_data};
        if( $rt_cut < 0 || $query_rt->{$query_name} == -1 ){
            my ($best_hit, $best_score) = get_anno($query_name, $method, \@dbs);
            $best_score = sprintf( "%.10f", $best_score);
            if( $best_score >= $cutoff){
                my $query_rt_ = "NA";
                unless( $query_rt->{$query_name} == -1){
                    $query_rt_ = $query_rt->{$query_name};
                }
                my $query_ms_ = get_query_ms($query_name);
                my $db_rt_ = "NA";
                unless( $db_rt->{$best_hit} == -1){
                    $db_rt_ = $db_rt->{$best_hit};
                }
                my $db_ms_ = get_db_ms($best_hit);
                my $db_name = $id2db{$best_hit};
                $best_hit = $id2name{$best_hit};
                
                print OUT "$best_hit\t$info_line\n";
                print LOG "$ar[0]\t$best_hit\t$best_score\tNA\t$query_rt_\t$query_ms_\t$db_rt_\t$db_ms_\t$db_name\n";

            }else{
                my $query_rt_ = "NA";
                unless( $query_rt->{$query_name} == -1){
                    $query_rt_ = $query_rt->{$query_name};
                }
                my $query_ms_ = get_query_ms($query_name);
                print OUT "$ar[0]\t$info_line\n";
                print LOG "$ar[0]\tNo Hit\tNA\tNA\t$query_rt_\t$query_ms_\tNA\tNA\tNA\n";          
            }
        }else{
            my @candi;
            foreach ( @dbs ) {
                if ( $db_rt->{$_} == -1 ) {
                    push @candi, $_;
                }else{
                    if( abs( $query_rt->{$query_name} - $db_rt->{$_} ) <= $rt_cut ){
                        push @candi, $_;
                    }
                }
            }
            if (@candi <= 0) {
                my $query_rt_ = "NA";
                unless( $query_rt->{$query_name} == -1){
                    $query_rt_ = $query_rt->{$query_name};
                }
                my $query_ms_ = get_query_ms($query_name);
                print OUT "$ar[0]\t$info_line\n";
                print LOG "$ar[0]\tNo Hit\tNA\tNA\t$query_rt_\t$query_ms_\tNA\tNA\tNA\n";  
            }else{
                my ($best_hit, $best_score) = get_anno($query_name, $method, \@candi);
                my $diff = "NA";
                unless ( $db_rt->{$best_hit} == -1 ) {
                      $diff = abs( $query_rt->{$query_name} - $db_rt->{$best_hit} );
                }
              
                
                ##print "$cutoff\n";
                $best_score = sprintf( "%.10f", $best_score);
                if( $best_score >= $cutoff){
                    my $query_rt_ = "NA";
                    unless( $query_rt->{$query_name} == -1){
                        $query_rt_ = $query_rt->{$query_name};
                    }
                    my $query_ms_ = get_query_ms($query_name);
                    my $db_rt_ = "NA";
                    unless( $db_rt->{$best_hit} == -1){
                        $db_rt_ = $db_rt->{$best_hit};
                    }
                    my $db_ms_ = get_db_ms($best_hit);
                    my $db_name = $id2db{$best_hit};
                    $best_hit = $id2name{$best_hit};
                    print OUT "$best_hit\t$info_line\n";
                    print LOG "$ar[0]\t$best_hit\t$best_score\t$diff\t$query_rt_\t$query_ms_\t$db_rt_\t$db_ms_\t$db_name\n";
                }else{
                    my $query_rt_ = "NA";
                    unless( $query_rt->{$query_name} == -1){
                        $query_rt_ = $query_rt->{$query_name};
                    }
                    my $query_ms_ = get_query_ms($query_name);
                    ##
                    ##print "$ar[0]\t$best_hit\t$best_score\t$diff\n";
                    print OUT "$ar[0]\t$info_line\n";
                    print LOG "$ar[0]\tNo Hit\tNA\tNA\t$query_rt_\t$query_ms_\tNA\tNA\tNA\n";          
                }
            }
        }
    }
}

close IN;
close OUT;
close LOG;

#&get_html_link($html, "Library Search", $output, $log);
get_uniq_table($output, $output2);
&get_html_link_v2($html, "GC-MS Library Search Results", $output, "tabH", $output2, "tabH", $log, "tabH");










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
            $sum += $_;
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


sub get_anno{
    my $query_name = shift;
    my $method = shift;
    my $db_names = shift;
    my @dbs = @{$db_names};
    
    my $q_sos = $query_sos -> {$query_name};
    my %q_ms = %{$query_data -> {$query_name}};
    ##print LOG "q_sos:  $q_sos\n";
    my @score;
    foreach my $db_name ( @dbs) {
        my $d_sos = $db_sos -> {$db_name};
        my %d_ms = %{$db_data -> {$db_name}};
        ##print LOG "d_sos:  $d_sos\n";
        if ( $q_sos == 0 || $d_sos == 0 ){
            if( $q_sos == 0 && $d_sos == 0 ){
                push @score, [$db_name, 1];
            } else {
                push @score, [$db_name, 0];
            }
        } else {
            if ( $method eq "dot") {
                my %ms;
                my @q_ms = keys %q_ms;
                foreach (@q_ms) {
                    $ms{$_} = 1;
                }
                my @d_ms = keys %d_ms;
                foreach( @d_ms){
                    $ms{$_} = 1;
                }
                my $pro = 0;
                foreach (keys %ms) {
                    ##print LOG "ms: $_\n";
                    ##print LOG "qms-v: $q_ms{$_}\n";
                    ##print LOG "dms-v: $d_ms{$_}\n";
                    if( defined($q_ms{$_})  && defined($d_ms{$_}) ){
                        $pro = $pro + $q_ms{$_} * $d_ms{$_} 
                    }
                }
                my $mf = $pro / (sqrt($q_sos) * sqrt($d_sos));
                push @score, [$db_name, $mf];
                ##print "$query_name\t$db_name\t$mf\n";
            }else{
                my %ms;
                my @q_ms = keys %q_ms;
                foreach (@q_ms) {
                    $ms{$_} = 1;
                }
                my @d_ms = keys %d_ms;
                foreach( @d_ms){
                    $ms{$_} = 1;
                }
                my $pro = 0;
                foreach (keys %ms) {
                    ##print LOG "ms: $_\n";
                    ##print LOG "qms-v: $q_ms{$_}\n";
                    ##print LOG "dms-v: $d_ms{$_}\n";
                    if( defined($q_ms{$_})  && defined($d_ms{$_}) ){
                        $pro = $pro + ($q_ms{$_}/sqrt($q_sos) - $d_ms{$_}/sqrt($d_sos))**2;
                        ##print LOG "run 1\n";
                    }elsif( defined($q_ms{$_})  && !defined($d_ms{$_})  ){
                        $pro = $pro + ( $q_ms{$_}/sqrt($q_sos) )**2;
                        ##print LOG "run 2\n";
                    }elsif( !defined($q_ms{$_})  && defined($d_ms{$_})  ){
                        $pro = $pro + ( 0 - $d_ms{$_}/sqrt($d_sos))**2;
                        ##print LOG "run 3\n";
                    }
                     
                }
                my $mf = 1 / ( 1 + $pro);
                ##print "$query_name\t$db_name\t$mf\n";
                push @score, [$db_name, $mf];
            }
        }
    }
    my @sort = sort { $b->[1] <=> $a->[1] } @score;
    #print "$query_name\t$sort[0][0]\t$sort[0][1]\n";
    return ($sort[0][0], $sort[0][1]);
}











sub get_sum_of_squares{
    my $data = shift;
    my %sum;
    foreach my $name ( keys %{$data}) {
        my @intens = values %{$data -> {$name}};
        my $sum;
        foreach  (@intens) {
            my $product = $_ * $_;
            $sum += $product;
        }
        $sum{$name} = $sum;
    }
    return \%sum;
}




sub readMSP{
    my ( $msp_file, $mzRes ) = @_ ;
    my %data;
    my %rt;
    open (MSP , "<" , $msp_file) or die $! ;
    local $/ = 'Name:';
    <MSP>;
    while(<MSP>){
        chomp;
        my @infos = split /\n/ , $_;
        my $name = $infos[0];
        $name=~s/^\s+//;
        $name=~s/\s+$//;
        my $rt = "-1";
        my @temp_mzs;
        my @temp_intensities;
        for(my $i=0 ; $i<@infos ; $i++) {
            if ($infos[$i] =~ /(\d+\.?\d*)\s+(\d+\.?\d*)\s*;\s*/) {
                my @ions = split ( /;/ , $infos[$i] );
                foreach my $ion (@ions) {
                    if ($ion =~ /^\s*(\d+\.?\d*)\s+(\d+\.?\d*)$/) {
                        my $mz = $1;
                        my $intensity =$2 ;
                        push ( @temp_intensities , $intensity ) ;
                        # Truncate/round mzs depending on $mzRes wanted
                        if ($mzRes == 0) {
                            my $mz_rounded = sprintf("%.".$mzRes."f", $mz) ;
                            push (@temp_mzs , $mz_rounded) ;
                        }
                        # Check that $mzRes is not greater than the number of digits after comma
		    			elsif ($mzRes > 0) {
                            if ($mz !~ /^\d+\.\d+$/) { die "*********\n\nYou are trying to specify $mzRes significant decimals, but one or more masses in the input file are unitary masses.\nYou should try again with mzRes = 0\n\n\n"; }
                            elsif($mzRes > length(( $mz =~ /.+\.(.*)/)[0] )) {
                                $mz = sprintf("%.".$mzRes."f" , $mz) ;
		    				}
                            my $mz_rounded = _round_num($mz,$mzRes) ;
                            push (@temp_mzs , $$mz_rounded) ;
		    			}
		    		}
		    	}
            }
            if ($infos[$i]=~/^rt:/) {
                $rt = (split /:/, $infos[$i])[1];
                $rt=~s/^\s+//;
                $rt=~s/\s+$//;
            }
        }
        for(my $i = 0; $i <= $#temp_mzs; $i++){
            $data{$name}{$temp_mzs[$i]} = $temp_intensities[$i];
            #$data{$id}{$temp_mzs[$i]} = $temp_intensities[$i];
        }
        $rt{$name} = $rt;
        #$rt{$id} = $rt;
        #$id++;
        #$id2name{$id} = $name;
    }
    close MSP;
    return (\%data, \%rt);
}


sub readMSP2{
    my ( $msp_file, $mzRes, $db_name ) = @_ ;
    my %data;
    my %rt;
    open (MSP , "<" , $msp_file) or die $! ;
    local $/ = 'Name:';
    <MSP>;
    while(<MSP>){
        chomp;
        my @infos = split /\n/ , $_;
        my $name = $infos[0];
        $name=~s/^\s+//;
        $name=~s/\s+$//;
        my $rt = "-1";
        my @temp_mzs;
        my @temp_intensities;
        for(my $i=0 ; $i<@infos ; $i++) {
            if ($infos[$i] =~ /(\d+\.?\d*)\s+(\d+\.?\d*)\s*;\s*/) {
                my @ions = split ( /;/ , $infos[$i] );
                foreach my $ion (@ions) {
                    if ($ion =~ /^\s*(\d+\.?\d*)\s+(\d+\.?\d*)$/) {
                        my $mz = $1;
                        my $intensity =$2 ;
                        push ( @temp_intensities , $intensity ) ;
                        # Truncate/round mzs depending on $mzRes wanted
                        if ($mzRes == 0) {
                            my $mz_rounded = sprintf("%.".$mzRes."f", $mz) ;
                            push (@temp_mzs , $mz_rounded) ;
                        }
                        # Check that $mzRes is not greater than the number of digits after comma
		    			elsif ($mzRes > 0) {
                            if ($mz !~ /^\d+\.\d+$/) { die "*********\n\nYou are trying to specify $mzRes significant decimals, but one or more masses in the input file are unitary masses.\nYou should try again with mzRes = 0\n\n\n"; }
                            elsif($mzRes > length(( $mz =~ /.+\.(.*)/)[0] )) {
                                $mz = sprintf("%.".$mzRes."f" , $mz) ;
		    				}
                            my $mz_rounded = _round_num($mz,$mzRes) ;
                            push (@temp_mzs , $$mz_rounded) ;
		    			}
		    		}
		    	}
            }
            if ($infos[$i]=~/^rt:/) {
                $rt = (split /:/, $infos[$i])[1];
                $rt=~s/^\s+//;
                $rt=~s/\s+$//;
            }
        }
        for(my $i = 0; $i <= $#temp_mzs; $i++){
            #$data{$name}{$temp_mzs[$i]} = $temp_intensities[$i];
            $data{$id}{$temp_mzs[$i]} = $temp_intensities[$i];
        }
        #$rt{$name} = $rt;
        $rt{$id} = $rt;
        $id2name{$id} = $name;
        $id2db{$id} = $db_name;
        $id++;
    }
    close MSP;
    return (\%data, \%rt);
}

sub _round_num {
    ## Retrieve Values
    my ( $number, $decimal ) = @_ ;
    my $round_num = 0 ;
    
	if ( ( defined $decimal ) and ( $decimal > 0 ) and ( defined $number ) and ( $number > 0 ) ) {
        $round_num = sprintf("%.".$decimal."f", $number);	## a rounding is used : 5.3 -> 5 and 5.5 -> 6
	}
	else {
		die "Can't round any number : missing value or decimal\n" ;
	}
    
    return(\$round_num) ;
}





sub get_html_link{
	my @param = @_;
	my $html = shift @param;
	my $title = shift @param;
	my @files = @param;
	open HTML, ">$html";
	print HTML "<html><head><title>$title</title></head><body><h3>Output Files:</h3><p><ul>\n";
	for (my $i = 0; $i <= $#files; $i++){
		my $name = basename $files[$i];
	    print HTML "<li><a href=\"$name\">$name</a></li>\n";
	}
	print HTML "</ul></p>\n";


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