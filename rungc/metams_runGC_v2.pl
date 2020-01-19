#!/usr/bin/perl
use strict;
use warnings;
use File::Copy;
use File::Path;
use File::Basename;
use File::Copy::Recursive qw(dircopy);
use FindBin qw($Bin);
use Getopt::Long;
my $r_path = "$Bin/../../bin/R/bin";
$ENV{'PATH'} = "$r_path;$ENV{'PATH'}";


my %opts;
GetOptions (\%opts, "i=s", "peaktable=s", "peaktable2=s", "peakspectra=s", "html=s", "rtrange=s", "mzrange=s", 
"db_input=s", "fwhm=s", "rtdiff=s", "minfeat=s", "simthreshold=s", "minclassfraction=s", "minclasssize=s", "param=s",  "normpeakspectra=s", "nSlaves=s",
);

my $usage = <<"USAGE";
       Program : $0
       Version : 1.0
       Discription: 
       Usage :perl $0 [options]
                   -i  input cdf file list
                   -o  output html file
USAGE

die $usage if ( !$opts{i} );

sub get_term{
    my $file = shift;
    open IN,"$file";
    my $term = <IN>;
    chomp $term;
    close IN;
    return $term;
}




my $input = get_term($opts{i});
#my $peaktable = "pkTable_tmp.txt";
my $peakspectra = $opts{peakspectra};
my $html = $opts{html};
my $rtrange = $opts{rtrange};
my $mzrange = $opts{mzrange};
my $db_input = $opts{db_input};
my $fwhm = $opts{fwhm};
my $rtdiff = $opts{rtdiff};
my $minfeat = $opts{minfeat};
my $simthreshold = $opts{simthreshold};
my $minclassfraction = $opts{minclassfraction};
my $minclasssize = $opts{minclasssize};
my $param = $opts{param};
my $peaktable2 = $opts{peaktable2};
my $normpeakspectra = $opts{normpeakspectra};

my @cdf_files;

my @files = split /__SEP__/, $input;
foreach ( @files ) {
    my $basename = basename $_;

    unless( $basename=~/^[a-zA-Z][\w\.]+$/ ){
        die "File Name Error: \"$basename\", Valid name consists of letters, numbers and the dot or underline characters and starts with a letter not followed by a number or dot.\n";
    }

    #unlink $basename;
    mklink($_, $basename);
    push @cdf_files, $basename;
}

die "Error: input samples must be >= 2\n" if @cdf_files <2;

my $cdffiles = join "\",\"", @cdf_files;
$cdffiles = "c(\"". $cdffiles . "\")";

my $lib_path = "$Bin/runGC_lib.r";
copy( $lib_path, "./__libs.r__" );

mkdir "EICs";

open RCMD, ">cmd.r";
print RCMD <<__EORSCRIPT__;

library(metaMS)
data(FEMsettings)
source("./__libs.r__")
cdffiles=$cdffiles

__EORSCRIPT__

if ($param eq "default") {
    print RCMD "GALAXY.GC=TSQXLS.GC\n";
} else {
    print RCMD <<__EORSCRIPT__;
	GALAXY.GC <- metaMSsettings("protocolName" = "GALAXY.GC",
								"chrom" = "GC",
								PeakPicking = list(
								  method = "matchedFilter",
								  step = 0.5,
								  steps = 2,
								  mzdiff = .5,
								  fwhm = $fwhm,
								  snthresh = 2,
								  max = 500),
							   CAMERA = list(perfwhm = 1))
	metaSetting(GALAXY.GC, "DBconstruction") <- list(
				minintens = 0.0,
				rttol = $rtdiff,
				intensityMeasure = "maxo",
				DBthreshold = .80, 
				minfeat = $minfeat)
	metaSetting(GALAXY.GC, "match2DB") <- list(
				simthresh = $simthreshold,
				timeComparison = "rt",
				rtdiff = $rtdiff,
				RIdiff = 5,
				minfeat = $minfeat)
	metaSetting(GALAXY.GC, "betweenSamples") <- list(
				min.class.fraction = $minclassfraction,
				min.class.size = $minclasssize,
				timeComparison = "rt",
				rtdiff = $rtdiff,
				RIdiff = 2,    
				simthresh = $simthreshold)
__EORSCRIPT__
}

if (defined($rtrange)) {
    $rtrange = "c(" . $rtrange . ")";
} else {
    $rtrange = "NULL";
}
if (defined($mzrange)) {
    $mzrange = "c(" . $mzrange . ")";
} else {
    $mzrange = "NULL";
}
if (defined($db_input)) {
    mklink($db_input, "db.msp");
    print RCMD "manual <- read.msp(\"db.msp\")\n";
    print RCMD "DBarg <- createSTDdbGC(stdInfo = NULL, settings = GALAXY.GC, manualDB = manual)\n";
} else {
    print RCMD "DBarg = NULL\n";
}


print RCMD <<__EORSCRIPT__;

resGC <- runGC(files = cdffiles, settings=GALAXY.GC, rtrange=$rtrange, removeArtefacts = TRUE, findUnknowns = TRUE, DB= DBarg, nSlaves = $opts{nSlaves}, returnXset = TRUE)
peaktable = resGC\$PeakTable
spectra = resGC\$PseudoSpectra

getBPC2s(files = cdffiles, xset = NULL, rt="raw", pdfname="BPCs.pdf")
getTIC2s(files = cdffiles, xset = NULL, rt="raw", pdfname="TICs.pdf")
unknarg <- c(1:nrow(peaktable))
plotUnknowns(resGC=resGC, unkn=unknarg, DB=DBarg, fileFrom="singlefile")

N = dim(peaktable)[1]
for( i in 1:N){
    peaktable[[1]][i] = paste("MC", i, sep="")
    spectra[[i]]\$Name = paste("MC", i, sep="")
}

#order = order(peaktable[, "rt"])
#peaktable = peaktable[order,]

write.table(peaktable, file="peaktable.tsv", sep="\t", row.names=FALSE, quote=FALSE)
write.msp(spectra, file="peakspectra.msp", newFile = TRUE)

__EORSCRIPT__


close RCMD;
&run_cmd("R --restore --no-save<cmd.r");



get_clean_pktable("peaktable.tsv", $peaktable2);
copy("peakspectra.msp", $peakspectra);
norm_msp($peakspectra, $normpeakspectra);

#my $html_dir = dirname($html);
#dircopy("./EICs", "$html_dir/EICs") or die $!;

&get_html_link_v3($html, "GC-MS Peak Table Profiles",  $peaktable2, "tabH", $peakspectra, "txt", $normpeakspectra, "txt",
"TICs.pdf", "pdf", 
"BPCs.pdf", "pdf", 
"EICs", "dir",
);



sub get_clean_pktable{
    my $in = shift;
    my $out = shift;
    open IN,"$in";
    open OUT,">$out";
    my $h = <IN>;
    chomp $h;
    my @s = split /\t/, $h;
    my $index;
    for(my $i =1; $i <= $#s; $i++){
        if ($s[$i] eq "rt") {
            $index = $i + 1;
            last;
        }
    }
    my @l = @s[0,$index..$#s];
    my $l  = join "\t", @l;
    print OUT "$l\n";
    while (<IN>) {
        chomp;
        @s = split /\t/, $_;
        @l = @s[0,$index..$#s];
        $l  = join "\t", @l;
        print OUT "$l\n";
    }

    close IN;
    close OUT;
}


# in.table in.msp out.table out.msp
sub rename_pktable_msp{
    my $in_table = shift;
    my $in_msp = shift;
    my $out_table = shift;
    my $out_msp = shift;
    open IN, "$in_table";
    open OUT,">$out_table";
    my $h = <IN>;
    print OUT "$h";
    my $num = 1;
    my %hs;
    while(<IN>){
        chomp;
        my @s = split /\t/, $_;
        my $new_name = "MC" . $num;
        $hs{$s[0]} = $new_name;
        $s[0] = $new_name;
        my $line = join "\t", @s;
        print OUT "$line\n";
        $num++;
    }
    close IN;
    close OUT;

    open IN, "$in_msp";
    open OUT,">$out_msp";
    while(<IN>){
        if($_=~/Name:/){
            my $name = $_;
            $name=~s/Name://;
            $name = trim($name);
            $name = $hs{$name};
            print OUT "Name: $name\n";
        }else{
            print OUT "$_";
        }
    }

    close IN;
    close OUT;
}


sub sort_metams{
    my $in = shift;
    my $out = shift;
    open IN,"$in";
    open OUT,">$out";
    my $h = <IN>;
    print OUT "$h";
    my @aa;
    while(<IN>){
        my $name = (split /\t/, $_)[0];
        $name =~s/Unknown//;
        $name = trim($name);
        push @aa, [$name, $_];
    }
    my @sort = sort { $a->[0] <=> $b->[0] } @aa;
    for(my $i = 0; $i <= $#sort; $i++){
        my $line = $sort[$i][1];
        print OUT  "$line";
    }

    close IN;
    close OUT;
}




sub trim
{
        my $string = shift;
        $string =~ s/^\s+//;
        $string =~ s/\s+$//;
        return $string;
}



sub norm_msp{

my $peakspectra = shift;
my $normpeakspectra = shift;
open (MSP , "<" , $peakspectra) or die $! ;
open OUT, ">$normpeakspectra";
local $/ = 'Name:';
<MSP>;
while(<MSP>){
    chomp;
    my @infos = split /\n/ , $_;
    my @temp_mzs;
    my @temp_intensities;
    print OUT "Name:";
    for(my $i=0 ; $i<@infos ; $i++) {
        if ($infos[$i] =~ /(\d+\.?\d*)\s+(\d+\.?\d*)\s*;\s*/) {
            my @ions = split ( /;/ , $infos[$i] );
            foreach my $ion (@ions) {
                if ($ion =~ /^\s*(\d+\.?\d*)\s+(\d+\.?\d*)$/) {
                    my $mz = $1;
                    my $intensity =$2 ;
                    push (@temp_intensities , $intensity ) ;
                    push (@temp_mzs , $mz) ;
		        }
            }  	
        } else {
            unless ($infos[$i]=~/^\s*$/) {
                print OUT "$infos[$i]\n";
            }
        }
    }
    my $sum = 0;
    foreach  (@temp_intensities) {
        $sum += $_;
    }
    foreach  (@temp_intensities) {
        $_ = ($_ * 999)/$sum;
        $_ = sprintf("%.2f", $_) ;
    }
    print OUT " ";
    for(my $i = 0; $i<=$#temp_mzs; $i++){
        print OUT "$temp_mzs[$i] $temp_intensities[$i]; ";
    }
    print OUT "\n\n";
}
close MSP;
close OUT;

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

    close HTML;
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