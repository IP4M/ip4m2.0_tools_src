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
GetOptions (\%opts, "i=s", "analysis_time=s", "min_peak_width=s", "min_peak_height=s", "noise_threshold=s", "avoid_processing_mz=s", "min_spectra_cor=s", "max_time_dist=s", "mz_range=s", "min_samples=s", "blocks_size=s", "peaktable=s", "peakspectra=s", "html=s", "normpeakspectra=s",
);

my $usage = <<"USAGE";
    Error!
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



my $peaktable = $opts{peaktable};
my $peakspectra = $opts{peakspectra};
my $html = $opts{html};
my $normpeakspectra = $opts{normpeakspectra};

my $input = get_term($opts{i});
my @cdf_files;
my @ids;
my @files = split /__SEP__/, $input;
foreach ( @files ) {
    my $basename = basename $_;

    unless( $basename=~/^[a-zA-Z][\w\.]+$/ ){
        die "File Name Error: \"$basename\", Valid name consists of letters, numbers and the dot or underline characters and starts with a letter not followed by a number or dot.\n";
    }

    mklink($_, $basename);
    unless($basename=~/mzML$/i){
        push @cdf_files, $basename;
    }else{
        my $output = ml2xml($basename);
        push @cdf_files, $output;
    }
    my @aa =split /\./, $basename;
    pop @aa;
    my $id = join ".", @aa;
    push @ids, $id;
}

die "Error: input samples must be >= 2\n" if @cdf_files <2;

my $cdffiles = join "\",\"", @cdf_files;
$cdffiles = "c(\"". $cdffiles . "\")";

open OUT,">inst.csv";
open PHE,">phe.csv";
print OUT "sampleID;filename;date;time\n";
print PHE "sampleID;class\n";
for(my $i = 0; $i <= $#ids; $i++){
    print OUT "$ids[$i];$cdf_files[$i];;\n";
    print PHE "$ids[$i];\n";
}
close OUT;


my $lib_path = "$Bin/erah_fix.r";
copy( $lib_path, "./__libs.r__" );

mkdir "EICs";

open RCMD, ">cmd.r";
print RCMD <<__EORSCRIPT__;

library(erah)
library(metaMS)

#fix the new version export2MSP bugs by using the old one!
source("./__libs.r__")

cdffiles=$cdffiles
getBPC2s(files = cdffiles, xset = NULL, rt="raw", pdfname="BPCs.pdf")
getTIC2s(files = cdffiles, xset = NULL, rt="raw", pdfname="TICs.pdf")

ex <- newExp(instrumental = "inst.csv", phenotype="phe.csv" )

ex.dec.par <- setDecPar(min.peak.width = $opts{min_peak_width}, min.peak.height=$opts{min_peak_height}, noise.threshold=$opts{noise_threshold},  avoid.processing.mz=$opts{avoid_processing_mz}, analysis.time=$opts{analysis_time})
ex <- deconvolveComp(ex, ex.dec.par )
ex.al.par <- setAlPar(min.spectra.cor = $opts{min_spectra_cor}, max.time.dist = $opts{max_time_dist}, mz.range = $opts{mz_range})
ex <- alignComp(ex, alParameters = ex.al.par, blocks.size=$opts{blocks_size})
ex <- recMissComp(ex, min.samples = $opts{min_samples})

align=alignList(ex)
write.table(align, file="peaktable.tsv", sep="\t", row.names=FALSE, quote=FALSE)
old_export2MSP(ex)

ids = align[[1]]
for( id in ids){
     plotSpectra = plotSpectra(ex, id, compare=F, return.spectra=T )
     max_index = NULL
     for (i in  length(plotSpectra):1){
         if( plotSpectra[i] > 0 ){
             max_index = i
             break
         }
     }
     pdf(width=26, height=10, file=paste("EICs/EC", id, ".pdf", sep="")) #A4 pdf						
     plotAlign(ex, id , per.class=F)
     par(mfrow = c(1, 1))
     plotSpectra(ex, id, compare=F, return.spectra=T, xlim=c(0, max_index + 20) )
     dev.off()
}

__EORSCRIPT__
close RCMD;
&run_cmd("R --restore --no-save<cmd.r");

my $pkTable = "peaktable.tsv";
my $msp_file = "ExportMSP/ExportedMSP.msp";


#move($pkTable, $peaktable);
#move($msp_file, $peakspectra);

## deal peak table
open IN,"$pkTable";
open OUT,">peaktable";
my @foundin;
my @rts;
my $h = <IN>;
chomp $h;
my @s = split /\t/, $h;
my $tmp = join "\t", @s[0, 4..$#s];
print OUT "$tmp\n";
while(<IN>){
    chomp;
    my @s = split /\t/, $_;
    my $tmp = join "\t", @s[0, 4..$#s];
    print OUT "$tmp\n";
    push @foundin, $s[3];
    push @rts, $s[2];
}
close IN;
close OUT;

## deal msp
open IN,"$msp_file";
open OUT,">peakspectra";
my $num = 1;
while(<IN>){
    if ($_=~/^Name/) {
        my $fi = $foundin[$num -1 ];
        my $rt = $rts[$num - 1];
        print OUT "Name: $num\n";
        print OUT "rt: $rt\n";
        print OUT "FoundIn: $fi\n";
        $num ++;
    }else {
        print OUT "$_";
    }
}


rename_pktable_msp( "peaktable", "peakspectra", $peaktable, $peakspectra );
#&get_html_link($html, "Peak Table and MS Profile", $peaktable,  $peakspectra);
norm_msp($peakspectra, $normpeakspectra);
&get_html_link_v3($html, "GC-MS Peak Table Profiles", $peaktable, "tabH", $peakspectra, "txt", $normpeakspectra, "txt",
"TICs.pdf", "pdf", 
"BPCs.pdf", "pdf", 
"EICs", "dir",
);




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





sub trim
{
        my $string = shift;
        $string =~ s/^\s+//;
        $string =~ s/\s+$//;
        return $string;
}


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
        my $new_name = "EC" . $num;
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


sub ml2xml{
    my $in = shift;
    my @aa =split /\./, $in;
    pop @aa;
    my $id = join ".", @aa;
    my $out = "$id.mzXML";
    open RCMD, ">cmd2.r";
print RCMD <<__EORSCRIPT__;
library(mzR)
file = "$in"
out_file = "$out"
ms_fl <- openMSfile(file)
pks <- spectra(ms_fl)
hdr <- header(ms_fl)
writeMSData(object = pks, file = out_file, outformat = "mzxml", header = hdr)

__EORSCRIPT__
    close RCMD;
    &run_cmd("R --restore --no-save<cmd2.r");
    return $out;

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