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
GetOptions (\%opts, "i=s",    "html=s", "rtrange=s", "mzrange=s", "settings=s",  "output1=s",  "peakpicking=s", "nSlaves=s", "output2=s", 
"ppm=s", "peakwidth=s", "prefilter=s", "step=s", "fwhm=s", "max=s", "snthresh=s", 
"minfrac=s", "minsamp=s", "mzwid=s", "bws=s", "missingratio=s", "extraratio=s", 
"retcor=s", "family=s", "profStep=s", "fillPeaks=s", 
"camera_perfwhm=s", "camera_cor_eic_th=s", "camera_ppm=s",
"intensity=s", "polarity=s",
);

my $usage = <<"USAGE";
       Program : $0
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
my $peaktable = $opts{output1};
my $cleantable =  $opts{output2};
my $html = $opts{html};


my $rtrange = $opts{rtrange};
my $mzrange = $opts{mzrange};
my $settings = $opts{settings};


my @cdf_files;

my @files = split /__SEP__/, $input;
foreach ( @files ) {
    my $basename = basename $_;

    unless( $basename=~/^[a-zA-Z][\w\.]+$/ ){
        die "File Name Error: \"$basename\", Valid name consists of letters, numbers and the dot or underline characters and starts with a letter not followed by a number or dot.\n";
    }

    mklink($_, $basename);
    push @cdf_files, $basename;
}

die "Error: input samples must be >= 2\n" if @cdf_files <2;

my $cdffiles = join "\",\"", @cdf_files;
$cdffiles = "c(\"". $cdffiles . "\")";


my $lib_path = "$Bin/runLC_lib.r";
copy( $lib_path, "./__libs.r__" );

mkdir "EICs";

open RCMD, ">cmd.r";
print RCMD <<__EORSCRIPT__;

library(metaMS)
data(FEMsettings)
source("./__libs.r__")
cdffiles=$cdffiles

__EORSCRIPT__




if( $settings eq "user" ){

    print RCMD "IP4M.LC = metaMSsettings( protocolName = \"USER\", chrom = \"LC\" )\n";

    if( $opts{peakpicking} eq "matchedFilter" ){
        print RCMD "metaSetting(IP4M.LC, \"PeakPicking\") = list( method = \"matchedFilter\",  step = $opts{step}, fwhm = $opts{fwhm},  snthresh = $opts{snthresh},    max = $opts{max}   )\n";
    }else{

         print RCMD "metaSetting(IP4M.LC , \"PeakPicking\") = list(  method = \"centWave\",   ppm = $opts{ppm},    prefilter = c($opts{prefilter}),    peakwidth = c($opts{peakwidth})   )\n";

    }


    if( $opts{retcor} eq "linear"){
        print RCMD "metaSetting(IP4M.LC, \"Alignment\") = list(    min.class.fraction = $opts{minfrac},   min.class.size = $opts{minsamp},  mzwid = $opts{mzwid},  bws = c($opts{bws}),   missingratio = $opts{missingratio},   extraratio =$opts{extraratio},  Retcor = list(  method = \"linear\", plottype=\"deviation\" , family = \"$opts{family}\"),    fillPeaks = $opts{fillPeaks}   )\n";
    }else{
        print RCMD "metaSetting(IP4M.LC, \"Alignment\") = list(    min.class.fraction = $opts{minfrac},   min.class.size = $opts{minsamp},  mzwid = $opts{mzwid},  bws = c($opts{bws}),   missingratio = $opts{missingratio},   extraratio =$opts{extraratio},  Retcor = list(  method = \"obiwarp\", plottype=\"deviation\" , profStep = $opts{profStep} ),    fillPeaks = $opts{fillPeaks}   )\n";
    }

    print RCMD  "metaSetting(IP4M.LC, \"CAMERA\") = list( perfwhm = $opts{camera_perfwhm}, cor_eic_th = $opts{camera_cor_eic_th},  ppm = $opts{camera_ppm}  )\n";

}else{
    print RCMD "IP4M.LC=$settings\n";
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


print RCMD <<__EORSCRIPT__;


pdf("rt_deviation_plot.pdf", width=16, height=10)
resLC<- runLC(files = cdffiles, settings=IP4M.LC, rtrange=$rtrange, mzrange=$mzrange, nSlaves = $opts{nSlaves},  returnXset=TRUE, intensity="$opts{intensity}", polarity="$opts{polarity}")
dev.off()

#peaktable<-resLC\$PeakTable<-resLC\$PeakTable[order(resLC\$PeakTable[,"rt"]),]
#write.table(peaktable, file="peaktable.tsv", sep="\t", row.names=FALSE, quote=FALSE)

xsAnnotate = resLC\$xset
xcmsSet=xsAnnotate\@xcmsSet
getTICs(xcmsSet=xcmsSet,pdfname="raw_tics.pdf",rt="raw")
getBPCs(xcmsSet=xcmsSet,pdfname="raw_bpcs.pdf",rt="raw")
getTICs2(xcmsSet=xcmsSet, pdfname="rtcorrected_tics.pdf",rt="corrected")
getBPCs2(xcmsSet=xcmsSet,pdfname="rtcorrected_bpcs.pdf",rt="corrected")
#getDeviation(xcmsSet=xcmsSet,pdfname="RT_Deviation_Plot2.pdf")

pt <- getPeaklist(xsAnnotate, intval="$opts{intensity}")
old_pt = pt
sortorder <- order(pt[,"rt"])
pt <- pt[sortorder,]
dataMatrix <- pt[,make.names(sampnames(xsAnnotate\@xcmsSet))] 
peakTable <- data.frame("pcgroup" = as.numeric(pt\$pcgroup),
                            "adduct" = pt\$adduct,
                            "isotopes" = pt\$isotopes,
                            "mz" = pt\$mz,
                            "rt" = pt\$rt/60,
                            dataMatrix,
                            row.names = NULL,stringsAsFactors = FALSE)
write.table(peakTable, file="peaktable.tsv", sep="\t", row.names=FALSE, quote=FALSE)
write.table(pt, file="peaktable_v2.tsv", sep="\t", row.names=FALSE, quote=FALSE)
getEICs_v2(xcmsSet=xcmsSet, pt=old_pt, sortorder=sortorder );

__EORSCRIPT__


close RCMD;
&run_cmd("R --restore --no-save<cmd.r");




deal("peaktable.tsv", $peaktable);
get_clean_table($peaktable, $cleantable );


&get_html_link_v3($html, "LC-MS Peak Table Profiles", $peaktable, "tabH", "$cleantable", "tabH", "raw_tics.pdf", "pdf", 
"raw_bpcs.pdf", "pdf",
"rtcorrected_tics.pdf", "pdf",
"rtcorrected_bpcs.pdf", "pdf",
"rt_deviation_plot.pdf", "pdf",
"EICs", "dir",
);



sub get_clean_table {
    my $in = shift;
    my $out = shift; 
	open IN,"$in";
	open OUT,">$out";
    my $h = <IN>;
    chomp $h;
    my @s = split /\t/, $h;
    my $line = join "\t", @s[0, 2..$#s];
    print OUT "$line\n";
    while(<IN>){
        chomp;
        my @s = split /\t/, $_;
        my $line = join "\t", @s[0, 2..$#s];
        if ( $s[1]=~/^\s*$/) {
            print OUT "$line\n";
        }elsif( $s[1]=~/\[M\]/ ){
            print OUT "$line\n";
        }
    }


    close IN;
    close OUT;
}



sub deal {
    my $in = shift;
    my $out = shift;
	open IN,"$in";
	open OUT,">$out";
    my $h = <IN>;
    chomp $h;
    my @s = split /\t/, $h;
    shift @s;
    shift @s;
    #shift @s;
    unshift @s, "ID";
    my $line = join "\t", @s;
    print OUT "$line\n";
    my $num = 1;
    while(<IN>){
        chomp;
        my @s = split /\t/, $_;
        shift @s;
        shift @s;
        #shift @s;
        unshift @s, "LC" . $num;
        my $line = join "\t", @s;
        print OUT "$line\n";
        $num ++;
    }

    close IN;
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
    close IN;
    close OUT;
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


