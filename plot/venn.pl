#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use File::Copy;
use FindBin qw($Bin);
my $r_path = "$Bin/../../bin/R/bin";
$ENV{'PATH'} = "$r_path;$ENV{'PATH'}";

my %opts;
GetOptions (\%opts,"f=s","l=s","dat=s","w=s","h=s","o=s", "html=s");

my $usage = <<"USAGE";
        Program : $0
        Discription: plot venn
        Usage:perl $0 [options]
                -f      files  a,b,c
                -l      labels x,y,z 
                -w      rank width de:9
                -h      rank height  de:9
                -dat    intersect data
                -o       out pdf
        example:
USAGE

die $usage if ( !( $opts{f} && $opts{l} && $opts{o} && $opts{dat} ) );


$opts{w}=$opts{w}?$opts{w}:9;
$opts{h}=$opts{h}?$opts{h}:9;


$opts{w} = $opts{w} / 2.54;
$opts{h} = $opts{h} / 2.54;

my $script_dir = "$Bin";

open RCMD, ">cmd.r";

########input files #####
my @files=split(/__FS__/,$opts{f});
my @labels=split(/__FS__/,$opts{l});

if( scalar(@files)!=scalar(@labels) ){
	die "number of files must be the same with labels\n";
}

print RCMD "
pdf(file=\"__out.pdf__\",width=$opts{w},height=$opts{h},pointsize=16)
lst <-list()";
for(my $i=0;$i<@files;$i++){
	my $copy_file = '__file__'.$i;
	copy( $files[$i], $copy_file );
	print RCMD "
lst[[\"$labels[$i]\"]]<-scan(\"$copy_file\",\"character\", sep=\"\\n\")
";
}


######## plot venn ######
        print RCMD "
	source(\"$script_dir/overLapper.new.r\")
	OLlist <- overLapper(setlist=lst, sep=\"\", type=\"vennsets\",keepdups=FALSE)
	#OLlist
	#names(OLlist)
	counts <- sapply(OLlist\$Venn_List, length)
	vennPlot(counts=counts,mymain=\"\")
	dev.off()
";

######## get intersect and unique ########
        print RCMD "
intersect_n <-function(x) {
        n <-length(x)
        x2 <- x[[1]]
        for(i in 2:n){
               x2 <- intersect(x[[i]],x2)
        } 
        return(x2)        
}

getset <- function(lst){   #####lst is a list
    sn <-length(names(lst))
    sets <- list()
    sls <- list()
    maxl <-0
    #### get all intersect ####  
    for(i in sn:2){
             sl <- combn (names(lst),i,simplify=FALSE)
             inter <-lapply(sl,function(x) intersect_n(lst[x]))
             names(inter) <-lapply(sl,function(x) paste(x,collapse =\" & \"))
             sets <- c(sets,inter)
             sls <- c(sls,sl)
                          
    }
    #### get all unique ####
    for(i in 1:sn){
             uniq <- list(setdiff(unlist(lst[names(lst)[i]]),unlist(lst[names(lst)[-i]])))
             names(uniq) <- paste(names(lst)[i],\" only\")
             sets <- c(sets,uniq)
             
    } 
    return(sets)

}
sets <- getset(lst)

###### write sets to file ######
maxl <-max(sapply(sets,length))
sets_table <- sapply(sets,function(x) c(x,rep(\"\",maxl-length(x))))

otuset <-\"__dat__\"
write.table(sets_table,otuset,sep = \"\t\",eol=\"\n\",row.names=FALSE,quote=FALSE)
";



close RCMD;
my $ret = system ('R --restore --no-save < cmd.r');
if ( $ret ) {
        die "Error, died with $ret";
}

move( '__out.pdf__', $opts{o} );
move( '__dat__', $opts{dat});
unlink 'cmd.r';



#&get_html_link($opts{html}, "Venn Analysis", $opts{dat}, $opts{o} );
get_html_link_v2($opts{html}, "Venn diagram Results", $opts{dat}, "tabH", $opts{o}, "pdf");








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
            }else{
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


















