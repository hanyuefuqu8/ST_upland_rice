#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use PerlIO::gzip;

my $GeneName;
my $size;
my ($Match,$nMatch);
GetOptions(
	"GN:s"	=>	\$GeneName,
	"size:i"=>	\$size,

	"Match:s"=>	\$Match,
	"nMatch:s"=>	\$nMatch,
);
$size ||= 50;

@ARGV || die "Usage:perl $0 <file> <size>\n";

$size || die "Check the input size\n";

my %hash;

open IN,"$ARGV[0]" || die "Can't open such file:$!";
my $head=<IN>;
close IN;

#format1
#x       y       MIDCounts       geneID
#
#format2
#geneID  x       y       MIDCounts

my $format;

if($head=~/^x\ty\tMIDCounts\tgeneID/){
	$format = 1;
}elsif($head=~/^geneID\tx\ty\tMIDCounts/){
	$format = 2;
	$head = "x\ty\tMIDCounts\tgeneID\n";
}else{
	die "Check the format...\n";
}

print "$head";

foreach my $file(@ARGV){
($file && -e $file && -s $file) || die "Check the input file\n";
if($file =~ /\.gz$/){
	open IN,"<::gzip",$file;
}else{
	open IN,$file;
}
while (<IN>){
	/MIDCounts/ && next;

	chomp;
	my @a=split /\s+/,$_;

	my ($x,$y,$MIDCounts,$geneID);

	if($format==1){
		($x,$y,$MIDCounts,$geneID) = @a[0,1,2,3];
	}elsif($format==2){
		($x,$y,$MIDCounts,$geneID) = @a[1,2,3,0];
	}

	if($Match){
		$geneID =~ /$Match/ || next;
	}

	if($nMatch){
		$geneID =~ /$nMatch/ && next;
	}

#	if($format==1){
#		$x=int($a[0]/$size);
#		$y=int($a[1]/$size);
#		$hash{$x}{$y}{$a[3]} +=$a[2];
#	}elsif($format==2){
#		$x=int($a[1]/$size);
#		$y=int($a[2]/$size);
#		$hash{$x}{$y}{$a[0]} +=$a[3];
#	}

	my $X=int($x/$size);
	my $Y=int($y/$size);
	$hash{$X}{$Y}{$geneID} +=$MIDCounts;

}
close IN;
}

if($GeneName){
	foreach my $hen(sort {$a<=>$b} keys %hash){
		foreach my $shu(sort {$a<=>$b} keys %{$hash{$hen}}){
			my $geneID = $GeneName;
			$hash{$hen}{$shu}{$geneID} ||= 0;
			print "$hen\t$shu\t$hash{$hen}{$shu}{$geneID}\t$geneID\n";
		}
	}
}else{
	foreach my $hen(sort {$a<=>$b} keys %hash){
		foreach my $shu(sort {$a<=>$b} keys %{$hash{$hen}}){
			foreach my $geneID(keys %{$hash{$hen}{$shu}}){
				print "$hen\t$shu\t$hash{$hen}{$shu}{$geneID}\t$geneID\n";
			}
		}
	}
}
