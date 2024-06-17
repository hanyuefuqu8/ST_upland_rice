#!/usr/bin/perl -w

use strict;

my $file=shift;
my (%hasha,%hashb);
my %ar;
open IN,"$file" || die "Can't open such file:$!";
<IN>;
while (<IN>){
	chomp;
	my @a=split /\s+/,$_;
	$ar{$a[2]}{$a[0]}=$a[1];
     $hasha{$a[2]}=1;
	 $hashb{$a[0]}=1;
	}
close IN;
my @a=sort keys %hashb;
my $line=join(",",@a);
print ",$line\n";
foreach my $a(sort keys %hasha){
	my @temp;
	push @temp,$a;
	foreach (sort keys %hashb){
		if (!exists $ar{$a}{$_}){
			my $g=0;
			push @temp,$g;
			}else{
			push @temp,$ar{$a}{$_}	;
				}
		}
	 my $g2=join(",",@temp);
	 print "$g2\n";
	}
