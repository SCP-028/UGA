#!/usr/bin/perl -w
use strict;
use warnings;
my $dir = "/home/yi/data/DNA_Methylation/Methy_ZhouYi/RData/Methy_GX/";

my $file1 = "hypo_up.txt";
my $file2 = "GXU.txt";
#hyper & up-regulate
open(HYPERUP,$dir.$file1) or die "Could not open $file1: $!\n";
<HYPERUP>;
my @hyper_up_list = <HYPERUP>;
close(HYPERUP);

open(INPUT,$dir.$file2) or die "No $file2: $!\n";
<INPUT>;
my @gxu_list = <INPUT>;
close(INPUT);

open(OUTPUT,">".$dir."result/hypo_GXup") or die "$!\n";
print OUTPUT "IlmnID"."\t"."Methy_P"."\t"."GX_P"."\t"."Entrez"."\t"."Refseq"."\t"."Gene_group"."\t"."Relation_to_CpG"."\n";

foreach my $row(@hyper_up_list)
{
	my @a = split("\t",$row);
	foreach my $line(@gxu_list)
	{
		my @b = split("\t",$line);
		if($a[2] == $b[0])
		{
			print OUTPUT $a[0]."\t".$a[1]."\t".$b[1]."\t".$a[2]."\t".$a[3]."\t".$a[4]."\t".$a[5];
		}
	}
}
close(OUTPUT);