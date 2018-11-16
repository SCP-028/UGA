#!/usr/bin/perl -w
use strict;
use warnings;
my $dir = "/home/yi/data/DNA_Methylation/Methy_ZhouYi/RData/Glmnet/";

my $file1 = "annot";
my $file2 = "annot_removed";

open(ANNOT,$dir.$file1) or die "Could not open $file1: $!\n";
<ANNOT>;
my @annotation_list = <ANNOT>;
close(ANNOT);

open(OUTPUT,">".$dir.$file2) or die "Could not open $file2: $!\n";
print OUTPUT "IlmnID"."\t"."Gene_Name"."\t"."Refgene_Accession"."\t"."Gene_Group"."\t"."CpG_Islands_Name"."\t"."Relation_to_CpG_Island"."\n";

foreach my $row(@annotation_list)
{
	my @a = split("\t",$row); ## annot cols
	my @b = split(";",$a[1]); ## Gene name
	my @c = split(";",$a[2]); ## Refgene
	my @d = split(";",$a[3]); ## Gene group
	my $length = scalar @b;

	print OUTPUT $a[0]."\t".$b[0]."\t".$c[0]."\t".$d[0]."\t".$a[4]."\t".$a[5]."\n";
	for(my $i=1;$i<$length;$i++)
	{
		if($b[$i] ne $b[$i-1])
		{
			print OUTPUT $a[0]."\t".$b[$i]."\t".$c[$i]."\t".$d[$i]."\t".$a[4]."\t".$a[5]."\n";
		}
		elsif($d[$i] ne $d[$i-1])
		{
			print OUTPUT $a[0]."\t".$b[$i]."\t".$c[$i]."\t".$d[$i]."\t".$a[4]."\t".$a[5]."\n";
		}
	}
}
close (OUTPUT);