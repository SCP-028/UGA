#!/usr/bin/perl -w
use strict;
use warnings;
my $dir = "/home/yi/data/DNA_Methylation/Methy_ZhouYi/RData/Glmnet/";
my $file1 = "annot_removed.txt";
my $file2 = "methyn";
my $file3 = "methyn_name";

open(ANNOT,$dir.$file1) or die "Could not open $file1: $!\n";
<ANNOT>;
my @annotation_list = <ANNOT>;
close(ANNOT);

open(METHY,$dir.$file2) or die "Could not open $file2: $!\n";
<METHY>;
my @methy_list = <METHY>;
close(METHY);

open(OUTPUT,">".$dir.$file3) or die "Could not open $file3: $!\n";
print OUTPUT "IlmnID"."\t"."Gene_name"."\t"."TCGA_AZ_6598"."\t"."TCGA_AA_3697"."\t"."TCGA_A6_2684"."\t"."TCGA_AZ_6601"."\t"."TCGA_A6_2680"."\t"."TCGA_AA_3712"."\t"."TCGA_A6_2671"."\t"."TCGA_A6_2682"."\t"."TCGA_AA_3655"."\t"."TCGA_A6_2685"."\t"."TCGA_AZ_6599"."\t"."TCGA_AA_3660"."\t"."TCGA_AA_3713"."\t"."TCGA_A6_2686"."\t"."TCGA_A6_2679"."\t"."TCGA_AA_3663"."\t"."TCGA_A6_5667"."\t"."TCGA_AZ_6600"."\t"."TCGA_A6_2675"."\n";
# iterate through the cases.
foreach my $row(@methy_list)
{
	my @a = split("\t",$row);
	foreach my $line(@annotation_list)
	{
		my @b = split("\t",$line);
		if($a[0] eq $b[0])
		{
			shift(@a);
			print OUTPUT $b[0].$b[1]."\t".join("\t",@a);
		}
	}
}
close(OUTPUT);