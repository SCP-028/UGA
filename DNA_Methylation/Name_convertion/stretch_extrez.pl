#!/usr/bin/perl -w
use strict;
use warnings;
my $dir = "/home/yi/data/DNA_Methylation/Methy_ZhouYi/RData/t_test/";

open(HYPO,$dir."hypo.txt") or die "$!\n";
<HYPO>;
my @hypo_list = <HYPO>;
close(HYPO);

open(OUTPUT,">".$dir."hypo_new") or die "$!\n";
print OUTPUT "IlmnID"."\t"."Methy_p"."\t"."Entrez"."\t"."Refseq"."\t"."Gene group"."\t"."Relation to CpG"."\n";

foreach my $row(@hypo_list)
{
	my @a = split("\t", $row);
	#chomp $a[2];
	#chomp $a[3];
	#chomp $a[4];
	my @Entrez = split(";",$a[2]);
	my @Refseq = split(";",$a[3]);
	my @Gene_group = split(";",$a[4]);

	if($a[2] =~ m/;/)
	{
		my $length = scalar @Entrez;
		for(my $i=0;$i<$length;$i++)
		{
			print OUTPUT $a[0]."\t".$a[1]."\t".$Entrez[$i]."\t".$Refseq[$i]."\t".$Gene_group[$i]."\t".$a[5]."\n";
		}
	}
	else
	{
		print OUTPUT $a[0]."\t".$a[1]."\t".$a[2]."\t".$a[3]."\t".$a[4]."\t".$a[5];
	}
}
close(OUTPUT);