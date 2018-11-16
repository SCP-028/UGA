#!/usr/bin/perl -w 
use strict;
my $dir = "/home/yi/data/DNA_Methylation/Methy_ZhouYi/";

open(INPUT, $dir."Accession_Num.csv") or die "$!\n";
<INPUT>;
my @Accession_list = <INPUT>;
close(INPUT);

open(HYPER, $dir."hyper_7.csv") or die "$!\n";
<HYPER>;
my @cg_list = <HYPER>;
close(HYPER);

open(OUTPUT, ">".$dir."identifier_7") or die "$!\n";
print OUTPUT "IlmnID"."\t"."RefGene_Name"."\t"."Fold"."\n";

foreach my $line(@cg_list)
{
	my @cg = split(",",$line);
	chomp $cg[0];
	#chomp $cg[2];
	foreach my $row(@Accession_list)
	{
		my @ref = split(",",$row);
		#chomp $ref[0];
		chomp $ref[1];
		if($ref[0] =~ m/^($cg[0])/)
		{
			print $ref[1];
			#print OUTPUT $cg[0]."\t".$ref[1]."\n";
		}
	}
}
close(OUTPUT);
