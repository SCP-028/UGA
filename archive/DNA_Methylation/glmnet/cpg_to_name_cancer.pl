#!/usr/bin/perl -w
use strict;
use warnings;
my $dir = "/home/yz73026/Glmnet/";
my $file1 = "annot_removed.txt";
my $file2 = "methyt";
my $file3 = "methyt_name";

open(ANNOT,$dir.$file1) or die "Could not open $file1: $!\n";
<ANNOT>;
my @annotation_list = <ANNOT>;
close(ANNOT);

open(METHY,$dir.$file2) or die "Could not open $file2: $!\n";
<METHY>;
my @methy_list = <METHY>;
close(METHY);

open(OUTPUT,">".$dir.$file3) or die "Could not open $file3: $!\n";
print OUTPUT "IlmnID"."\t"."Gene_name"."\t"."TCGA_D5_6530"."\t"."TCGA_G4_6323"."\t"."TCGA_A6_6781"."\t"."TCGA_AZ_4616"."\t"."TCGA_CK_6746"."\t"."TCGA_AZ_6607"."\t"."TCGA_CM_6677"."\t"."TCGA_AD_6888"."\t"."TCGA_D5_6929"."\t"."TCGA_AD_6895"."\t"."TCGA_D5_6537"."\t"."TCGA_G4_6309"."\t"."TCGA_AZ_6599"."\t"."TCGA_CM_5868"."\t"."TCGA_CM_4743"."\t"."TCGA_AZ_4682"."\t"."TCGA_CM_5861"."\t"."TCGA_CK_6748"."\t"."TCGA_CM_6679"."\t"."TCGA_AA_3662"."\t"."TCGA_CK_5914"."\t"."TCGA_A6_2684"."\t"."TCGA_A6_4105"."\t"."TCGA_AU_3779"."\t"."TCGA_G4_6297"."\t"."TCGA_G4_6299"."\t"."TCGA_D5_7000"."\t"."TCGA_A6_6138"."\t"."TCGA_A6_6782"."\t"."TCGA_CK_5913"."\t"."TCGA_G4_6317"."\t"."TCGA_D5_6539"."\t"."TCGA_CM_5348"."\t"."TCGA_G4_6302"."\t"."TCGA_F4_6460"."\t"."TCGA_AM_5821"."\t"."TCGA_F4_6855"."\t"."TCGA_A6_2684"."\t"."TCGA_AY_6196"."\t"."TCGA_A6_6648"."\t"."TCGA_A6_5660"."\t"."TCGA_CM_6165"."\t"."TCGA_D5_6930"."\t"."TCGA_A6_6142"."\t"."TCGA_F4_6570"."\t"."TCGA_A6_2684"."\t"."TCGA_A6_6652"."\t"."TCGA_CM_6162"."\t"."TCGA_CA_6717"."\t"."TCGA_A6_6780"."\t"."TCGA_A6_6781"."\t"."TCGA_CA_5796"."\t"."TCGA_A6_6650"."\t"."TCGA_D5_5539"."\t"."TCGA_D5_6924"."\t"."TCGA_AA_3506"."\t"."TCGA_G4_6311"."\t"."TCGA_CA_6719"."\t"."TCGA_CM_6171"."\t"."TCGA_AZ_4315"."\t"."TCGA_D5_6923"."\t"."TCGA_A6_5659"."\t"."TCGA_G4_6626"."\t"."TCGA_AZ_4614"."\t"."TCGA_AA_3496"."\t"."TCGA_CK_6751"."\t"."TCGA_CK_4948"."\t"."TCGA_AU_6004"."\t"."TCGA_D5_6529"."\t"."TCGA_CM_5863"."\t"."TCGA_AD_6548"."\t"."TCGA_F4_6807"."\t"."TCGA_CM_5862"."\t"."TCGA_AZ_6600"."\t"."TCGA_G4_6293"."\t"."TCGA_CA_6715"."\t"."TCGA_AY_5543"."\t"."TCGA_AA_3697"."\t"."TCGA_AD_6963"."\t"."TCGA_G4_6294"."\t"."TCGA_AZ_6603"."\t"."TCGA_AA_3712"."\t"."TCGA_A6_6654"."\t"."TCGA_G4_6320"."\t"."TCGA_D5_6533"."\t"."TCGA_D5_5540"."\t"."TCGA_CK_4952"."\t"."TCGA_F4_6809"."\t"."TCGA_CM_6674"."\t"."TCGA_AD_6964"."\t"."TCGA_CM_6172"."\t"."TCGA_AY_6386"."\t"."TCGA_AZ_5407"."\t"."TCGA_G4_6315"."\t"."TCGA_G4_6306"."\t"."TCGA_D5_6920"."\t"."TCGA_AA_3502"."\t"."TCGA_CM_6168"."\t"."TCGA_G4_6628"."\t"."TCGA_F4_6463"."\t"."TCGA_A6_5659"."\t"."TCGA_A6_6141"."\t"."TCGA_D5_6540"."\t"."TCGA_CM_6161"."\t"."TCGA_F4_6569"."\t"."TCGA_F4_6856"."\t"."TCGA_CM_6166"."\t"."TCGA_CM_6170"."\t"."TCGA_CA_6718"."\t"."TCGA_CK_4951"."\t"."TCGA_A6_5661"."\t"."TCGA_G4_6627"."\t"."TCGA_G4_6317"."\t"."TCGA_D5_6922"."\t"."TCGA_A6_5656"."\t"."TCGA_D5_6898"."\t"."TCGA_AZ_4313"."\t"."TCGA_AD_6901"."\t"."TCGA_CA_5797"."\t"."TCGA_G4_6310"."\t"."TCGA_G4_6303"."\t"."TCGA_D5_5538"."\t"."TCGA_A6_2682"."\t"."TCGA_D5_6931"."\t"."TCGA_A6_6780"."\t"."TCGA_A6_6653"."\t"."TCGA_CM_6163"."\t"."TCGA_CA_6716"."\t"."TCGA_A6_5666"."\t"."TCGA_CA_5256"."\t"."TCGA_CM_5349"."\t"."TCGA_AA_3509"."\t"."TCGA_A6_5656"."\t"."TCGA_AM_5820"."\t"."TCGA_F4_6461"."\t"."TCGA_F4_6854"."\t"."TCGA_A6_2685"."\t"."TCGA_CM_4744"."\t"."TCGA_AY_6197"."\t"."TCGA_CM_6164"."\t"."TCGA_A6_6649"."\t"."TCGA_A6_5661"."\t"."TCGA_A6_2675"."\t"."TCGA_G4_6304"."\t"."TCGA_AA_3489"."\t"."TCGA_A6_6650"."\t"."TCGA_G4_6588"."\t"."TCGA_AZ_5403"."\t"."TCGA_D5_6538"."\t"."TCGA_CK_5912"."\t"."TCGA_CM_6680"."\t"."TCGA_AZ_6608"."\t"."TCGA_CM_5860"."\t"."TCGA_F4_6805"."\t"."TCGA_AA_3663"."\t"."TCGA_CM_6678"."\t"."TCGA_CK_5915"."\t"."TCGA_G4_6586"."\t"."TCGA_AD_6889"."\t"."TCGA_A6_6137"."\t"."TCGA_A6_6651"."\t"."TCGA_A6_6781"."\t"."TCGA_D5_6536"."\t"."TCGA_G4_6298"."\t"."TCGA_AA_3492"."\t"."TCGA_AZ_6598"."\t"."TCGA_G4_6322"."\t"."TCGA_D5_6531"."\t"."TCGA_CK_6747"."\t"."TCGA_CK_4950"."\t"."TCGA_AZ_6606"."\t"."TCGA_AA_3495"."\t"."TCGA_CM_6676"."\t"."TCGA_A6_5664"."\t"."TCGA_A6_2686"."\t"."TCGA_D5_6928"."\t"."TCGA_CM_4747"."\t"."TCGA_F4_6704"."\t"."TCGA_A6_5662"."\t"."TCGA_CM_6167"."\t"."TCGA_A6_6140"."\t"."TCGA_D5_6932"."\t"."TCGA_G4_6307"."\t"."TCGA_D5_6541"."\t"."TCGA_A6_5667"."\t"."TCGA_A6_5665"."\t"."TCGA_A6_6650"."\t"."TCGA_F4_6703"."\t"."TCGA_CM_6169"."\t"."TCGA_CM_5344"."\t"."TCGA_AA_3511"."\t"."TCGA_AZ_4684"."\t"."TCGA_D5_6926"."\t"."TCGA_AD_5900"."\t"."TCGA_AA_3655"."\t"."TCGA_G4_6314"."\t"."TCGA_CA_5254"."\t"."TCGA_D5_6532"."\t"."TCGA_G4_6321"."\t"."TCGA_D5_6534"."\t"."TCGA_AD_6890"."\t"."TCGA_D5_5541"."\t"."TCGA_AZ_6605"."\t"."TCGA_F4_6808"."\t"."TCGA_CA_5255"."\t"."TCGA_AD_6965"."\t"."TCGA_G4_6625"."\t"."TCGA_CM_6675"."\t"."TCGA_CK_4947"."\t"."TCGA_F4_6459"."\t"."TCGA_D5_6535"."\t"."TCGA_D5_6927"."\t"."TCGA_AA_3713"."\t"."TCGA_A6_5659"."\t"."TCGA_F4_6806"."\t"."TCGA_AZ_6601"."\t"."TCGA_AA_3660"."\t"."TCGA_CK_5916"."\t"."TCGA_A6_5656"."\t"."TCGA_D5_5537"."\t"."TCGA_G4_6295"."\t"."TCGA_AZ_4615"."\t"."TCGA_A6_5657"."\t"."TCGA_CM_5864"."\t"."TCGA_AZ_4323"."\t"."TCGA_CM_4751"."\t"."TCGA_AD_6899"."\t"."TCGA_A6_6780"."\t"."TCGA_A6_5665"."\n";
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