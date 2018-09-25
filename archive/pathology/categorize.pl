#!/usr/bin/perl -w
use strict;
my $dir = "E:/bioInfo/Perl/";
chdir($dir."Pathology/data"); # change the directory.
my @final_list = grep(/^(TCGA)-(\w){2}-(\w){4}.*(\.pdf)$/, `ls -R`); # Only get TCGA***.pdf files

open(INPUT, $dir."case_id_project") or die "$!\n";
<INPUT>;
my @file = <INPUT>;
close(INPUT); # Store names in @file.

foreach my $pdf_id(@final_list)
{
	chomp $pdf_id;
	foreach my $line(@file)
	{
		my @a = split("\t", $line);
		chomp $a[1];
		$a[1] =~ s/[\n\r]*//;
		chomp $a[0];
		if($pdf_id =~ m/^($a[1])/)
		{
			my $path_old = `find $dir"Pathology/data" -name \"$pdf_id\"`;# Set dir to data, mind . and /
			chomp $path_old;
			my $path_new = $dir.$a[0];
			system("mv $path_old $path_new");
		}
	}
}