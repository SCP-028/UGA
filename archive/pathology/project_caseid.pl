#!/usr/bin/perl -w  #be sure to include this line for every per script
use strict; #to make perl stricter on your code, so it will be safer.
my $dir = "/Users/mikeaalv/Documents/TCGA";
chdir($dir."/data/"); # change the directory.
my @project_list = `ls -d TCGA*`; #`` is for system command like ls, you could also use system(). These two are different as the former will return the result, try to baidu or google it
open(OUTPUT, ">".$dir."/pathology/case_id_project") or die"$!\n"; #open file or quit the program. > is for output. string are connected by ".". OUTPUT is called file handle.
foreach my $project_id (@project_list)  #recursion, @ is a array $ is a normal variable. this means work for every element in @project_list
{
  chomp $project_id;  # get rid of \n or \r in the end
  open(INPUT,$dir."/data/".$project_id."/".$project_id."_barcode_uuid") or die "$!\n";
  <INPUT>;  #this give you one line in the file
  my @file = <INPUT>; #this give every remaining line of the file to the array. a normal variable will get one line in this situaion. try to understand "push", "pop", "shift", "unshift"
  close(INPUT); #don't forget to close the file handle after using it.
  my %case_id;  #it is a hash. baidu or goole to understand it. it is an important concept. try to understand "keys" and "exists".
  foreach my $line (@file)
  {
    my @a = split("\t", $line); #baidu or google to understand "split" and "join". The former seperate a string into an array, while the latter join an array into a string.
    my $thre_sep;
    if($a[0] =~ m/^(?<mat>(TCGA-[\w]{2}-[\w]{4}))-/) # this is an example of regex, but there is more about regex. try to baid perl 正则 and learn it
    {
      $thre_sep = $+{mat}; #this means the matching part.
    }
    else
    {
      die "$!\n";
    }
    if(!(exists $case_id{$thre_sep})) #! is negate, this means if in the hash case_id, there isn't a keys==$a[0], the progress will go in the {}
    {
      $case_id{$thre_sep} = 1;
    }
  }
  foreach my $case (keys %case_id)
  {
    print OUTPUT $project_id."\t".$case."\n";  #OUTPUT
  }
}
close(OUTPUT);
#also remember to understand usage of hash, array, and try to understand "上下文" in perl
#also be sure to understand regex
