#!/usr/bin/perl

$m = $ARGV[0]; #ms output file
$c = $ARGV[1]; #mean coverage 

open(FILE, "<$m") or die;

$head=<FILE>; chomp($head);
@spaces=split(/ /,$head);
$n=$spaces[1];
$s=$c*$spaces[2];

$l=0;
while (<FILE>) {
	$line = $_; chomp($line);
	if ($line eq "//") {
		$seg = <FILE>; chomp($seg); ($a, $b) = split(/: /,$seg);
		if ($b == 0) {
			$x=0;
			while ($x < $n) {
				$hap="0";
				$data[$l][$x]=$hap;
				$x++;
			}	
		} else {
			<FILE>; #skip position line
			$x=0;
			while ($x < $n) {
				$hap=<FILE>; chomp($hap);
				$data[$l][$x]=$hap;
				$x++;
			}
		}
		$l++;
	}
}
close FILE;

$x=1;
while ($x <= $s) {
	$sl=int(rand($l));
	$sn=int(rand($n));
	$hash[$sl]{$data[$sl][$sn]}++;
	$x++;
}

$x=0;
while ($x < $l) {
	$nhaps=keys(%{$hash[$x]});
	print "$nhaps\n";
	$x++;
}


