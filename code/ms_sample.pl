#!/usr/bin/perl

$m = $ARGV[0]; #ms output file
$c = $ARGV[1]; #mean coverage 
$e = $ARGV[2]; #error rate 0-1, 0.1 for 10%
$t = $ARGV[3]; #threshold to count haplotypes, 2 if haplotype needs to be seen twice
$ecount=0; #error count

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
	$seq=$data[$sl][$sn];
	$r=rand(1);
	if ($r <= $e) { $seq = $seq . "_" . $ecount; $ecount++; }
	$hash[$sl]{$seq}++;
	$x++;
}

$x=0;
while ($x < $l) {
	$nhaps=0;
	foreach $key ( keys(%{$hash[$x]}) ) {
		if ($hash[$x]{$key} >= $t) { $nhaps++; } 
	}
	print "$nhaps\n";
	$x++;
}


