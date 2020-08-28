#!/usr/bin/perl

$m = $ARGV[0]; #ms output file
$c = $ARGV[1]; #mean coverage 
$e = $ARGV[2]; #error rate 0-1, 0.1 for 10%
$t = $ARGV[3]; #threshold to count haplotypes, 2 if haplotype needs to be seen twice
$d = $ARGV[4]; #individual distribution file
$ecount=0; #error count


open(FILE, "<$d") or die;
$line=<FILE>; close(FILE); chomp($line);
@commas=split(/,/,$line);
$sum=0;
$sum2=0;
foreach (@commas) {
	$sum+=$_;
}
foreach (@commas) {
	$sum2+=$_;
	push(@dist,($sum2/$sum)); #make new normalized distribution array
}


open(FILE, "<$m") or die;

$head=<FILE>; chomp($head);
@spaces=split(/ /,$head);
$n=$spaces[1];
if ($n != (($#dist+1)*2)) {
	die; #number of haplotypes is not two times number of individuals in dist file
}
$s=$c*$spaces[2]; #s is number of sequences
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
	$sl=int(rand($l)); #choose random locus

	$si = rand(1);
	$y = 0;
	while ($y <= $#dist) {
		if ($si < $dist[$y]) {
			$si = $y; #sample from specific individual
			$y = $#dist+1; #exit loop
		}
		$y++;
	}	
	$sh=int(rand(2)); #first or second haplotype from that individual (will equal 0 or 1)
	$sn = ($si*2)+$sh; 
	
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


