#!/usr/bin/perl -w

use warnings;
use strict;
use Getopt::Long;

my $communityFile;
my $cluFile;
my $nNode;

GetOptions ("c=s" => \$communityFile,
			"f=s" => \$cluFile,
			"n=i" => \$nNode);

die "[USAGE]>perl calculate_normalized_mutual_info.pl -c <communityFile> -f <cluFile> -n <numNode>\n"
unless defined($communityFile) and defined($cluFile) and defined($nNode);


my %origCommunity;		### key = communityID, value = member list.
my %clusterCommunity;	### key = moduleID, value = member list.


my @Community = `cut -f 2 $communityFile`;
my @tempModules = `cat $cluFile`;
chomp(@Community);
chomp(@tempModules);

my @Modules = @tempModules[1..$#tempModules];	## First line is '*Vertices #numbers'

## Check the number of elements of @Community and @Modules.
die "The Community array and Modules array are different in size.\n"
if (scalar(@Community) != scalar(@Modules));

for (my $i = 0; $i <= $#Community; $i++) {
	if (not exists($origCommunity{$Community[$i]})) {
		$origCommunity{$Community[$i]} = [$i];
	}
	else {
		push(@{$origCommunity{$Community[$i]}}, $i);
	}
}

for (my $i = 0; $i <= $#Modules; $i++) {
	if (not exists($clusterCommunity{$Modules[$i]})) {
		$clusterCommunity{$Modules[$i]} = [$i];
	}
	else {
		push(@{$clusterCommunity{$Modules[$i]}}, $i);
	}
}


my %probX;
my %probY;

foreach my $key (keys %origCommunity) {
	$probX{$key} = scalar(@{$origCommunity{$key}}) / $nNode;
}

foreach my $key (keys %clusterCommunity) {
	$probY{$key} = scalar(@{$clusterCommunity{$key}}) / $nNode;
}

my $mutualInfo = 0.0;

## calculate mutual information.
foreach my $keyX (keys %origCommunity) {
	foreach my $keyY (keys %clusterCommunity) {
		my $nXY = commonMembers($origCommunity{$keyX}, $clusterCommunity{$keyY});

		if ($nXY > 0) {
			my $probXY = $nXY / $nNode;
			#$mutualInfo += $probXY * log($probXY / ($probX{$keyX}*$probY{$keyY}));
			$mutualInfo += $probXY * (log($probXY) - log($probX{$keyX}) - log($probY{$keyY}));
		}
	}
}


## calculate H(X) and H(Y).
my $entropyX = 0.0;
my $entropyY = 0.0;

foreach my $key (keys %probX) {
	$entropyX -= $probX{$key} * log($probX{$key});
}

foreach my $key (keys %probY) {
	$entropyY -= $probY{$key} * log($probY{$key});
}

my $normalizedMI = 2 * $mutualInfo / ($entropyX + $entropyY);

print "$normalizedMI\n";





## Find common elements of two given arrays, 
## and return the number of the common elements.
sub commonMembers
{
	my ($arrRef1, $arrRef2) = @_;

	my @Array1 = sort {$a <=> $b} @$arrRef1;
	my @Array2 = sort {$a <=> $b} @$arrRef2;

	my $idx1 = 0;
	my $idx2 = 0;

	my $count = 0;
	while ($idx1 < scalar(@Array1) and $idx2 < scalar(@Array2)) {
		if ($Array1[$idx1] == $Array2[$idx2]) {
			$count++;
			$idx1++;
			$idx2++;
		}
		elsif($Array1[$idx1] > $Array2[$idx2]) {
			$idx2++;
		}
		else {
			$idx1++;
		}
	}

	return $count;
}
