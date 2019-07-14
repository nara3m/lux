use strict;
use warnings;

my ($sum, $least_sum, $centroid);

open A, "$ARGV[0]";

# Input file must look like this
 
# A,0
# B,0.9,0
# C,0.8,0.6,0
# D,0.7,0.5,0.4,0

my @input = <A>;
close A;
	
my $size = @input;

# Converting similairty matrix into 2-D array
my @temporary_array; my @sheet; my $row; my $col;

for($row = 0; $row < $size; $row++) 
{ 
	chomp ($input[$row]);
	@temporary_array = ();
	@temporary_array = split (/\,/, $input[$row]);
	
	$sum=0;

	for($col = 0; $col < $row+1; $col++) 
	{ 
		$sheet[$row][$col] = $temporary_array[$col+1];
		$sheet[$col][$row] = $sheet[$row][$col];
	}
}

$sum=0;

for($col = 0; $col < $size; $col++){ 	
		$sum = $sum+$sheet[0][$col];
}

$least_sum=$sum;
$centroid=0;

for($row = 1; $row < $size; $row++) 
{ 
	$sum =0;
	for($col = 0; $col < $size; $col++){ 	
		$sum = $sum+$sheet[$row][$col];
	}
	if ($sum < $least_sum) { $least_sum = $sum;$centroid=$row;}
}

# Print Result

#@temporary_array = ();
#@temporary_array = split (/\,/, $input[$centroid]);
#print "\nCentroid molecule is $temporary_array[0] \(Molecule number ",$centroid+1," in Input file\)\n\n";
print "$centroid";
__END__
