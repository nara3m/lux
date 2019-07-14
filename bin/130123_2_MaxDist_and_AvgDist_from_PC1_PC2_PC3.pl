use strict;
use warnings;

# This Program Computes Hausdroff Distance (denoted as Hdist) between Two Sets. Since Hdist is assymetrical [(Set1 -> Set2) != (Set2 -> Set1)], 
# maximum of the two is choosen as Hdist. 

# Along with Hdist, Average Minimum Distance between sets is also Calculated. This is denoted as AvgHdist.

# Both Hdist and AvgHdist are calculated in two ways 1) WithOut Concentrations 2) With Concentrations. A total of 4 Scores are returned for each pair.

# Input for this program are 1) List of Files for which pair-wise Hdist needs to be computed 2) Combined PCA file 3) Folder containing List of Files.
# It is assumed that directory structure is preserved !

# Output is written to file Hdist

#_______________________________________________________________________________________________

my $PCAfilename = $ARGV[0]; 
my $FileNameList = $ARGV[1];
my $InputFolderName = $ARGV[2];

# Reading PCA file
my %XYcors; # Hash table for saving XY coordinates for all molecules
my @TempArray1;
my $MinDist;
	
open PCA, "./$InputFolderName/$PCAfilename" or die "\n\nUnable to open $PCAfilename file\n";
while (<PCA>){
	push @TempArray1, $_;
}
close PCA;
	
# Editing PCA file
shift(@TempArray1); # Removing header line from PCA file
foreach (@TempArray1){
	chomp($_);
	$_ =~ s/\"//g;
	@_ = split(" ", $_);
	$XYcors{"$_[0]"} = "$_[1] $_[2] $_[3]";
}

# Writing Results to Files
# open OUT1, "> ./$InputFolderName/WithConcMaxdist_from_PC1_PC2_PC3.csv";
# open OUT2, "> ./$InputFolderName/WithOutConcMaxdist_from_PC1_PC2_PC3.csv";
# open OUT3, "> ./$InputFolderName/WithConcAvgDist_from_PC1_PC2_PC3.csv";
 open OUT4, "> ./$InputFolderName/LUX_PC1_PC2_PC3.csv";

# Calculating Average distance and Maximum distance for all pairs of Input files
my @FileNameList = split ("XXXX", $FileNameList);
my $NumOfInputFiles = scalar(@FileNameList);
my ($Var1, $Var2, $File1, $File2);

# print OUT1 "$FileNameList[0]\t0\n";
# print OUT2 "$FileNameList[0]\t0\n";
# print OUT3 "$FileNameList[0]\t0\n";
print OUT4 "$FileNameList[0]\t0\n";

for ($Var1 = 1; $Var1 < $NumOfInputFiles; $Var1++){

	$File1 = $FileNameList[$Var1];

#	print OUT1 "$File1\t";
#	print OUT2 "$File1\t";
#	print OUT3 "$File1\t";
	print OUT4 "$File1\t";

	for ($Var2 = 0; $Var2 < $Var1; $Var2++){

		$File2 = $FileNameList[$Var2];
		
		my ($WOC_Median_1, $WOC_Mean_1, $WOC_Max_1, $WOC_IQR_1, $WC_Median_1, $WC_Mean_1, $WC_Max_1, $WC_IQR_1) = RunDist ($File1, $File2);
		my ($WOC_Median_2, $WOC_Mean_2, $WOC_Max_2, $WOC_IQR_2, $WC_Median_2, $WC_Mean_2, $WC_Max_2, $WC_IQR_2) = RunDist ($File2, $File1);
		
		my $WOC_Median = max($WOC_Median_1, $WOC_Median_2);
		my $WOC_Mean = max($WOC_Mean_1, $WOC_Mean_2);
		my $WOC_Max = max($WOC_Max_1, $WOC_Max_2);
		my $WOC_IQR = max($WOC_IQR_1, $WOC_IQR_2);

		my $WC_Median = max($WC_Median_1, $WC_Median_2);
		my $WC_Mean = max($WC_Mean_1, $WC_Mean_2);
		my $WC_Max = max($WC_Max_1, $WC_Max_2);
		my $WC_IQR = max($WC_IQR_1, $WC_IQR_2);

#		print OUT1 "$WC_Max\t";
#		print OUT2 "$WOC_Max\t";
#		print OUT3 "$WC_Mean\t";
		print OUT4 "$WOC_Mean\t";
		
	}

#	print OUT1 "0\n";
#	print OUT2 "0\n";
#	print OUT3 "0\n";
	print OUT4 "0\n";
	
}

# close OUT1; close OUT2; close OUT3; 
close OUT4;

#######
sub RunDist {

	my ($File1, $File2) = ($_[0], $_[1]);
	
	# Reading Sample1 file
	
	my %Sample1;
	my @TempArray2 = ();
	
	open Sample1, "./$InputFolderName/$File1" or die "\n\nUnable to open $File1 file\n";
	while (<Sample1>){
		chomp($_);
		@TempArray2 = split("\t", $_);
		$TempArray2[0] =~ s/\s/\_/g;
		$TempArray2[0] =~ s/\:/\_/g;
		$TempArray2[0] =~ s/\;/\_/g;
		#$TempArray2[0] =~ s/\//\_/g;
		
		$Sample1{$TempArray2[0]} = "$TempArray2[1]";
	}

	close Sample1;
	
	# Reading Sample 1 file
	
	my %Sample2;
	
	open Sample2, "./$InputFolderName/$File2" or die "\n\nUnable to open $File2 file\n";
	while (<Sample2>){
		chomp($_);
		@TempArray2 = split("\t", $_);
		$TempArray2[0] =~ s/\s/\_/g;
		$TempArray2[0] =~ s/\:/\_/g;
		$TempArray2[0] =~ s/\;/\_/g;
		#$TempArray2[0] =~ s/\//\_/g;
		
		$Sample2{$TempArray2[0]} = "$TempArray2[1]";
	}

	close Sample2;

	# Normalizing Concentrations by Maximum Concentration
	
	my @TempArray1 = ();
	push (@TempArray1, values(%Sample1), values(%Sample2));
	my $MaxConc = max(@TempArray1); # print "MaxConc = $MaxConc\n";

	foreach (keys (%Sample1)){ 
		$Sample1{$_} = ($Sample1{$_}/$MaxConc); # print "$Sample1{$_}\n";
	}
	foreach (keys (%Sample2)){
		$Sample2{$_} = ($Sample2{$_}/$MaxConc); # print "$Sample2{$_}\n";
	}
	
	# Iterating through all pairs
	
	my $WithOutConcDist;

	my $WithOutConcMinDist; # Minimum distance between this molecule (in Sample1) to every other molecule in Sample2 
	my $WithConcMinDist;
	my $MoleculeThatGaveMinDist;

	my @WithOutConcMinDist = (); # Array to store All min distances - will be used to gather Hausdroff distance and other metrics
	my @WithConcMinDist = ();

	my $WithConcHdist = 0; # Hausdroff Distance for Sample1 and Sample2
	my $WithOutConcHdist = 0;
	
	my $WithConcAvgMinDist = 0; # Average Minimum Distance between all pairs of molecules in Sample1 and Sample2
	my $WithOutConcAvgMinDist = 0;
	
	foreach (keys (%Sample1)){ # Keys are Molecule names
		
		my $p1Name = $_;
		if (defined($XYcors{$p1Name})){

			$WithOutConcMinDist = 100000000; # Minimum distance between this molecule (in Sample1) to every other molecule in Sample2 

			my ($xp1, $yp1, $zp1) = split (" ", $XYcors{$p1Name});
			my $cp1 = $Sample1{$p1Name}; # Concentration of Molecule
		
			EXIT1: {
		
			foreach (keys (%Sample2)){
		
				my $p2Name = $_;
				if (defined($XYcors{$p2Name})){

					my ($xp2, $yp2, $zp2) = split (" ", $XYcors{$p2Name});
			
					# Calculating Minimum Distance WithOut considering Concentrations
					$WithOutConcDist = &CalcDist ($p1Name, $p2Name, $xp1, $xp2, $yp1, $yp2, $zp1, $zp2);

					if ($WithOutConcDist < $WithOutConcMinDist) {
						$WithOutConcMinDist = $WithOutConcDist;
						$MoleculeThatGaveMinDist = $p2Name;
					}
			
					if ($WithOutConcDist == 0 ) { # print "Minimum Distance found... not searching for minimum distance any more\n";
						# print "Exiting at line 152\n"; 
						last EXIT1;
					}
				}

			}
			} # print "\n\n\n";


			# Incorporating Concentrations into the distance measure
			$WithConcMinDist = (($WithOutConcMinDist*10)*(abs($Sample2{$MoleculeThatGaveMinDist} - $cp1)));

			# Storing MinDistances for deriving Hausdroff distance and other metrics
			push (@WithOutConcMinDist, $WithOutConcMinDist);
			push (@WithConcMinDist, $WithConcMinDist);
		}
	}

	# Calculating other basic stats on these minimum distances using R
	open R, "> ./$InputFolderName/CalcStats";
	print R join ("\n", @WithOutConcMinDist);
	close R;
	
	my @BasicStats = `Rscript ./bin/121030_BasicStats.R ./$InputFolderName/CalcStats`;
	@_ = split(":", $BasicStats[0]); chomp($_[1]); $_[1] =~ s/\s+$//; my $WOC_Min = $_[1];
	@_ = split(":", $BasicStats[1]); chomp($_[1]); $_[1] =~ s/\s+$//; my $WOC_Q1 = $_[1];
	@_ = split(":", $BasicStats[2]); chomp($_[1]); $_[1] =~ s/\s+$//; my $WOC_Median = $_[1];
	@_ = split(":", $BasicStats[3]); chomp($_[1]); $_[1] =~ s/\s+$//; my $WOC_Mean = $_[1];
	@_ = split(":", $BasicStats[4]); chomp($_[1]); $_[1] =~ s/\s+$//; my $WOC_Q3 = $_[1];
	@_ = split(":", $BasicStats[5]); chomp($_[1]); $_[1] =~ s/\s+$//; my $WOC_Max = $_[1];
	
	my $WOC_IQR = $WOC_Q3 - $WOC_Q1; 
	
	# Calculating other basic stats on these minimum distances using R
	open R, "> ./$InputFolderName/CalcStats";
	print R join ("\n", @WithConcMinDist);
	close R;
	
	@BasicStats = `Rscript ./bin/121030_BasicStats.R ./$InputFolderName/CalcStats`;
	@_ = split(":", $BasicStats[0]); chomp($_[1]); $_[1] =~ s/\s+$//; my $WC_Min = $_[1];
	@_ = split(":", $BasicStats[1]); chomp($_[1]); $_[1] =~ s/\s+$//; my $WC_Q1 = $_[1];
	@_ = split(":", $BasicStats[2]); chomp($_[1]); $_[1] =~ s/\s+$//; my $WC_Median = $_[1];
	@_ = split(":", $BasicStats[3]); chomp($_[1]); $_[1] =~ s/\s+$//; my $WC_Mean = $_[1];
	@_ = split(":", $BasicStats[4]); chomp($_[1]); $_[1] =~ s/\s+$//; my $WC_Q3 = $_[1];
	@_ = split(":", $BasicStats[5]); chomp($_[1]); $_[1] =~ s/\s+$//; my $WC_Max = $_[1];

	my $WC_IQR = $WC_Q3 - $WC_Q1; 

	return ($WOC_Median, $WOC_Mean, $WOC_Max, $WOC_IQR, $WC_Median, $WC_Mean, $WC_Max, $WC_IQR);

}

#####
sub CalcDist {

	my ($p1Name, $p2Name, $xp1, $xp2, $yp1, $yp2, $zp1, $zp2) = ($_[0], $_[1], $_[2], $_[3], $_[4], $_[5], $_[6], $_[7]);
	my $dist; 

	if ($p1Name =~ m/\b$p2Name\b/) { # print "Entered line 169";<stdin>;
		# means two species are identical by name itself
		$dist=0;
	}

	else { # print "Entered line 174";<stdin>;
			$dist = sqrt((($xp2-$xp1)**2)+(($yp2-$yp1)**2)+(($zp2-$zp1)**2));
	}

	return $dist;
}

#######
sub max {
    my ($max, @vars) = @_;
    for (@vars) {
        $max = $_ if $_ > $max;
    }
    return $max;
}

#######
sub min {
    my ($min, @vars) = @_;
    for (@vars) {
        $min = $_ if $_ < $min;
    }
    return $min; 
}

__END__

	# Taking log10 values of concentrations and scaling them. Scaling is necessary because log10 gives negative values for concentrations. 
	# So, all log values must be moved to positive x-axis.
	# Scaling is done in the following manner : Lowest concentration among Sample1 & Sample2 is selected and log10 value of this Conc is added to 
	# all other ConcLog10 values, effectively moving the scale to positive x-axis.
	
	@TempArray1 = ();
	push (@TempArray1, values(%Sample1), values(%Sample2));
	my $MinConc = min(@TempArray1); # print "MinConc = $MinConc\n";
	my $Log10ofMinConc = ((log($MinConc)/log(20))/100); # print "Log10ofMinConc = $Log10ofMinConc\n";

	foreach (keys (%Sample1)){ # print "var1 in sample1 = $Sample1{$_}\n";
		$Sample1{$_} = ((log($Sample1{$_})/log(20))/100); # print "log10 of var1 in sample1 = $Sample1{$_}\n";
		$Sample1{$_} = $Sample1{$_} + abs($Log10ofMinConc); # print "$Sample1{$_}\n";
	}
	foreach (keys (%Sample2)){
		$Sample2{$_} = ((log($Sample2{$_})/log(20))/100);
		$Sample2{$_} = $Sample2{$_} + abs($Log10ofMinConc); # print "$Sample2{$_}\n";
	}


# After K-means clustering...
	# Input : 1) PCA file after K-means clustering 2)Sample1 3) Sample2

#	$ConcFlag = 0;

#	my $HdistSum = 0;
	# for each cluster ... have to write code for extracting cluster from PCA file ---- depends on output of k-means clustering output.
#		$HdistSum = ($HdistSum + CalcHdist($ConcFlag));
#	my $AvgHdistWoutConc = $HdistSum/$k;

#	$ConcFlag = 1;

#	$HdistSum = 0;

	# for each cluster... have to write code for extracting cluster from PCA file ---- depends on output of k-means clustering output.
#		$HdistSum = ($HdistSum + CalcHdist($ConcFlag));

#	my $AvgHdistWConc = $HdistSum/$k;

