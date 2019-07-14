use warnings;
use strict;
use List::Util qw(min max); # For Finding Max and Min in a list
use Data::Dumper qw(Dumper);

# Syntax for using this program :
# perl 121012_isomers.pl <FolderName>
# Example : perl 121012_isomers.pl 121012_YeastTest

my $MethodChoice = 4;

# By default, Levenshtein Distance is used to determine similarity between molecules. 
# It is possible to use other methods also. To do so, syntax must be changed to..

# perl 121012_isomers.pl <FolderName> <MethodNumber>
# Example : perl 121012_isomers.pl 121012_YeastTest 4

# If you are changing method, un-comment next line
# $MethodChoice = $ARGV[1];

# MethodNumber Choices :
# 1. Fingerprint method
# 2. Word Frequency (LINGO) method
# 3. Bioisosteric method
# 4. Levenshtein distance
# 5. (Smith-Waterman) Local Alignment
# 6. (MUSCLE) Multiple Sequence Alignment

# Saving Input Folder Name
my $InputFolderName = $ARGV[0];

# Along with SMILES containing files, input folder must contain three other files :
# 1. with name "UserSMILES" : A plain text file that must contain SMILES for those lipid species that are not supported by LipidMapsTools. If you are not sure which Lipid species are not supported by LipidMapsTools, run this program with empty "UserSMILES" file and this will program will tell you which lipid species are not supported by LipidMapsTools
# 2. with name "FA" : A plain text file that will tell the program which fatty acids to use.
# 3. with name "Cer" : A plain text file that will tell the program which set of sphingosine(s) to use for drawing Ceramides
# These two files must follow the same formatting as that of example files.

# Changes from previous version of this program : 
# 1. Will take multiple files as input (All files must be placed in a folder and folder name is taken as input from user)
# 2. Will read concentrations from input files
# EDIT 15.04.2015
# Also takes a single tab/comma/space separated file as input

# Caveats :
# 1. This program uses bash commands, so this program works only on linux OS
# 2. The format of molecule abbreviations in input file must be exactly as they are in example files. For instance, PI 18:1-18:1 has to be written as 
#    PI 18:1/18:1. Similarly other lipid species also have certain fixed abbreviation formats. This program will read a molecule only if it is given
#    in correct format !
# 3. This program is meant to generate structure possibilities for lipid species that have defined SN1 and SN2 Chains (Eg. PI 18:1/16:2)
#    For instance, this program cannot generate isomer possibilities for PI 36:2. 
#    If you donot have SN1 and SN2 information, use 120418_isomers.pl instead !
#    /media/second HDD 1/chakravarthy_backup/DesktopDump120717/Dm_lipidome/120418_isomers.pl

# Giving path to LipidMaps Structure Drawing Tools
`export PATH="\$PATH:./bin/*_lipidmapstools/bin"`; # In bash/Linux

# Declaring Variables
my $VarString1;
my $VarString2;
my @VarList1;
my @VarList2;
my ($c1,$c2,$c3,$c4,$d1,$d2,$d3,$d4,$h1,$h2,$h_position,$db_positions1,$db_positions2,$db_positions3,$db_positions4,$i, $c_atoms, $d_bonds);

#c1 -> number of carbon atoms in Fatty Acid 1
#c2 -> number of carbon atoms in Fatty Acid 2
#c3 -> number of carbon atoms in Fatty Acid 3
#c4 -> number of carbon atoms in Fatty Acid 4
#d1 -> number of double bonds in Fatty Acid 1
#d2 -> number of double bonds in Fatty Acid 2
#d3 -> number of double bonds in Fatty Acid 3
#d4 -> number of double bonds in Fatty Acid 4
#h1 -> number of hydroxyl groups in Fatty Acid 1
#h2 -> number of hydroxyl groups in Fatty Acid 2

my @ClassList=(); # Initializing an array to store lipid classes
my @CentroidList=();
my @FAVarList2=(); # Probably not required. Included for Debugging
my $LipidName;
my $LipidNameCopy;
my $centroid;
my $SMILESlistFileName;
my $CentroidListFileName;
my @LipidNamesList=();
my $LipidClass;
my $LipidMapsProgramToUse;
my @Chains; # A Temporary list to store fatty acid chains
my @UnsupportedSpecies; # Array to store unused lipid abbreviations 
my $TempFileName;
my $HighestConc; # Variable to store highest concentration observed in given set of files. Useful for plotting.
my @Lipid; # A list to store all lipid species
my @TempLipidsVarList1;
my @TempLipidsVarList2;
my @TempFileNameList;


my $LysoFlag=0; # If 0, its NOT a lyso species. If 1, its a lyso species.

# Obtaining Names of files present in InputFolder
my @FileNameList = `ls ./$InputFolderName`;

if (scalar(@FileNameList) == 0){ die "\nSomething is wrong, didnt find any files in $InputFolderName folder. Program terminating\n\n";}

# Removing known FileNames ("UserSMILES" and "FA" and "Cer") from @FileNameList
my @TempArray = (); my $file;
foreach $file (@FileNameList) { 
	chomp($file);
	push (@TempArray, $file) if $file !~ m/(UserSMILES|FA|Cer)\.t(xt|ab)/;
}
@FileNameList = ();@FileNameList = @TempArray; 

# Checking for duplicate FileNames in InputFolder
my @UniqFileNameList = uniq(@FileNameList); 
if (scalar(@UniqFileNameList) != scalar(@FileNameList)) { 
	die "\nThere are ",  scalar(@FileNameList) - scalar(@UniqFileNameList), " duplicate file names in your input folder.... program is terminating\n\n";
}

# Saving number of Input files
my $NumOfInputFiles = scalar(@FileNameList);

# Reading UserSMILES file
my %UserSMILES; # Hash table to store UserSMILES
my $count = 0 ;
open UserSMILES_FILE, "./$InputFolderName/UserSMILES.tab";
while (<UserSMILES_FILE>){ $count++;
	$_ =~ s/\R//g; # Removing carraige return and/or new line character
	@_ = split (/\t/,$_); if ( defined($_[0]) or defined($_[1]) ) {} else { die "Unable to parse UserSMILES file at line $count";}
	chomp($_[1]);
	$UserSMILES{$_[1]}="$_[0]";
}
close UserSMILES_FILE; 

# Loading Fatty Acids from File
my @FAlist; # A list containing Fatty Acids
open FA_FILE, "./$InputFolderName/FA.tab" or die "Can\'t open \.\/$InputFolderName\/FA.tab \: $!";
while (<FA_FILE>){
	$_ =~ s/\R//g; # Removing carraige return and/or new line character
	chomp($_);
	@_ = split("\t",$_);
	push (@FAlist, [@_]); 
}
close FA_FILE;

# Loading Sphingosine(s) from File
my @CerList; # A list containing Sphingosoine composition
open CER_FILE, "./$InputFolderName/Cer.txt" or die "Can\'t open \.\/$InputFolderName\/Cer.txt \: $!";
while (<CER_FILE>){
	chomp($_);
	$_ =~ s/\R//g;
	push (@CerList, $_); 
}
close CER_FILE; 

#################################################
#1----------Reading Lipidome Data----------------
#################################################
my @SampleList;
my %DataTable = ();

if ((scalar (@FileNameList)) > 1) {# If multiple input files, read each file in InputFolder and make a combined list of lipid species

	$HighestConc=0; # Variable to store highest concentration observed in given set of files. Useful for plotting.
	foreach $TempFileName (@FileNameList){
		@TempLipidsVarList1=();
		@TempLipidsVarList2=();
		open InFile, "./$InputFolderName/$TempFileName" or die "Can't open file $TempFileName in folder $InputFolderName : $!";
		while (<InFile>){
			$_ =~ s/\R//g; # Removing carraige return and/or new line character
			@_ = split(/\t/,$_); # Separating LipidSpecies Names from their Concentrations
			push @TempLipidsVarList1, $_[0];
			# $_[0] =~ s/(\s|;|:|,)/_/g;
			$DataTable{"$_[0]"}{"$TempFileName"} = $_[1]; # Hash of hash. Autovivification. $DataTable{Molecule Name}{Sample} = Molecule Concentration
			chomp($_[1]);
			if ($HighestConc < $_[1]) { $HighestConc = $_[1];}
		}
		close InFile;

		# Checking if there are any duplicate entries in a given file.
		@TempLipidsVarList2 = uniq(@TempLipidsVarList1);
		if (scalar(@TempLipidsVarList2) != scalar(@TempLipidsVarList1)) { 
			die "\nThere are ",  scalar(@TempLipidsVarList1) - scalar(@TempLipidsVarList2), " duplicate molecule names in file $TempFileName. Program terminating\n\n";
		}
		
		push @Lipid,@TempLipidsVarList1;
	}
	
	
	# Combining lipid species from all files into one single Duplicate-Free list.
	@TempLipidsVarList2 = ();
	@TempLipidsVarList2 = uniq(@Lipid);
	@Lipid = ();
	@Lipid = sort @TempLipidsVarList2; print "No. of unique lipid species = ", scalar(@Lipid), "\n"; 
	my $fh; open $fh, "> ./$InputFolderName/LipidsList.txt"; print $fh join ("\n", @Lipid); close $fh;

	# Writing combined data to a single file
	open $fh, "> ./$InputFolderName/$InputFolderName.tab";
	print $fh join ("\t", @FileNameList);
	my ($MolName, $FileName);
	foreach $MolName (@Lipid){
		print $fh ("\n", $MolName, "\t");
		@TempLipidsVarList2 = ();
		foreach $FileName (@FileNameList) {
			if (defined($DataTable{$MolName}{$FileName})) { push @TempLipidsVarList2, $DataTable{$MolName}{$FileName}; }
			else { push @TempLipidsVarList2, "NA"; }
		}
		print $fh join ("\t", @TempLipidsVarList2);
	}
	close $fh;
	# print "Check files"; <stdin>;
}

elsif ((scalar (@FileNameList)) == 1) { # If only one tab/space/comma delimited input file, read data file and make a list of lipid species, sample names

	$TempFileName = $FileNameList[0];
	$HighestConc=0; # Variable to store highest concentration observed in given set of files. Useful for plotting.

	open my $fh, "./$InputFolderName/$TempFileName" or die "Can't open file $TempFileName in folder $InputFolderName : $!";
	my @lines = map [ split(/[\t,\,]+/) ], <$fh>;
	close $fh;
	
	my $NumSamples = scalar @{$lines[0]}; # print "$NumSamples"; <stdin>;
	my $NumLipids = ((scalar @lines)-1); # print "$NumLipids"; <stdin>;
	for (my $i = 1; $i <= $NumLipids; $i++){
		push @Lipid, @{$lines[$i]}[0];
	} # print $Lipid[0], "\t", $Lipid[-1];<stdin>;
	@{$lines[0]}[-1] =~ s/\R//g;
	@SampleList = @{$lines[0]}; # print $FileNameList[0],"\t", $FileNameList[-1]; <stdin>;

	# Checking if there are any duplicate lipid entries in a given file.
	@TempLipidsVarList1 = uniq(@Lipid);
		if (scalar(@TempLipidsVarList1) != scalar(@Lipid)) { 
			die "\nThere are ", scalar(@TempLipidsVarList1) - scalar(@Lipid), "\nDuplicate molecule names in file $TempFileName. Program terminating\n\n";
		}
	@Lipid = @TempLipidsVarList1;
	
	# Checking if there are any duplicate sample entries in a given file.
	@TempLipidsVarList1 = uniq(@SampleList);
		if (scalar(@TempLipidsVarList1) != scalar(@SampleList)) { 
			die "\nThere are ",  scalar(@TempLipidsVarList1) - scalar(@SampleList), "\nDuplicate sample/experiment names in file $TempFileName. Program terminating\n\n";
		}
	@SampleList = @TempLipidsVarList1;
	
	# Obtaining highest concentration
	my $j;
	for ($i=1; $i<=$NumLipids; $i++) { 
		for ($j=1; $j<=$NumSamples; $j++) {
			if (@{$lines[$i]}[$j] !~ m/NA/){
				if (@{$lines[$i]}[$j] > $HighestConc){$HighestConc = @{$lines[$i]}[$j]; } 
			}
		}
	}
#	print "\n\n\$HighestConc = $HighestConc"; <stdin>;

	# Writing Lipidome data from each sample to a separate file
	$j = 1; my $write2file; my $Filename;
	foreach $Filename (@SampleList){
		open $write2file, "> ./$InputFolderName/$Filename";
		for ($i = 1; $i <= $NumLipids; $i++){
			if (@{$lines[$i]}[$j] !~ /NA/){
				@{$lines[$i]}[$j] =~ s/\R//g;
				print $write2file @{$lines[$i]}[0], "\t", @{$lines[$i]}[$j];
				print $write2file "\n" if ($i < $NumLipids);
			}
		} 
		close $write2file;
		$j++;
	} # print "\nCheck Files.."; <stdin>;
	
}

print "Highest conc = $HighestConc\n";

#################################################
#2---------Generate Lipid Structures-------------
#################################################
my ($PlasmogenFlag, $SN_Flag, $OH_Flag );
my ($dbp1, $dbp2);

open OutFile2, "> ./$InputFolderName/Isomers.txt" or die "Unable to open ./$InputFolderName/Isomers.txt file"; 

foreach (@Lipid){
exit_foreach_lipid_loop: {

	$LipidName=$_;
	chomp($LipidName);
	print "Generating possible structures for $LipidName\n";

	# Opening a file to write isomer possibilities for centroid determination
	open OutFile1, "> ./$InputFolderName/OutFile1.txt";

	$LipidClass = &Get_LipidClass ; 
	# print "Lipid Class = $LipidClass\n";# <stdin>;
	
	$PlasmogenFlag = &CheckPlasmogen;
	$SN_Flag = &Check_SN_composition;
	$OH_Flag = &Get_OH_info_for_Ceramides; # Get information on hydroxylation state for Ceramides. Applicable only if SN1 and SN2 FA composition is not specified.

	my $user_SMILES_found = &Search_UserSMILES;

	if ($user_SMILES_found == 1){
		close OutFile1;
		last exit_foreach_lipid_loop;
	}

	elsif ($LipidClass !~ m/\b(PC|PA|PS|PI|PE|PG|DG|TG|TAG|DAG|CL|IPC|MIPC|MIP2C|Cer|CerP|PICer|HexCer|PECer|CerPE|LCB|LCBP|MIPCer|MIP2Cer|LPC|LPA|LPS|LPI|LPE|LPG|LIPC)\b/){
		push @UnsupportedSpecies, $LipidName;
		last exit_foreach_lipid_loop;
		# print "\nCannot generate structure for this molecule..Press any key to coninue"; <stdin>;
	}

	push (@ClassList, $LipidClass);

	$LysoFlag = &CheckForLyso; # This sub routine also changes the string $LipidClass
	

	# print "\nLipid Class = $LipidClass\n\$PlasmogenFlag = $PlasmogenFlag\n\$SN_Flag = $SN_Flag\n\$LysoFlag = $LysoFlag\n\$OH_Flag = $OH_Flag\n";<stdin>;
	# print "\nPlasmogenFlag = $PlasmogenFlag\n";

	if (($SN_Flag == 1) or ($LysoFlag == 1)) {
		&Parse_and_draw_Isomers;
	}

	else {
		# Predict SN composition and draw PhosphoLipids and DG
		if ($LipidClass =~ m/\b(PC|PA|PS|PI|PE|PG|DG)\b/){ # print "\nPredicting SN composition and drawing PC_DG...\n";
			($c1, $c2, $d1, $d2, $dbp1, $dbp2) = &Predict_SN_and_draw_PL_DG($LipidName, $LipidClass);
		}
		
		# Predict SN composition and draw Ceramides
		elsif ($LipidClass =~ m/\b(Cer|CerP|PICer|HexCer|PECer|MIPCer|MIP2Cer)\b/){ # print "\nPredicting SN composition and drawing Ceramides...\n";
			($c1, $c2, $d1, $d2, $dbp1, $dbp2) = &Predict_SN_and_draw_Cer($LipidName, $LipidClass);
		}

		# Draw TG
		elsif ($LipidClass =~ m/\bTG\b/){ # print "\nPredicting SN composition and drawing TG...\n";
			($c1, $c2, $d1, $d2, $dbp1, $dbp2) = &Predict_SN_and_draw_TG($LipidName, $LipidClass);
		}
		
		# Draw CL
		elsif ($LipidClass =~ m/\bCL\b/){ # print "\nPredicting SN composition and drawing CL...\n";
			($c1, $c2, $d1, $d2, $dbp1, $dbp2) = &Predict_SN_and_draw_CL($LipidName, $LipidClass);
		}
	}

	# Program design generates duplicates for Ceramides, PC/PE TAG etc. (TODO : Modify design so that duplicates are not generated)
	# Deleting duplicates...
	close OutFile1; # print "\nCheck OutFile1 and press any key to continue..."; <stdin>;
	&DelDuplicates;
	&FindCentroidIsomer;
	
}
}
&WriteLists2Files;
close OutFile2;
	
################################################################################
#3------------Calculating pair wise distances between lipid species-------------
################################################################################

my $CSVfilename = "$InputFolderName"."\.csv"; # print "./$InputFolderName/$CSVfilename";<stdin>;
my $LOGfilename = "$InputFolderName"."\.log";
my $ALNfilename = "$InputFolderName"."\.aln";
	
if    ($MethodChoice==1){`perl		./bin/*m1.pl	./$InputFolderName/$SMILESlistFileName				./$InputFolderName/$CSVfilename	./$InputFolderName/$LOGfilename`;}
elsif ($MethodChoice==2){`python	./bin/*m2.py	./$InputFolderName/$SMILESlistFileName				./$InputFolderName/$CSVfilename	./$InputFolderName/$LOGfilename`;}
#elsif ($MethodChoice==3){`perl	./bin/*m3.pl	./$InputFolderName/$SMILESlistFileName				./$InputFolderName/$CSVfilename	./$InputFolderName/$LOGfilename`;}
elsif ($MethodChoice==4){`python	./bin/*m4.py	./$InputFolderName/$SMILESlistFileName				./$InputFolderName/$CSVfilename	./$InputFolderName/$LOGfilename`;}
elsif ($MethodChoice==5){`python	./bin/*m5.py	./$InputFolderName/$SMILESlistFileName	./bin/mat4	./$InputFolderName/$CSVfilename	./$InputFolderName/$LOGfilename`;}
elsif ($MethodChoice==6){`perl		./bin/*m6.pl	./$InputFolderName/$SMILESlistFileName	./$InputFolderName/O2.aln	./$InputFolderName/$CSVfilename	./$InputFolderName/$LOGfilename`;
					`rm ./$InputFolderName/O2.aln ./$InputFolderName/O2.aln_joined`;`mv ./$InputFolderName/O2.aln_joined_rearranged ./$InputFolderName/$ALNfilename`;
}

################################################################################
#4--------------------Performing Principal Component Analysis ------------------
################################################################################

my $PCAfilename = "$InputFolderName"."\.pca";
my $VariancePlotFileName = "$InputFolderName"."_VariancePlot"."\.pdf";
my @PCAcoordinates = `Rscript ./bin/121012_PCAfromCSV.R ./$InputFolderName/$CSVfilename ./$InputFolderName/$VariancePlotFileName ./$InputFolderName/$PCAfilename`;

chomp($PCAcoordinates[0]); #print "$PCAcoordinates[0] \n";
chomp($PCAcoordinates[1]); #print "$PCAcoordinates[1] \n";
chomp($PCAcoordinates[2]); #print "$PCAcoordinates[2] \n";
chomp($PCAcoordinates[3]); #print "$PCAcoordinates[3] \n $HighestConc \n"; <stdin>;

# print "\n\nCompleted PCA.. press any key to continue"; <stdin>;
################################################
#5---------Plot 2D and 3D Maps------------------
################################################

my @MultiPlotOutput=();

if (scalar (@FileNameList) == 1) { @FileNameList = (); @FileNameList = @SampleList; }

foreach $TempFileName (@FileNameList){
	@MultiPlotOutput = `Rscript ./bin/130115_PCA_plot_150507Edit2.R "./$InputFolderName/$TempFileName" "./$InputFolderName/$InputFolderName\.pca" "./$InputFolderName/ClassList.txt" "./$InputFolderName/CentroidList.txt" "$PCAcoordinates[0]" "$PCAcoordinates[1]" "$PCAcoordinates[2]" "$PCAcoordinates[3]" "$HighestConc" "./$InputFolderName/$TempFileName"`;
	# print @MultiPlotOutput;
}

# Combined Lipidome Plot
my $SVGfilename_2D = "$InputFolderName"."_2D\.svg";
my $SVGfilename_3D = "$InputFolderName"."_3D\.svg";
`Rscript ./bin/130124_CombinedLipidomePlot_150507Edit2.R "./$InputFolderName/$SVGfilename_2D" "./$InputFolderName/$PCAfilename" "./$InputFolderName/ClassList.txt" "./$InputFolderName/$SVGfilename_3D"`;

# Writing un-analyzed lipid species to <stdout>
print "\nUnused lipid abbreviations\n\n";
print join ("\n", @UnsupportedSpecies);
print "\n\n";

open my $fh, "> ./$InputFolderName/Unused_LipidAbbreviations.txt";
print $fh join ("\n", @UnsupportedSpecies);
close $fh;

################################################################################################
#6----------------Calculating Hausdroff distance and ploting Historgrams------------------------
################################################################################################

my $FileNameListRef = \@FileNameList;
my $FileNameList = join("XXXX", @FileNameList); # print "$FileNameList\n";

# Calculating Hausdroff distance and Average distance from PC1 PC2 Coordinates
 `perl ./bin/130117_3_MaxDist_and_AvgDist_from_PC1_PC2.pl $InputFolderName $PCAfilename $FileNameList`; 

# Calculating Hausdroff distance and Average distance from PC1 PC2 PC3 Coordinates
# $FileNameList = join("XXXX", @FileNameList);
 `perl ./bin/130123_2_MaxDist_and_AvgDist_from_PC1_PC2_PC3.pl $PCAfilename $FileNameList $InputFolderName`;

# Calculating Hausdroff distance and Average distance from Similarity Scores
# $FileNameList = join("XXXX", @FileNameList);
 `perl ./bin/140805_2_MaxDist_and_AvgDist_from_SimilarityScore.pl $CSVfilename $FileNameList $InputFolderName`;

# Creating Dendrograms as pdf - one per page
 `Rscript ./bin/130117_Hdist2Dendro_150629_Edit.R $InputFolderName "./$InputFolderName/LUX_PC1_PC2.csv" "./$InputFolderName/LUX_PC1_PC2_PC3.csv" "./$InputFolderName/LUX_SimilarityScore.csv"`;

`rm ./$InputFolderName/CalcStats`;
###########################################
#----------- Sub-routine's ----------------
###########################################

sub uniq {

	my @InArray = @_;
	my @unique = ();
	my %seen   = ();
	foreach my $elem ( @InArray ) {
		next if $seen{ "$elem" }++;
		push @unique, $elem;
	}

	return (@unique);
	# return keys %{{ map { $_ => 1 } @_ }};  
}   

sub ModifySMILES {

	$_ = $_[0];

	$_ =~ s/\[C@\@H\]/C/g;			# Removing steriochemistry for chiral carbon
	$_ =~ s/\[C\@H\]/C/g;
	$_ =~ s/\[C\@\@]/C/g;
	$_ =~ s/\[C\@]/C/g;

	$_ =~ s/\\//g;					# Forward and backward slashes (meaning - cis - trans) are removed
	$_ =~ s/\///g;

	$_ =~ s/\[O\-\]/O/g;			# Removing charges for Oxygen and Nitrogen atoms
	$_ =~ s/\[N\+\]/N/g;	

	# Carbon, Nitrogen, Phosphorus, Sulphur, Fluorine and Iodine retain their characters
#	$_ =~ s/O|\[O\-\]/A/g;			# Since character 'O' is not a defined animo acid, Oxyzen is given 'E' character
#	$_ =~ s/Na|\[Na\]|\[Na\+\]/D/g;	# Two character atoms (Na, Br, Cl) are given single character
#	$_ =~ s/Br|\[Br\]/E/g;			
#	$_ =~ s/Cl|\[Cl\]|\[Cl\-\]/G/g;
#	$_ =~ tr/=#/HK/;				# Double and triple bonds are given G and H characters respectively
#	$_ =~ tr/()/LM/;				# Round brackets (meaning - branch opening and closing) are given K and L characters
#	$_ =~ tr/cn/CN/;				# Carbon, Nitrogen in ring structures (represented in small letters) and converted to capital letters
#	$_ =~ tr/H//d;					# Hydrogen atoms are removed
#	$_ =~ tr/[+,\-,*,.]//d;			# Charges (+,-), wild character (*) and join character (.) are removed
#	$_ =~ tr/1-9//d;				# Numbers 1-9 are removed
#	$_ =~ tr/[]//d;				# Square brackets are removed; These are used for differentiating special atoms
#	$_ =~ tr/@//d;					# @ (specifies position in with respect to chiral centre) is removed

	return ($_);
}

sub Search_UserSMILES {
	my ($key,$LipidNameCopy);
	my $count = 0;
	$LipidNameCopy = $LipidName;
	foreach $key (keys (%UserSMILES)) { 
		if ($key =~ m/\b\Q$LipidName\E\b/){ $count++;
			# Add2Lists
			$LipidNameCopy =~ s/(\s|;|:|,)/_/g;
			push (@CentroidList, "$LipidNameCopy\tUserDefinedSMILES\t0"); 
			push (@LipidNamesList, "$UserSMILES{$LipidName} $LipidNameCopy");
			print "User defined SMILES taken for this molecule\n";
		}
	}
	
	if ($count > 1){
		print "\nFound duplicate molecule name in user defined SMILES list\n" and die "Program Terminating";
	}
	
	return ($count);
}

sub Get_LipidClass {
	if ($LipidName =~ m/\s/){
		@_ = split /(\s|_)/ , "$LipidName";
		$LipidClass = $_[0]; #print "\$LipidClass at line 438  = $LipidClass";<stdin>;
	}
}


sub CheckPlasmogen {
	$PlasmogenFlag = 0;
	if ($LipidName =~ m/_O/){ # print "Encountered O-";<stdin>;
		$PlasmogenFlag = 1;
		# $LipidName =~ s/\sO-/\_O/; # print $LipidClassDup;<stdin>;
	}
	if ($LipidName =~ m/_P/){ 
		$PlasmogenFlag=2;
		# $LipidName =~ s/\sP-/\_P/;
	}
	return ($PlasmogenFlag);
}

sub Check_SN_composition {

	@_ = split /\s/ , "$LipidName";
	my @SN_Chains = split /\// , $_[1] if (defined($_[1]));

	if (defined($SN_Chains[1])){
		$SN_Flag = 1;
	}
	else {
		$SN_Flag = 0;
	}
	
	return ($SN_Flag);
}

sub CheckForLyso { 	# Finding out of a lipid species is a lyso molecule or not.

	$LysoFlag = 0;
	if ($LipidClass =~ m/^TAG/) {$LipidClass =~ s/TAG/TG/;}
	elsif ($LipidClass =~ m/^DAG/) {$LipidClass =~ s/DAG/DG/;}
	elsif ($LipidClass =~ m/^IPC/) {$LipidClass =~ s/IPC/PICer/;}
	elsif ($LipidClass =~ m/^CerPE/) {$LipidClass =~ s/CerPE/PECer/;}
	elsif ($LipidClass =~ m/^MIPC/) {$LipidClass =~ s/^MIPC/MIPCer/;}
	elsif ($LipidClass =~ m/^M\(IP\)2C/) {$LipidClass =~ s/^M\(IP\)2C/MIP2Cer/;}
	elsif ($LipidClass =~ m/^MIP2C/) {$LipidClass =~ s/^MIP2C/MIP2Cer/;}
	elsif ($LipidClass =~ m/^LCB$/) {$LipidClass =~ s/LCB/Cer/; $LysoFlag=1;} # Its a lyso molecule
	elsif ($LipidClass =~ m/^LCBP$/) {$LipidClass =~ s/LCBP/CerP/; $LysoFlag=1;} # Its a lyso molecule
	elsif ($LipidClass =~ m/^LPC/) {$LipidClass =~ s/LPC/PC/; $LysoFlag=1;} # Its a lyso molecule
	elsif ($LipidClass =~ m/^LPA/) {$LipidClass =~ s/LPA/PA/; $LysoFlag=1;} # Its a lyso molecule
	elsif ($LipidClass =~ m/^LPS/) {$LipidClass =~ s/LPS/PS/; $LysoFlag=1;} # Its a lyso molecule
	elsif ($LipidClass =~ m/^LPI/) {$LipidClass =~ s/LPI/PI/; $LysoFlag=1;} # Its a lyso molecule
	elsif ($LipidClass =~ m/^LPE/) {$LipidClass =~ s/LPE/PE/; $LysoFlag=1;} # Its a lyso molecule
	elsif ($LipidClass =~ m/^LIPC/) {$LipidClass =~ s/LIPC/PICer/; $LysoFlag=1;} # Its a lyso molecule
	
	return ($LysoFlag);
}

sub Get_OH_info_for_Ceramides{
	$OH_Flag = 0;
	if ($LipidName =~ m/C_dt/) { 
		$OH_Flag = 3;
		# $LipidClass =~ s/C_dt/C/;print "\n\$LipidClass = $LipidClass\n";
	}
	elsif ($LipidName =~ m/C_d/) {
		$OH_Flag = 1;
		# $LipidClass =~ s/C_d/C/;print "\n\$LipidClass = $LipidClass\n";
	}
	elsif ($LipidName =~ m/C_t/) {
		$OH_Flag = 2;
		# $LipidClass =~ s/C_t/C/;print "\n\$LipidClass = $LipidClass\n";
	}
	
	return ($OH_Flag);
}

sub Parse_and_draw_Isomers {
	
	chomp($LipidName);
	my $LipidNameCopy = $LipidName; # making a copy of lipid name to use it later (while generating CentroidList and SMILESlist files )

	# Breaking lipid name to extract fatty acid chain information
	my @LipidNameSplit;
	if ($LipidName =~ m/\s/){ # Any lipid except sterols
		@LipidNameSplit = split /\s/ , "$LipidName";
		@Chains = split /\// , $LipidNameSplit[1];
	}
	if (($LipidClass =~ m/Cer/) and defined($LipidNameSplit[2])){ # print "\n \$LipidNameSplit\[2\] = $LipidNameSplit[2]\n"; 
		if ($LipidNameSplit[2] == 1){$h2 = 1;}
		elsif ($LipidNameSplit[2] == 2){$h2 = 2;}
		else {print "\nUnknown text found in Lipid name\n" and die "Program Terminating !!";}
	}

	# Identifying number of c'atoms and d'bonds for (L)P(C|A|S|I|E|G) and DaG and Cer(P) CerPE HexCer
	if ($LipidClass =~ m/\b(PC|PA|PS|PI|PE|PG|DG|Cer|CerP|PICer|HexCer|PECer|MIPCer|MIP2Cer)\b/){

		# Splitting string "chain_1" to obtain 1) No. of carbon atoms 2) No. of double bonds
		@_ = split("\:|\;", $Chains[0]);
		($c1,$d1) = ($_[0], $_[1]); # print "\nc1 = $c1 and d1 = $d1\n";<stdin>;
		if (defined($_[2]) and $LipidClass =~ m/Cer/) { 
			if ($_[2] == 1){$h1 = 0;}
			elsif ($_[2] == 2){$h1 = 1;}
			elsif ($_[2] == 3){$h1 = 2;}
			else {print "\nOH position in LongChainBase is not parsed for molecule $LipidName\n" and die "Program Terminating\n";}
			# print "\$_\[2\] = $_[2]   \$h1 = $h1\n";
		} # Applicable to Cermaides
		
		if ($LysoFlag == 1) {
			$c2 = 0;
			$d2 = 0;
			if ($LipidClass =~ m/Cer/){ $h2 = 0;}
		}
		
		# Splitting string "chain_2" to obtain 1) No. of carbon atoms 2) No. of double bonds
		elsif ($LysoFlag == 0) {
			@_ = split("\:|\;", $Chains[1]);
			($c2,$d2) = ($_[0], $_[1]);
			if (($LipidClass =~ m/Cer/) and defined($_[2])){ $h2 = $_[2];}
		}
		
		foreach (@FAlist){
			if (($_->[0] == $c1) and ($_->[1] == $d1)) {
				if($d1>0){ 
					if (defined($_->[2])){ 
						$db_positions1= $_->[2];
					}
					else {
						$db_positions1="";
					}
				}else{
					$db_positions1="";
				}
				foreach (@FAlist){
					if (($_->[0] == $c2) and ($_->[1] == $d2)) {
						if($d2>0){
							if (defined($_->[2])){
								$db_positions2= $_->[2];
							}
							else {
								$db_positions2="";
							}
						}else{
							$db_positions2="";
						}
						if ($LipidClass =~ m/DG/){
							print OutFile1 "$LipidClass\($c1\:$d1$db_positions1\/$c2\:$d2$db_positions2\/0\:0\)\n";
							print OutFile1 "$LipidClass\($c2\:$d2$db_positions2\/$c1\:$d1$db_positions1\/0\:0\)\n"; # Xchanging SN1 and SN2 positions
						}elsif ($LipidClass =~ m/Cer/){
							# print "c1 = $c1, d1 = $d1,db_positions1 = $db_positions1, c2 = $c2,d2 = $d2, db_positions2 = $db_positions2";
							&GoToCer;
						}else{
							print OutFile1 "$LipidClass\($c1\:$d1$db_positions1\/$c2\:$d2$db_positions2\)\n";
							# If its not a Lyso molecule, then exchange SN1 and SN2 positions
							if ($LysoFlag != 1){ 
								# Xchanging SN1 and SN2 positions
								print OutFile1 "$LipidClass\($c2\:$d2$db_positions2\/$c1\:$d1$db_positions1\)\n"; 
							}
						}
					}
				}
			}
		}
	}
	
	# For TaG's only
	elsif ($LipidClass =~ m/TG/){
		
		# Splitting string "chain_1" to obtain 1) No. of carbon atoms 2) No. of double bonds
		@_ = split("\:|\;", $Chains[0]);	
		($c1,$d1) = ($_[0], $_[1]);
		
		# Splitting string "chain_2" to obtain 1) No. of carbon atoms 2) No. of double bonds
		@_ = split("\:|\;", $Chains[1]);
		($c2,$d2) = ($_[0], $_[1]);
		
		# Splitting string "chain_2" to obtain 1) No. of carbon atoms 2) No. of double bonds
		@_ = split("\:|\;", $Chains[2]);
		($c3,$d3) = ($_[0], $_[1]);

		foreach (@FAlist){
			if (($_->[0] == $c1) and ($_->[1] == $d1)) {
				if($d1>0){$db_positions1= $_->[2];}else{$db_positions1="";}
				foreach (@FAlist){
					if (($_->[0] == $c2) and ($_->[1] == $d2)) {
						if($d2>0){$db_positions2= $_->[2];}else{$db_positions2="";}
						foreach (@FAlist){
							if (($_->[0] == $c3) and ($_->[1] == $d3)) {
							if($d3>0){$db_positions3= $_->[2];}else{$db_positions3="";}
								print OutFile1 "$LipidClass\($c1\:$d1$db_positions1\/$c2\:$d2$db_positions2\/$c3\:$d3$db_positions3\)\n";
								print OutFile1 "$LipidClass\($c1\:$d1$db_positions1\/$c3\:$d3$db_positions3\/$c2\:$d2$db_positions2\)\n";
								print OutFile1 "$LipidClass\($c2\:$d2$db_positions2\/$c1\:$d1$db_positions1\/$c3\:$d3$db_positions3\)\n";
								print OutFile1 "$LipidClass\($c2\:$d2$db_positions2\/$c3\:$d3$db_positions3\/$c1\:$d1$db_positions1\)\n";
								print OutFile1 "$LipidClass\($c3\:$d3$db_positions3\/$c1\:$d1$db_positions1\/$c2\:$d2$db_positions2\)\n";
								print OutFile1 "$LipidClass\($c3\:$d3$db_positions3\/$c2\:$d2$db_positions2\/$c1\:$d1$db_positions1\)\n";
							}
						}
					}
				}
			}
		}
	}									

	# For CL's only
	elsif ($LipidClass =~ m/CL/){

		# Splitting string "chain_1" to obtain 1) No. of carbon atoms 2) No. of double bonds
		@_ = split("\:|\;", $Chains[0]);
		($c1,$d1) = ($_[0], $_[1]);
		
		# Splitting string "chain_2" to obtain 1) No. of carbon atoms 2) No. of double bonds
		@_ = split("\:|\;", $Chains[1]);
		($c2,$d2) = ($_[0], $_[1]);
		
		# Splitting string "chain_2" to obtain 1) No. of carbon atoms 2) No. of double bonds
		@_ = split("\:|\;", $Chains[2]);
		($c3,$d3) = ($_[0], $_[1]);

		# Splitting string "chain_2" to obtain 1) No. of carbon atoms 2) No. of double bonds
		@_ = split("\:|\;", $Chains[3]);
		($c4,$d4) = ($_[0], $_[1]);

		foreach (@FAlist){
			if (($_->[0] == $c1) and ($_->[1] == $d1)) {
				if($d1>0){$db_positions1= $_->[2];}else{$db_positions1="";}
				foreach (@FAlist){
					if (($_->[0] == $c2) and ($_->[1] == $d2)) {
					if($d2>0){$db_positions2= $_->[2];}else{$db_positions2="";}
					foreach (@FAlist){
							if (($_->[0] == $c3) and ($_->[1] == $d3)) {
							if($d3>0){$db_positions3= $_->[2];}else{$db_positions3="";}
								foreach (@FAlist){
									if (($_->[0] == $c4) and ($_->[1] == $d4)) {
									if($d4>0){$db_positions4= $_->[2];}else{$db_positions4="";}

										print OutFile1 "$LipidClass\(1\'\-\[$c1\:$d1$db_positions1\/$c2\:$d2$db_positions2\]\, 3\'\-\[$c3\:$d3$db_positions3\/$c4\:$d4$db_positions4\]\)\n";
										print OutFile1 "$LipidClass\(1\'\-\[$c1\:$d1$db_positions1\/$c2\:$d2$db_positions2\]\, 3\'\-\[$c4\:$d4$db_positions4\/$c3\:$d3$db_positions3\]\)\n";
										print OutFile1 "$LipidClass\(1\'\-\[$c1\:$d1$db_positions1\/$c3\:$d3$db_positions3\]\, 3\'\-\[$c2\:$d2$db_positions2\/$c4\:$d4$db_positions4\]\)\n";
										print OutFile1 "$LipidClass\(1\'\-\[$c1\:$d1$db_positions1\/$c3\:$d3$db_positions3\]\, 3\'\-\[$c4\:$d4$db_positions4\/$c2\:$d2$db_positions2\]\)\n";
										print OutFile1 "$LipidClass\(1\'\-\[$c1\:$d1$db_positions1\/$c4\:$d4$db_positions4\]\, 3\'\-\[$c2\:$d2$db_positions2\/$c3\:$d3$db_positions3\]\)\n";
										print OutFile1 "$LipidClass\(1\'\-\[$c1\:$d1$db_positions1\/$c4\:$d4$db_positions4\]\, 3\'\-\[$c3\:$d3$db_positions3\/$c2\:$d2$db_positions2\]\)\n";

										print OutFile1 "$LipidClass\(1\'\-\[$c2\:$d2$db_positions2\/$c1\:$d1$db_positions1\]\, 3\'\-\[$c3\:$d3$db_positions3\/$c4\:$d4$db_positions4\]\)\n";
										print OutFile1 "$LipidClass\(1\'\-\[$c2\:$d2$db_positions2\/$c1\:$d1$db_positions1\]\, 3\'\-\[$c4\:$d4$db_positions4\/$c3\:$d3$db_positions3\]\)\n";
										print OutFile1 "$LipidClass\(1\'\-\[$c2\:$d2$db_positions2\/$c3\:$d3$db_positions3\]\, 3\'\-\[$c1\:$d1$db_positions1\/$c4\:$d4$db_positions4\]\)\n";
										print OutFile1 "$LipidClass\(1\'\-\[$c2\:$d2$db_positions2\/$c3\:$d3$db_positions3\]\, 3\'\-\[$c4\:$d4$db_positions4\/$c1\:$d1$db_positions1\]\)\n";
										print OutFile1 "$LipidClass\(1\'\-\[$c2\:$d2$db_positions2\/$c4\:$d4$db_positions4\]\, 3\'\-\[$c1\:$d1$db_positions1\/$c3\:$d3$db_positions3\]\)\n";
										print OutFile1 "$LipidClass\(1\'\-\[$c2\:$d2$db_positions2\/$c4\:$d4$db_positions4\]\, 3\'\-\[$c3\:$d3$db_positions3\/$c1\:$d1$db_positions1\]\)\n";

										print OutFile1 "$LipidClass\(1\'\-\[$c3\:$d3$db_positions3\/$c1\:$d1$db_positions1\]\, 3\'\-\[$c2\:$d2$db_positions2\/$c4\:$d4$db_positions4\]\)\n";
										print OutFile1 "$LipidClass\(1\'\-\[$c3\:$d3$db_positions3\/$c1\:$d1$db_positions1\]\, 3\'\-\[$c4\:$d4$db_positions4\/$c2\:$d2$db_positions2\]\)\n";
										print OutFile1 "$LipidClass\(1\'\-\[$c3\:$d3$db_positions3\/$c2\:$d2$db_positions2\]\, 3\'\-\[$c1\:$d1$db_positions1\/$c4\:$d4$db_positions4\]\)\n";
										print OutFile1 "$LipidClass\(1\'\-\[$c3\:$d3$db_positions3\/$c2\:$d2$db_positions2\]\, 3\'\-\[$c4\:$d4$db_positions4\/$c1\:$d1$db_positions1\]\)\n";
										print OutFile1 "$LipidClass\(1\'\-\[$c3\:$d3$db_positions3\/$c4\:$d4$db_positions4\]\, 3\'\-\[$c1\:$d1$db_positions1\/$c2\:$d2$db_positions2\]\)\n";
										print OutFile1 "$LipidClass\(1\'\-\[$c3\:$d3$db_positions3\/$c4\:$d4$db_positions4\]\, 3\'\-\[$c2\:$d2$db_positions2\/$c1\:$d1$db_positions1\]\)\n";

										print OutFile1 "$LipidClass\(1\'\-\[$c4\:$d4$db_positions4\/$c1\:$d1$db_positions1\]\, 3\'\-\[$c2\:$d2$db_positions2\/$c3\:$d3$db_positions3\]\)\n";
										print OutFile1 "$LipidClass\(1\'\-\[$c4\:$d4$db_positions4\/$c1\:$d1$db_positions1\]\, 3\'\-\[$c3\:$d3$db_positions3\/$c2\:$d2$db_positions2\]\)\n";
										print OutFile1 "$LipidClass\(1\'\-\[$c4\:$d4$db_positions4\/$c2\:$d2$db_positions2\]\, 3\'\-\[$c1\:$d1$db_positions1\/$c3\:$d3$db_positions3\]\)\n";
										print OutFile1 "$LipidClass\(1\'\-\[$c4\:$d4$db_positions4\/$c2\:$d2$db_positions2\]\, 3\'\-\[$c3\:$d3$db_positions3\/$c1\:$d1$db_positions1\]\)\n";
										print OutFile1 "$LipidClass\(1\'\-\[$c4\:$d4$db_positions4\/$c3\:$d3$db_positions3\]\, 3\'\-\[$c1\:$d1$db_positions1\/$c2\:$d2$db_positions2\]\)\n";
										print OutFile1 "$LipidClass\(1\'\-\[$c4\:$d4$db_positions4\/$c3\:$d3$db_positions3\]\, 3\'\-\[$c2\:$d2$db_positions2\/$c1\:$d1$db_positions1\]\)\n";
									}
								}
							}
						}
					}
				}
			}
		}
	}
}

sub GoToCer { 
	if (defined($h2)){
		if ($h2==0){ $h_position = "";}
		elsif ($h2==1){ $h_position = "(2OH)";}
		elsif ($h2==2){ $h_position = "(2OH,4OH)";}
	}	
	else {$h_position = "";} # print "\$LipidClass = $LipidClass";<stdin>;

	if ($h1==0){ print OutFile1 "$LipidClass\(m$c1\:$d1$db_positions1\/$c2\:$d2$db_positions2$h_position\)\n";}
	elsif ($h1==1){ print OutFile1 "$LipidClass\(d$c1\:$d1$db_positions1\/$c2\:$d2$db_positions2$h_position\)\n";}
	elsif ($h1==2){ print OutFile1 "$LipidClass\(t$c1\:$d1$db_positions1\/$c2\:$d2$db_positions2$h_position\)\n";}
}

sub Predict_SN_and_draw_PL_DG{

	my $PlasmogenText;
	if ($PlasmogenFlag == 0){ $PlasmogenText = "";}
	elsif ($PlasmogenFlag == 1){ $PlasmogenText = "O-";}
	elsif ($PlasmogenFlag == 2){ $PlasmogenText = "P-";}

	# Breaking lipid name to extract fatty acid chain information
	if ($LipidName =~ m/\s/){
		@_ = split /\s/ , "$LipidName";
		@Chains = split /\// , $_[1];
	}

	# Splitting string "chain_1" to obtain 1) No. of carbon atoms 2) No. of double bonds
	@_ = split("\:|\;", $Chains[0]);
	($c_atoms,$d_bonds) = ($_[0], $_[1]);

	foreach (@FAlist){
		($c1, $d1, $db_positions1) = ($_->[0], $_->[1], $_->[2]); # NoC1 = No. of Carbon atoms in FA 1
		if (defined($db_positions1)){}else{$db_positions1="";}
		foreach (@FAlist){
			($c2, $d2, $db_positions2) = ($_->[0], $_->[1], $_->[2]);
			if (defined($db_positions2)){}else{$db_positions2="";}
			if ( (($c1+$c2)== $c_atoms) and (($d1+$d2) == $d_bonds) and ($c1!=0) and ($c1!=0) ) {
				if ($LipidClass =~ m/DG/){
					print OutFile1 "$LipidClass\($c1\:$d1$db_positions1\/$c2\:$d2$db_positions2\/0\:0\)\n";
					print OutFile1 "$LipidClass\($c2\:$d2$db_positions2\/$c1\:$d1$db_positions1\/0\:0\)\n"; # Xchanging SN1 and SN2 positions
				}else{
					print OutFile1 "$LipidClass\($PlasmogenText$c1\:$d1$db_positions1\/$c2\:$d2$db_positions2\)\n";
					# If its not a Lyso molecule, then exchange SN1 and SN2 positions
					if ($LysoFlag != 1){ 
						# Xchanging SN1 and SN2 positions
						print OutFile1 "$LipidClass\($PlasmogenText$c2\:$d2$db_positions2\/$c1\:$d1$db_positions1\)\n"; 
					}
				}
			}
		}
	}
}

sub Predict_SN_and_draw_TG{
	# Breaking lipid name to extract fatty acid chain information
	if ($LipidName =~ m/\s/){
		@_ = split /\s/ , "$LipidName";
		@Chains = split /\// , $_[1];
	}

	# Splitting string "chain_1" to obtain 1) No. of carbon atoms 2) No. of double bonds
	@_ = split("\:|\;", $Chains[0]);
	($c_atoms,$d_bonds) = ($_[0], $_[1]);

	foreach (@FAlist){
		($c1, $d1, $db_positions1) = ($_->[0], $_->[1], $_->[2]); # NoC1 = No. of Carbon atoms in FA 1
		if (defined($db_positions1)){}else{$db_positions1="";}
		foreach (@FAlist){
			($c2, $d2, $db_positions2) = ($_->[0], $_->[1], $_->[2]); 
			if (defined($db_positions2)){}else{$db_positions2="";}
			foreach (@FAlist){
				($c3, $d3, $db_positions3) = ($_->[0], $_->[1], $_->[2]); 
				if (defined($db_positions3)){}else{$db_positions3="";}
				if ( (($c1+$c2+$c3)== $c_atoms) and (($d1+$d2+$d3) == $d_bonds) and ($c1!=0) and ($c2!=0) and ($c3!=0) ) {
					print OutFile1 "$LipidClass\($c1\:$d1$db_positions1\/$c2\:$d2$db_positions2\/$c3\:$d3$db_positions3\)\n";
					print OutFile1 "$LipidClass\($c1\:$d1$db_positions1\/$c3\:$d3$db_positions3\/$c2\:$d2$db_positions2\)\n";
					print OutFile1 "$LipidClass\($c2\:$d2$db_positions2\/$c1\:$d1$db_positions1\/$c3\:$d3$db_positions3\)\n";
					print OutFile1 "$LipidClass\($c2\:$d2$db_positions2\/$c3\:$d3$db_positions3\/$c1\:$d1$db_positions1\)\n";
					print OutFile1 "$LipidClass\($c3\:$d3$db_positions3\/$c1\:$d1$db_positions1\/$c2\:$d2$db_positions2\)\n";
					print OutFile1 "$LipidClass\($c3\:$d3$db_positions3\/$c2\:$d2$db_positions2\/$c1\:$d1$db_positions1\)\n";
				}
			}
		}
	}
}	

sub Predict_SN_and_draw_CL{
	# Breaking lipid name to extract fatty acid chain information
	if ($LipidName =~ m/\s/){
		@_ = split /\s/ , "$LipidName";
		@Chains = split /\// , $_[1];
	}

	# Splitting string "chain_1" to obtain 1) No. of carbon atoms 2) No. of double bonds
	@_ = split("\:|\;", $Chains[0]);
	($c_atoms,$d_bonds) = ($_[0], $_[1]);

	foreach (@FAlist){
		($c1, $d1, $db_positions1) = ($_->[0], $_->[1], $_->[2]); # NoC1 = No. of Carbon atoms in FA 1
		if (defined($db_positions1)){}else{$db_positions1="";}
		foreach (@FAlist){
			($c2, $d2, $db_positions2) = ($_->[0], $_->[1], $_->[2]); 
			if (defined($db_positions2)){}else{$db_positions2="";}
			foreach (@FAlist){
				($c3, $d3, $db_positions3) = ($_->[0], $_->[1], $_->[2]); 
				if (defined($db_positions3)){}else{$db_positions3="";}
				foreach (@FAlist){
					($c4, $d4, $db_positions4) = ($_->[0], $_->[1], $_->[2]); 
					if (defined($db_positions4)){}else{$db_positions4="";}
					if ( (($c1+$c2+$c3+$c4)== $c_atoms) and (($d1+$d2+$d3+$d4) == $d_bonds) and ($c1!=0) and ($c2!=0) and ($c3!=0) and ($c4!=0)) {
						print OutFile1 "$LipidClass\(1\'\-\[$c1\:$d1$db_positions1\/$c2\:$d2$db_positions2\]\, 3\'\-\[$c3\:$d3$db_positions3\/$c4\:$d4$db_positions4\]\)\n";
						print OutFile1 "$LipidClass\(1\'\-\[$c1\:$d1$db_positions1\/$c2\:$d2$db_positions2\]\, 3\'\-\[$c4\:$d4$db_positions4\/$c3\:$d3$db_positions3\]\)\n";
						print OutFile1 "$LipidClass\(1\'\-\[$c1\:$d1$db_positions1\/$c3\:$d3$db_positions3\]\, 3\'\-\[$c2\:$d2$db_positions2\/$c4\:$d4$db_positions4\]\)\n";
						print OutFile1 "$LipidClass\(1\'\-\[$c1\:$d1$db_positions1\/$c3\:$d3$db_positions3\]\, 3\'\-\[$c4\:$d4$db_positions4\/$c2\:$d2$db_positions2\]\)\n";
						print OutFile1 "$LipidClass\(1\'\-\[$c1\:$d1$db_positions1\/$c4\:$d4$db_positions4\]\, 3\'\-\[$c2\:$d2$db_positions2\/$c3\:$d3$db_positions3\]\)\n";
						print OutFile1 "$LipidClass\(1\'\-\[$c1\:$d1$db_positions1\/$c4\:$d4$db_positions4\]\, 3\'\-\[$c3\:$d3$db_positions3\/$c2\:$d2$db_positions2\]\)\n";

						print OutFile1 "$LipidClass\(1\'\-\[$c2\:$d2$db_positions2\/$c1\:$d1$db_positions1\]\, 3\'\-\[$c3\:$d3$db_positions3\/$c4\:$d4$db_positions4\]\)\n";
						print OutFile1 "$LipidClass\(1\'\-\[$c2\:$d2$db_positions2\/$c1\:$d1$db_positions1\]\, 3\'\-\[$c4\:$d4$db_positions4\/$c3\:$d3$db_positions3\]\)\n";
						print OutFile1 "$LipidClass\(1\'\-\[$c2\:$d2$db_positions2\/$c3\:$d3$db_positions3\]\, 3\'\-\[$c1\:$d1$db_positions1\/$c4\:$d4$db_positions4\]\)\n";
						print OutFile1 "$LipidClass\(1\'\-\[$c2\:$d2$db_positions2\/$c3\:$d3$db_positions3\]\, 3\'\-\[$c4\:$d4$db_positions4\/$c1\:$d1$db_positions1\]\)\n";
						print OutFile1 "$LipidClass\(1\'\-\[$c2\:$d2$db_positions2\/$c4\:$d4$db_positions4\]\, 3\'\-\[$c1\:$d1$db_positions1\/$c3\:$d3$db_positions3\]\)\n";
						print OutFile1 "$LipidClass\(1\'\-\[$c2\:$d2$db_positions2\/$c4\:$d4$db_positions4\]\, 3\'\-\[$c3\:$d3$db_positions3\/$c1\:$d1$db_positions1\]\)\n";

						print OutFile1 "$LipidClass\(1\'\-\[$c3\:$d3$db_positions3\/$c1\:$d1$db_positions1\]\, 3\'\-\[$c2\:$d2$db_positions2\/$c4\:$d4$db_positions4\]\)\n";
						print OutFile1 "$LipidClass\(1\'\-\[$c3\:$d3$db_positions3\/$c1\:$d1$db_positions1\]\, 3\'\-\[$c4\:$d4$db_positions4\/$c2\:$d2$db_positions2\]\)\n";
						print OutFile1 "$LipidClass\(1\'\-\[$c3\:$d3$db_positions3\/$c2\:$d2$db_positions2\]\, 3\'\-\[$c1\:$d1$db_positions1\/$c4\:$d4$db_positions4\]\)\n";
						print OutFile1 "$LipidClass\(1\'\-\[$c3\:$d3$db_positions3\/$c2\:$d2$db_positions2\]\, 3\'\-\[$c4\:$d4$db_positions4\/$c1\:$d1$db_positions1\]\)\n";
						print OutFile1 "$LipidClass\(1\'\-\[$c3\:$d3$db_positions3\/$c4\:$d4$db_positions4\]\, 3\'\-\[$c1\:$d1$db_positions1\/$c2\:$d2$db_positions2\]\)\n";
						print OutFile1 "$LipidClass\(1\'\-\[$c3\:$d3$db_positions3\/$c4\:$d4$db_positions4\]\, 3\'\-\[$c2\:$d2$db_positions2\/$c1\:$d1$db_positions1\]\)\n";

						print OutFile1 "$LipidClass\(1\'\-\[$c4\:$d4$db_positions4\/$c1\:$d1$db_positions1\]\, 3\'\-\[$c2\:$d2$db_positions2\/$c3\:$d3$db_positions3\]\)\n";
						print OutFile1 "$LipidClass\(1\'\-\[$c4\:$d4$db_positions4\/$c1\:$d1$db_positions1\]\, 3\'\-\[$c3\:$d3$db_positions3\/$c2\:$d2$db_positions2\]\)\n";
						print OutFile1 "$LipidClass\(1\'\-\[$c4\:$d4$db_positions4\/$c2\:$d2$db_positions2\]\, 3\'\-\[$c1\:$d1$db_positions1\/$c3\:$d3$db_positions3\]\)\n";
						print OutFile1 "$LipidClass\(1\'\-\[$c4\:$d4$db_positions4\/$c2\:$d2$db_positions2\]\, 3\'\-\[$c3\:$d3$db_positions3\/$c1\:$d1$db_positions1\]\)\n";
						print OutFile1 "$LipidClass\(1\'\-\[$c4\:$d4$db_positions4\/$c3\:$d3$db_positions3\]\, 3\'\-\[$c1\:$d1$db_positions1\/$c2\:$d2$db_positions2\]\)\n";
						print OutFile1 "$LipidClass\(1\'\-\[$c4\:$d4$db_positions4\/$c3\:$d3$db_positions3\]\, 3\'\-\[$c2\:$d2$db_positions2\/$c1\:$d1$db_positions1\]\)\n";
					}
				}
			}
		}
	}
}	
	
sub Predict_SN_and_draw_Cer {

	my $total_OH_in_cer; # Total number of hydroxylations in the ceramide molecule. Andrej Drosophila Lipidome data uses this format
	# Breaking lipid name to extract fatty acid chain information
	my @LipidNameSplit;
	if ($LipidName =~ m/\s/){ # Any lipid except sterols
		@LipidNameSplit = split /\s/ , "$LipidName";
		@Chains = split /\// , $LipidNameSplit[1];
	}
	if (($LipidClass =~ m/Cer/) and defined($LipidNameSplit[2])){ # Applicable for Howard's yeastlipidome data. This parsing rule minimises complexity in nomencature. Not recommened 
	# print "\n \$LipidNameSplit\[2\] = $LipidNameSplit[2]\n"; 
		if ($LipidNameSplit[2] == 1){$h2 = 1;}
		elsif ($LipidNameSplit[2] == 2){$h2 = 2;}
		else {print "\nUnknown text found in Lipid name\n" and die "Program Terminating !!";}
	}

	# Splitting string "chain_1" to obtain 1) No. of carbon atoms 2) No. of double bonds
	@_ = split("\:|\;", $Chains[0]);
	($c_atoms,$d_bonds) = ($_[0], $_[1]);
	
	if (defined ($_[2])) { # Applicable for Andrej's Drosophila lpidome data ONLY	

		$total_OH_in_cer = $_[2]; # Total number of hydroxylations in the ceramide molecule. Andrej Drosophila Lipidome data uses this format

		# For Cer/CerPE/HexCer only... 
		if ($LipidClass =~ m/(Cer|PECer|HexCer)/){ # print "Writing possiblities for Cer/CerPE...";<stdin>;
			$LipidClass =~ s/HexCer/GalCer/; # Giving Galactosyl Ceramide as one possibility for HexCer

			foreach (@CerList){
				$c1 = $_; # print "\n\$c1 = $c1..\n";
				$d1 = "0";
				$db_positions1 = "";
				foreach (@FAlist){
					($c2, $d2, $db_positions2) = ($_->[0], $_->[1], $_->[2]); # NoC1 = No. of Carbon atoms in FA 1
					if (defined($db_positions2)) {} else {$db_positions2="";}

					if ((($c1+$c2)== $c_atoms) and ($d1+$d2 == $d_bonds) and ($c1!=0) and ($c2!=0))  {

						if ($total_OH_in_cer == 2) { 
							# If there are two hydroxylations, then two possibilities are taken into consideration
							# 1. Both in head group (3,4)
								print OutFile1 "$LipidClass\(d$c1\:$d1$db_positions1\/$c2\:$d2$db_positions2\)\n";
								# Giving Glucosyl Ceramide as second possibility for HexCer
								if ($LipidClass =~ m/GalCer/){
									print OutFile1 "GlcCer\(d$c1\:$d1$db_positions1\/$c2\:$d2$db_positions2\)\n";
								}
							# 2. One in head group (3); second in Fatty acid group (at 2 nd carbon atom from COOH end)
								print OutFile1 "$LipidClass\(m$c1\:$d1$db_positions1\/$c2\:$d2$db_positions2\(2OH\)\)\n";
								# Giving Glucosyl Ceramide as second possibility for HexCer
								if ($LipidClass =~ m/GalCer/){
									print OutFile1 "GlcCer\(m$c1\:$d1$db_positions1\/$c2\:$d2$db_positions2\(2OH\)\)\n";
								}
						}
						elsif ($total_OH_in_cer == 3) {
							# If there are two hydroxylations, then the positions of these hydroxylations are fixed
							# two in Head group ; one in F'Acid chain
								print OutFile1 "$LipidClass\(d$c1\:$d1$db_positions1\/$c2\:$d2$db_positions2\(2OH\)\)\n";
								# Giving Glucosyl Ceramide as second possibility for HexCer
								if ($LipidClass =~ m/GalCer/){
									print OutFile1 "GlcCer\(d$c1\:$d1$db_positions1\/$c2\:$d2$db_positions2\(2OH\)\)\n";
								}
						}
					}

				}
			}
		}
		else {
			print "\nUnknown LipidClass $LipidClass in molecule $LipidName" and die "Program Terminating";
		}

	}

	else { # print "I am at line 965\n";
		foreach (@CerList){
			$c1 = $_; # print "\n\$c1 = $c1..\n";
			$d1 = "0";
			$db_positions1 = "";
			foreach (@FAlist){
				($c2, $d2, $db_positions2) = ($_->[0], $_->[1], $_->[2]); # NoC1 = No. of Carbon atoms in FA 1
				if (defined($db_positions2)) {} else {$db_positions2="";}

				if ((($c1+$c2)== $c_atoms) and ($d1+$d2 == $d_bonds) and ($c1!=0) and ($c2!=0))  {
					if ($OH_Flag == 0){ # $OH_Flag here is similar to $h1. $OH_Flag and $h1 represent no. of hydroxylations in LongChainBase
						print OutFile1 "$LipidClass\(m$c1\:$d1$db_positions1\/$c2\:$d2$db_positions2$h_position\)\n";
					}
					elsif ($OH_Flag==1){ print OutFile1 "$LipidClass\(d$c1\:$d1$db_positions1\/$c2\:$d2$db_positions2$h_position\)\n";}
					elsif ($OH_Flag==2){ print OutFile1 "$LipidClass\(t$c1\:$d1$db_positions1\/$c2\:$d2$db_positions2$h_position\)\n";}
					elsif ($OH_Flag==3){ 
						print OutFile1 "$LipidClass\(t$c1\:$d1$db_positions1\/$c2\:$d2$db_positions2\)\n";
						print OutFile1 "$LipidClass\(d$c1\:$d1$db_positions1\/$c2\:$d2$db_positions2$h_position\)\n";
					}
				}
			}
		}
	}
}

sub DelDuplicates {
	open FILEA, "./$InputFolderName/OutFile1.txt" or die "\n$! Unable to open OutFile1.txt\n";
	my @list3=<FILEA>;
	close FILEA;
	@list3 = uniq(@list3);
	open FILEB, "> ./$InputFolderName/OutFile1.txt";
	print FILEB @list3;
	close FILEB;
}

sub FindCentroidIsomer {
	# Generating MOL-MOL file for all possible structures using LipidMaps Structure Drawing Tools
	if ($LipidClass =~ m/(TG|DG)/){ $LipidMapsProgramToUse = "GLStrGen.pl";}
	elsif ($LipidClass =~ m/(Cer|PICer|PECer|HexCer|MIPCer|MIP2Cer)/){ $LipidMapsProgramToUse = "SPStrGen.pl";}
	elsif ($LipidClass =~ m/(PC|PA|PS|PI|PE|PG)/){ $LipidMapsProgramToUse = "GPStrGen.pl";}
	elsif ($LipidClass =~ m/CL/){ $LipidMapsProgramToUse = "CLStrGen.pl";} # print "\n\$LipidMapsProgramToUse = $LipidMapsProgramToUse"; <stdin>;
	else {print "\nI could not figure out the LipidMaps Structure Drawing program to use\n" and die "Program Terminating\n";}
	
	`perl ./bin/*_lipidmapstools/bin/$LipidMapsProgramToUse -m AbbrevFileName ./$InputFolderName/OutFile1.txt -r ./$InputFolderName/OutFile1`;
	# `mv OutFile1.sdf ./$InputFolderName/OutFile1.sdf`;
	
	# Converting MOL-MOL structures to SMILES using "OpenBabel : Molecule Conversion package"
	`babel -isdf ./$InputFolderName/OutFile1.sdf -osmi ./$InputFolderName/OutFile1.smi -d -b `;  # print "check 01.smi file";<stdin>;

	# Finding Centroid structure
		# Changing molecule names (because of limitation in ASCII characters that can be used as Molecule names

		open OutFile3, "./$InputFolderName/OutFile1.smi";
		open O2, "> ./$InputFolderName/O2.smi";
		$i=0;
		@VarList1=(); # Array to store SMILES strings
		@VarList2=(); # Array to store Molecule names
		while (<OutFile3>){
			($VarString1, $VarString2) = split("\t", $_);
			$VarString1 = &ModifySMILES($VarString1);
			push (@VarList1, "$VarString1");
			chomp($VarString2);
			push (@VarList2, "$VarString2");
			print O2 "$VarString1 $i\n";
			$i++;
		}
		close O2;  # print "check 02.smi file";<stdin>;
		close OutFile3;	

	# If there are no possible sructures..
	
	if ($i == 0){ 
		print STDERR "\nUnable to generate structures from Fatty Acid possibilities defined in ./$InputFolderName/FA file\n" ; 
		push @UnsupportedSpecies, $LipidName;
#		`rm ./$InputFolderName/OutFile1.sdf ./$InputFolderName/OutFile1.smi ./$InputFolderName/OutFile1.txt ./$InputFolderName/O2.smi`;
	}
	
	# If there is only one possible structure for a given lipid species, take that structure as default centroid structure
	elsif ($i == 1) {$centroid = 0;} 
	
	# If there is only TWO possible structure for a given lipid species, choose first one as centroid structure
	elsif ($i == 2) {$centroid = 0;} 
	# If there is only TWO possible structure for a given lipid species, randomly choose one of them as centroid structure
	# elsif ($i == 2) {$centroid = int(rand(2));} 

	# If there are more than TWO possible structures, find centroid  structure using any of the six similarity scoring methods
	else {
	
		if    ($MethodChoice==1){`perl		./bin/*m1.pl	./$InputFolderName/O2.smi					./$InputFolderName/O2.csv	./$InputFolderName/O2.log`;}
		elsif ($MethodChoice==2){`python	./bin/*m2.py	./$InputFolderName/O2.smi					./$InputFolderName/O2.csv	./$InputFolderName/O2.log`;}
		#elsif ($MethodChoice==3){`perl	./bin/*m3.pl	./$InputFolderName/O2.smi					./$InputFolderName/O2.csv	./$InputFolderName/O2.log`;}
		elsif ($MethodChoice==4){`python	./bin/*m4.py	./$InputFolderName/O2.smi					./$InputFolderName/O2.csv	./$InputFolderName/O2.log`;}
		elsif ($MethodChoice==5){`python	./bin/*m5.py	./$InputFolderName/O2.smi	./bin/mat4		./$InputFolderName/O2.csv	./$InputFolderName/O2.log`;}
		elsif ($MethodChoice==6){`perl		./bin/*m6.pl	./$InputFolderName/O2.smi	./$InputFolderName/O2.aln	./$InputFolderName/O2.csv	./$InputFolderName/O2.log`;
						`rm ./$InputFolderName/O2.aln ./$InputFolderName/O2.aln_joined ./$InputFolderName/O2.aln_joined_rearranged`;}

		# print "check csv file";<stdin>;
		# input *.csv file into centroid.pl program
		# get the center molecule and it's SMILES
		# $centroid =`echo "./$InputFolderName/O2.csv" | perl ./bin/centroid2.pl;`; print "$centroid";<stdin>;
		$centroid =`perl ./bin/centroid2.pl ./$InputFolderName/O2.csv`; # print "\$centroid = $centroid";

		# Deleting temporary files..
		`rm ./$InputFolderName/O2.csv`;
	
	}

	# Deleting temporary files..
	`rm ./$InputFolderName/OutFile1.sdf ./$InputFolderName/OutFile1.smi ./$InputFolderName/OutFile1.txt ./$InputFolderName/O2.smi`;

	if ($i != 0) {&Add2Lists;}
}

sub Add2Lists {

	# To CentroidList
	$LipidName =~ s/(\s|;|:|,)/_/g;
	$_ = scalar(@VarList2);
	push (@CentroidList, "$LipidName\t$VarList2[$centroid]\t$_"); 
	print OutFile2 "$LipidName\t";
	print OutFile2 join ("\t", @VarList2); 
	print OutFile2 "\n";
	
	# To SMILESlist
	push (@LipidNamesList, "$VarList1[$centroid] $LipidName"); 
}

sub WriteLists2Files{

	# Write CentroidList to File
	$CentroidListFileName = "CentroidList\.txt";
	open O3, "> ./$InputFolderName/$CentroidListFileName" or die "Unable to open ./$InputFolderName/$CentroidListFileName file";
	print O3 join ("\n", @CentroidList); 
	close O3;

	# Write SMILESlist to File
	$SMILESlistFileName = "SMILESlist\.smi";
	open O4, "> ./$InputFolderName/$SMILESlistFileName";
	print O4 join ("\n", @LipidNamesList); 
	close O4; 

	# Write ClassList to File
	@ClassList = uniq(@ClassList); # removing duplicate classes
	open O5, "> ./$InputFolderName/ClassList.txt";
	print O5 join ("\n", @ClassList);
	print O5 "\n";
	close O5;
}

__END__
