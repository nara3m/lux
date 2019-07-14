-------------------------------------------------------------
Before you can start working with the scripts in this folder, download and install following software and libraries. All scripts were designed-for and tested on ubuntu (linux) operating system

1. R software from http://www.r-project.org/
2. 'graphics', 'gridExtra', 'gplots', 'RColorBrewer', 'scatterplot3d', 'gdata', 'RSVGTipsDevice' libraries for R from http://cran.r-project.org/web/packages/gdata/index.html
5. Latest version of perl software from http://www.perl.org/
6. Latest version of python software from https://www.python.org/
7. Latest version of openbabel software from http://openbabel.org

-------------------------------------------------------------

To Analyze a set of lipidomes, follow these 5 steps:

Step 1: place all you lipidomes in a folder. Two example lipidomes are provided with this package. Folder "DM" contains data from Drosophila lipidome (Carvalho et. al. Mol. Sys. Bio., Vol 1, Issue 8, 2012) and folder "SC" contains yeast lipidome (Ejsing et. al., PNAS, Vol 106, Issue 7, 2009). Follow organization of your data inside the folder similar to the example files.

As you can see in the example folder "SC" there are a total of 11 files with the following names :
(example folder "DM" contains 15 files)

1. BY4741_24
2. BY4741_37
3. Elo1_24
4. Elo1_37
5. Elo2_24
6. Elo2_37
7. Elo3_24
8. Elo3_37
9. FA.tab
10.UserSMILES.tab
11. Cer.tab

Files 1-8 are real lipidome data (Esjing et. al 2010 PNAS) included here for illustration. 
9, 10 and 11 are special files that have to included in all new folders you create.

Files 1-8 follow a tab delimited format. First column is the name of the lipid. Second column is the concentration of that lipid. If you do not have concentrations, put '0's in concentration column - DO NOT LEAVE THAT COLOUMN EMPTY.

Step 2: Check Lipid Names
Naming of each lipid follows a format. More than 20 classes of lipids are used in example lipidome datasets (files 1-8). Use those examples if you have to name a new lipid. Spacing and punctuation marks are important in writing lipid names.

Step 3: Check input files for formatting errors
Make sure there are no empty rows in any of the input files. Our programs do not check for empty rows/columns. In our experience, having empty rows/columns will lead to incorrect results.

Step 4: Check contents of Mandatory input files.
Contents of File 9 "FA.tab" will be used to draw fatty acid structures. In any given organism, including yeast, not all possible fatty acids are synthesized. User has to define what fatty acid possibilities must be used. Fatty Acid combinations we used for yeast data in the manuscript are provided in this file. Similarly, fatty acids used for drosophila data are provided in the folder "DM" under the file name "FA.tab". User is free to add/remove some of them as per their requirement. But note that if you are including some "unusual" fatty acids, LIPIDMAPS structure drawing programs also need to modified. User may contact LipidMaps team (www.lipidmaps.org) or us to know how to edit LIPIDMAPS structure drawing programs.

NOTE: In our case we specifically changed the program "121010_lipidmapstools\lib\ChainAbbrev.pm" to enable drawing specific yeast carbon chains. Modifications are clearly marked in the program, which can help you to recognize how you might customize LIPIDMAPS programs for your lipidome of choice.

Contents of File 11 "Cer.txt" will be used to assign the number of carbon atoms in long chain base (LCB) of Ceramide molecules. Notice the contents of this file in folders "SC" and "DM" to see the differences in Ceramide composition of yeast and Drosophila.

File number 10, with name "UserSMILES.tab" contains User_defined_SMILES for lipid species. Our program uses LIPIDMAPS structure drawing tools to generate structures. LIPIDMAPS tools do not support all classes of lipids.  It is highly possible that some of the lipids in your dataset/organism are NOT supported by LIPIDMAPS. In such cases, user can define their own SMILES. If you do not know which lipids are not supported, just run this program on your dataset(s) with default options. Lipids for which structures were not drawn will be written to file "Unused_LipidAbbreviations.txt". After adding SMILES, user can re-run the program.

make sure "bin" folder (provided with this package) is located in the same directory as this README file.

Last but not the least, there should be no other files in the folder. Our program treats all files in the folder (excluding files named "FA.tab", "Cer.txt" and "UserSMILES.tab") to be lipidome data.

Step 5: Run the program with following syntax

perl 150629_AnalyzeLipidomes.pl <folder_name> 

example 1: > perl 150629_AnalyzeLipidomes.pl DM
example 2: > perl 150629_AnalyzeLipidomes.pl SC

If there are any dependency issues, you will be notified through STDERR or STDOUT.

This program generates many result files, all of which will be located in the folder you have provided as input. In the example shown above, all results will be placed in the folder "SC" (or "DM"). Taking result files and "SC" folder as example, we will describe here what to expect in the output. 

Our program generates the following (24) files with above example "SC" (excluding 11 input files)

Set 1: Main Results - Shown in Figures of Manuscript
01. SC_2D.svg: Interactive PCA plot for first two principal compoents. For all lipidomes combined.
02. SC_3D.svg: Interactive PCA plot for first three principal compoents. Again, for all lipidomes combined.
03. Dendrograms.pdf: Dendrograms and heatmaps shown in Panels A, B and C in Fig6 of manuscript. 

Set 2: Main Data - Data used for generating figures in manuscript
04. SC.tab: Merged data from all lipidomes
05. SC.pca: Coordinates of first three principal components for all lipid species.
06. SC.csv: Pairwise distances between all lipid species using Levenshtein distance.
07. SC_VariancePlot.pdf: Distribution of variance in the first 10 principal components.
08. LUX_PC1_PC2.csv: Pairwise LUX Scores derived using PC1 and PC2 coordinates.
09. LUX_PC1_PC2_PC3.csv: Pairwise LUX Scores derived using PC1, PC2 and PC3 coordinates. 
10. LUX_SimilarityScore.csv: Pairwise LUX Scores derived using structural similarity values.

Set 3: Auxillary data related to Set 2 results
07. ClassList.txt: This file containts list of lipid class (abbreviated name) present in input dataset. 
08. CentroidList.txt: This file contains three columns. First column is the lipid species (input). Second column is the representative isomer for that lipid. Third column is the number of isomers considered in choosing representative isomer of second column.
09. SMILESlist.smi: List of Non Canonical SMILES for the representative isomer for each lipid.
10. LipidList.txt: List of all lipids in the given folder "SC"
11. Isomers.txt: List of all isomers for all lipids. Refer to manuscript and supplimentary methods for the description of isomers.
12. Unused_LipidAbbreviations.txt

Set 4: Data for further exploration: Interactive PCA plots for first two principal compoents for each lipidome. 
17. BY4741_24_1.svg
18. BY4741_37_1.svg
19. Elo1_24_1.svg
20. Elo1_37_1.svg
21. Elo2_24_1.svg
22. Elo2_37_1.svg
23. Elo3_24_1.svg
24. Elo3_37_1.svg