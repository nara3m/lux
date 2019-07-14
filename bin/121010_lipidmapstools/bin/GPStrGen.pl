#!/usr/bin/perl -w
#
# File: GPStrGen.pl
# Author: Manish Sud
# Contributor(s): Eoin Fahy
#
# Copyright (C) 2006-2012. The Regents of the University of California. All Rights Reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
# 3. Neither the name of the LIPID MAPS (TM) consortium nor the names of its
#    contributors may be used to endorse or promote products derived from this
#    software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
# ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# LIPID MAPS (TM) Tools software package is developed by Manish Sud and Eoin Fahy
# at San Diego Supercomputer Center at the University of California, San Diego as
#  a part of LIPID Metabolites and Pathways Strategy consortium efforts.
#
# The LIPID MAPS (TM) consortium is a multi-institutional multi-year effort
# involved in identification and characterization of existing and novel lipids,
# quantifications of changes in their metabolites, and development of biochemical
# pathways and interaction network maps involving lipids. Please acknowledge use
# of the LIPID MAPS (TM) Tools in your work by citing appropriate literature
# references. For further details, please visit www.lipidmaps.org.
#
# Support for this effort is provided by  by National Institutes of Health (NIH)
# and National Institute of General Medical Sciences (NIGMS) Glue Grant
# NIH/NIGMS Grant 1 U54 GM69338.
#

use strict;
use FindBin; use lib "$FindBin::Bin/../lib";
use Getopt::Long;
use File::Basename;
use Text::ParseWords;
use Benchmark;
use FileUtil;
use LMAPSStr;
use GPStr;

my($ScriptName, %Options, $StartTime, $EndTime, $TotalTime);

# Autoflush STDOUT
$| = 1;

# Starting message...
$ScriptName = basename($0);
print "\n$ScriptName: Starting...\n\n";
$StartTime = new Benchmark;

# Get the options and setup script...
SetupScriptUsage();
if ($Options{help} || @ARGV < 1) {
  die GetUsageFromPod("$FindBin::Bin/$ScriptName");
}

print "Processing options...\n";
my($SDFileName, $SpecifiedCmpdAbbrevsRef, $WriteSDFile, $AllowArbitraryChainAbbrev);
ProcessOptions();

print "Processing abbreviation(s)...\n";
GPStr::ProcessGPCmpdAbbrevs($SpecifiedCmpdAbbrevsRef, $AllowArbitraryChainAbbrev, $WriteSDFile, $SDFileName);

print "$ScriptName:Done...\n\n";

$EndTime = new Benchmark;
$TotalTime = timediff ($EndTime, $StartTime);
print "Total time: ", timestr($TotalTime), "\n";

###############################################################################

# Process command line options...
sub ProcessOptions {
  # Setup SD file name
  $SDFileName = LMAPSStr::SetupSDFileName('GP', \%Options);
  if (!$Options{overwrite}) {
      if (-e $SDFileName) {
	die "Error: The file $SDFileName already exists\n";
      }
  }
  $SpecifiedCmpdAbbrevsRef = LMAPSStr::SetupCmpdAbbrevs(\%Options);

  $WriteSDFile = ($Options{processmode} =~ /^WriteSDFile$/i) ? 1 : 0;
  $AllowArbitraryChainAbbrev = ($Options{chainabbrevmode} =~ /^Arbitrary$/i) ? 1 : 0;
}

# Setup script usage  and retrieve command line arguments specified using various options...
sub SetupScriptUsage {

  # Retrieve all the options...
  %Options = ();
  $Options{chainabbrevmode} = 'MostLikely';
  $Options{mode} = 'Abbrev';
  $Options{processmode} = 'WriteSDFile';

  if (!GetOptions(\%Options, "chainabbrevmode|c=s", "help|h", "mode|m=s", "processmode|p=s", "overwrite|o", "root|r=s", "workingdir|w=s")) {
    die "\nTo get a list of valid options and their values, use \"$ScriptName -h\" or\n\"perl -S $ScriptName -h\" command and try again...\n";
  }
  if ($Options{workingdir}) {
    if (! -d $Options{workingdir}) {
      die "Error: The value specified, $Options{workingdir}, for option \"-w --workingdir\" is not a directory name.\n";
    }
    chdir $Options{workingdir} or die "Error: Couldn't chdir $Options{workingdir}: $! \n";
  }
  if ($Options{chainabbrevmode} !~ /^(MostLikely|Arbitrary)$/i) {
    die "Error: The value specified, $Options{chainabbrevmode}, for option \"-c, --ChainAbbrevMode\" is not valid. Allowed values: MostLikely or Arbitrary\n";
  }
  if ($Options{mode} !~ /^(Abbrev|AbbrevFileName)$/i) {
    die "Error: The value specified, $Options{mode}, for option \"-m --mode\" is not valid. Allowed values: Abbrev or AbbrevFileName\n";
  }
  if ($Options{processmode} !~ /^(WriteSDFile|CountOnly)$/i) {
    die "Error: The value specified, $Options{processmode}, for option \"-p, --ProcessMode\" is not valid. Allowed values: WriteSDFile or CountOnly\n";
  }
}


__END__

=head1 NAME

GPStrGen.pl - Generate structures for Glycerophospholipids (GP)

=head1 SYNOPSIS

GPStrGen.pl  GPAbbrev|GPAbbrevFileName ...

GPStrGen.pl [B<-c, --ChainAbbrevMode> I<MostLikely | Arbitrary>]
[B<-h, --help>] [B<-m, --mode> I<Abbrev | AbbrevFileName>]
[B<-p, --ProcessMode> I<WriteSDFile | CountOnly>] [B<-o, --overwrite>]
[B<-r, --root> rootname] [B<-w, --workingdir> dirname] <arguments>...

=head1 DESCRIPTION

Generate Glyceriphospholipids (GP) structures using compound abbreviations specified on
a command line or in a CSV/TSV Text file. All the command line arguments represent either
compound abbreviations or file name containing abbreviations. Use mode option to control
the type of command line arguments.

A SD file, containing structures for all GP abbreviations along with ontological information, is
generated as an output.

=head1 SUPPORTED ABBREVIATIONS

Current support for GP structure generation include these main classes and sub classes:

o Glycerophosphocholines (PC)

    . Diacylglycerophosphocholines
    . 1-alkyl,2-acylglycerophosphocholines
    . 1Z-alkenyl,2-acylglycerophosphocholines
    . Dialkylglycerophosphocholines
    . Monoacylglycerophosphocholines
    . 1-alkyl glycerophosphocholines
    . 1Z-alkenylglycerophosphocholines

o Glycerophosphoethanolamines (PE)

    . Diacylglycerophosphoethanolamines
    . 1-alkyl,2-acylglycerophosphoethanolamines
    . 1Z-alkenyl,2-acylglycerophosphoethanolamines
    . Dialkylglycerophosphoethanolamines
    . Monoacylglycerophosphoethanolamines
    . 1-alkyl glycerophosphoethanolamines
    . 1Z-alkenylglycerophosphoethanolamines

o Glycerophosphoserines (PS)

    . Diacylglycerophosphoserines
    . 1-alkyl,2-acylglycerophosphoserines
    . 1Z-alkenyl,2-acylglycerophosphoserines
    . Dialkylglycerophosphoserines
    . Monoacylglycerophosphoserines
    . 1-alkyl glycerophosphoserines
    . 1Z-alkenylglycerophosphoserines

o Glycerophosphoglycerols (PG)

    . Diacylglycerophosphoglycerols
    . 1-alkyl,2-acylglycerophosphoglycerols
    . 1Z-alkenyl,2-acylglycerophosphoglycerols
    . Dialkylglycerophosphoglycerols
    . Monoacylglycerophosphoglycerols
    . 1-alkyl glycerophosphoglycerols
    . 1Z-alkenylglycerophosphoglycerols

o Glycerophosphoglycerophosphates (PGP)

    . Diacylglycerophosphoglycerophosphates
    . 1-alkyl,2-acylglycerophosphoglycerophosphates
    . 1Z-alkenyl,2-acylglycerophosphoglycerophosphates
    . Dialkylglycerophosphoglycerophosphates
    . Monoacylglycerophosphoglycerophosphates
    . 1-alkyl glycerophosphoglycerophosphates
    . 1Z-alkenylglycerophosphoglycerophosphates

o Glycerophosphoinositols (PI)

    . Diacylglycerophosphoinositols
    . 1-alkyl,2-acylglycerophosphoinositols
    . 1Z-alkenyl,2-acylglycerophosphoinositols
    . Dialkylglycerophosphoinositols
    . Monoacylglycerophosphoinositols
    . 1-alkyl glycerophosphoinositols
    . 1Z-alkenylglycerophosphoinositols

o Glycerophosphoinositol monophosphates (PIP)

    . Diacylglycerophosphoinositol monophosphates
    . 1-alkyl,2-acylglycerophosphoinositol monophosphates
    . 1Z-alkenyl,2-acylglycerophosphoinositol monophosphates
    . Dialkylglycerophosphoinositol monophosphates
    . Monoacylglycerophosphoinositol monophosphates
    . 1-alkyl glycerophosphoinositol monophosphates
    . 1Z-alkenylglycerophosphoinositol monophosphates

o Glycerophosphates (PA)

    . Diacylglycerophosphates
    . 1-alkyl,2-acylglycerophosphates
    . 1Z-alkenyl,2-acylglycerophosphates
    . Dialkylglycerophosphates
    . Monoacylglycerophosphates
    . 1-alkyl glycerophosphates
    . 1Z-alkenylglycerophosphates

o Glyceropyrophosphates (PPA)

    . Diacylglyceropyrophosphates
    . Monoacylglyceropyrophosphates

o Glycerophosphonocholines (PnC)

    . Diacylglycerophosphonocholines
    . 1-alkyl,2-acylglycerophosphonocholines
    . 1Z-alkenyl,2-acylglycerophosphonocholines
    . Dialkylglycerophosphonocholines
    . Monoacylglycerophosphonocholines
    . 1-alkyl glycerophosphonocholines
    . 1Z-alkenylglycerophosphonocholines

o Glycerophosphonoethanolamines (PnE)

    . Diacylglycerophosphonoethanolamines
    . 1-alkyl,2-acylglycerophosphonoethanolamines
    . 1Z-alkenyl,2-acylglycerophosphonoethanolamines
    . Dialkylglycerophosphonoethanolamines
    . Monoacylglycerophosphonoethanolamines
    . 1-alkyl glycerophosphonoethanolamines
    . 1Z-alkenylglycerophosphonoethanolamines


=head1 OPTIONS

=over 4

=item B<-c, --ChainAbbrevMode> I<MostLikely|Arbitrary>

Specify what types of acyl chain abbreviations are allowed during processing of complete
abbreviations: allow most likely chain abbreviations containing specific double bond geometry
specifications; allow any acyl chain abbreviation with valid chain length and double bond
geometry specificatios. Possible values: I<MostLikely or Arbitrary>. Default value: I<MostLikely>.

I<Arbitrary> value of B<-c, --ChainAbbrevMode> option is not allowed during processing of
abbreviations containing wild cards.

During I<MostLikely> value of B<-c, --ChainAbbrevMode> option, only the most likely acyl chain
abbreviations specified in ChainAbbrev.pm module are allowed. However, during I<Arbitrary> value
of B<-c, --ChainAbbrevMode> option, any acyl chain abbreviations with valid chain length and
double bond geometry can be specified. The current release of lipidmapstools support chain
lengths from 2 to 50 as specified in ChainAbbev.pm module.

In addition to double bond geometry specifications, valid substituents can be specified for in the acyl
chain abbreviations.

=item B<-h, --help>

Print this help message

=item B<-m, --mode> I<Abbrev|AbbrevFileName>

Controls interpretation of command line arguments. Two different methods are provided:
specify compound abbreviations or a file name containing compound abbreviations. Possible
values: I<Abbrev or AbbrevFileName>. Default: I<Abbrev>

In I<AbbrevFileName> mode, a single line in CSV/TSV files can contain multiple compound
abbreviations. The file extension determines delimiter used to process data lines: comma for
CSV and tab for TSV. For files with TXT extension, only one compound abbreviation per line
is allowed.

Wild card character, *, is also supported in compound abbreviations.

Examples:

    Specific structures: PC(12:0/13:0) PC(17:1(9Z)/0:0)
                         PA(13:0/0:0)
    Specific structures: PC(O-16:0/13:0) PC(P-16:0/0:0)
    Specific possibilities: PC(21:0/22:*) PA(17:*/0:0)
                            PE(O-18:0/*:*)
    All possibilites: *(*:*/*:*) or *(*/*)

With wild card character, +/- can also be used for chain lengths to indicate even and odd lengths at
sn1/sn2/sn3 positions; additionally > and < qualifiers are also allowed to specify length
requirements. Examples:

    Odd and even number chains at sn1 and sn2: *(*+:*/*-:*)
    Odd and even number chains at sn1 and sn2 with length longer than 10
       and 20: *(*+>10:*/*->20:*)

Default sn2 stereochemistry is R. However, abbreviation format also supports these additional stereochemistry
specifications for sn2 position: S; U - unknown; rac - racemic mixture. Examples:

    PC(12:0/13:0)[rac]
    PC(17:1(9Z)/14:0)[S]
    PA(13:0/12:0)[U]

=item B<-p, --ProcessMode> I<WriteSDFile|CountOnly>

Specify how abbreviations are processed: generate structures for specified abbreviations along
with generating a SD file or just count the number of structures corresponding to specified
abbreviations without generating any SD file. Possible values: I<WriteSDFile or CountOnly>.
Default: I<WriteSDFile>.

It can take substantial amount of time for generating all the structures and writing out a SD file
for abbreviations containing wild cards. I<CountOnly> value of B<--ProcessMode> option can
be used to get a quick count of number of structures to be generated without writing out any
SD file.

=item B<-o, --overwrite>

Overwrite existing files

=item B<-r, --root> I<rootname>

New file name is generated using the root: <Root>.sdf. Default for new file names: GPAbbrev.sdf,
<AbbrevFilenName>.sdf, or <FirstAbbrevFileName>1To<Count>.sdf.

=item B<-w, --workingdir> I<dirname>

Location of working directory. Default: current directory

=back

=head1 EXAMPLES

On some systems, command line scripts may need to be invoked using
I<perl -s GLStrGen.pl>; however, all the examples assume direct invocation
of command line script works.

To generate a GPStructures.sdf file containing a structure specified
by a command line GP abbreviation, type:

    % GPStrGen.pl -r GPStructures -o "PC(16:0/0:0)"

To generate a GPStructures.sdf file containing structures specified
by a command line GL abbreviations, type:

    % GPStrGen.pl -r GPStructures -o "PC(16:0/0:0)" "PE(18:1(11E)/16:0)"

To generate a GPStructures.sdf file containing structures specified
by a command line GP abbreviations with specific stereochemistry, type:

    % GPStrGen.pl -r GPStructures -o "PC(16:0/0:0)[U]"
    "PE(18:1(11E)/16:0)[S]"

To enumerate all possible GP structures and generate a GPStructures.sdf
file, type:

    % GPStrGen.pl -r GPStructures -o "*(*/*)"

or

    % GPStrGen.pl -r GPStructures -o "*(*:*/*:*)"

or

    % GPStrGen.pl -r GPStructures -o "*(*:*(*)/*:*(*))"

To enumerate all possible GP structures with a sn1 chain, and generate a
GPStructures.sdf file, type:

    % GPStrGen.pl -r GPStructures -o "*(*/0:0)"


To enumerate all possible GP structures with a sn1 chain containing one
double bond, and generate a GPStructures.sdf file, type:

    % GPStrGen.pl -r GPStructures -o "*(*:1/0:0)"

To enumerate all possible GP structures with even chain length larger than
10 at sn1 position, and generate and generate a GPStructures.sdf file, type:

    % GPStrGen.pl -r GPStructures -o "*(*+>10:*/0:0)"

To enumerate all possible GP structures with odd chains longer
than 10 at sn1 and even chains longer than 18 at sn2, and generate a
GPStructures.sdf file, type:

    % GPStrGen.pl -r GPStructures -o "*(*->10:*/*+>18:*)"

=head1 AUTHOR

Manish Sud

=head1 CONTRIBUTOR

Eoin Fahy

=head1 SEE ALSO

CLStrGen.pl, FAStrGen.pl, GLStrGen.pl, SPStrGen.pl, STStrGen.pl

=head1 COPYRIGHT

Copyright (C) 2006-2012. The Regents of the University of California. All Rights Reserved.

=head1 LICENSE

Modified BSD License

=cut
