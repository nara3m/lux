#!/usr/bin/perl -w
#
# File: FAStrGen.pl
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
use FAStr;

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
my($SDFileName, $SpecifiedCmpdAbbrevsRef, $WriteSDFile);
ProcessOptions();

print "Processing abbreviation(s)...\n";
FAStr::ProcessFACmpdAbbrevs($SpecifiedCmpdAbbrevsRef, $WriteSDFile, $SDFileName);

print "\n$ScriptName:Done...\n\n";

$EndTime = new Benchmark;
$TotalTime = timediff ($EndTime, $StartTime);
print "Total time: ", timestr($TotalTime), "\n";

###############################################################################

# Process command line options...
sub ProcessOptions {
  # Setup SD file name
  $SDFileName = LMAPSStr::SetupSDFileName('FA', \%Options);
  if (!$Options{overwrite}) {
      if (-e $SDFileName) {
	die "Error: The file $SDFileName already exists\n";
      }
  }
  $SpecifiedCmpdAbbrevsRef = LMAPSStr::SetupCmpdAbbrevs(\%Options);

  $WriteSDFile = ($Options{processmode} =~ /^WriteSDFile$/i) ? 1 : 0;
}

# Setup script usage  and retrieve command line arguments specified using various options...
sub SetupScriptUsage {

  # Retrieve all the options...
  %Options = ();
  $Options{mode} = 'Abbrev';
  $Options{processmode} = 'WriteSDFile';

  if ($Options{mode} !~ /^(Abbrev|AbbrevFileName)$/i) {
    die "Error: The value specified, $Options{mode}, for option \"-m --mode\" is not valid. Allowed values: Abbrev or AbbrevFileName\n";
  }
  if (!GetOptions(\%Options, "help|h", "mode|m=s", "processmode|p=s", "overwrite|o", "root|r=s", "workingdir|w=s")) {
    die "\nTo get a list of valid options and their values, use \"$ScriptName -h\" or\n\"perl -S $ScriptName -h\" command and try again...\n";
  }
  if ($Options{workingdir}) {
    if (! -d $Options{workingdir}) {
      die "Error: The value specified, $Options{workingdir}, for option \"-w --workingdir\" is not a directory name.\n";
    }
    chdir $Options{workingdir} or die "Error: Couldn't chdir $Options{workingdir}: $! \n";
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

FAStrGen.pl - Generate structures for Fatty Acyls (FA)

=head1 SYNOPSIS

FAStrGen.pl  FAAbbrev|FAAbbrevFileName ...

FAStrGen.pl [B<-h, --help>] [B<-m, --mode> I<Abbrev | AbbrevFileName>]
[B<-p, --ProcessMode> I<WriteSDFile | CountOnly>] [B<-o, --overwrite>]
[B<-r, --root> rootname] [B<-w, --workingdir> dirname] <arguments>...

=head1 DESCRIPTION

Generate Fatty Acyls (FA) structures using compound abbreviations specified on
a command line or in a CSV/TSV Text file. All the command line arguments represent either
compound abbreviations or file name containing abbreviations. Use mode option to control
the type of command line arguments.

A SD file, containing structures for all SP abbreviations along with ontological information, is
generated as an output.

=head1 SUPPORTED ABBREVIATIONS

Current support for FA structure generation include these main classes and sub classes:

o Fatty Acids and Conjugates

    . Straight chain fatty acids
    . Methyl branched fatty acids
    . Unsaturated fatty acids
    . Hydroperoxy fatty acids
    . Hydroxy fatty acids
    . Oxo fatty acids
    . Epoxy fatty acids
    . Methoxy fatty acids
    . Halogenated fatty acids
    . Amino fatty acids
    . Cyano fatty acids
    . Nitro fatty acids
    . Thia fatty acids

o Eicosanoids

    . Prostaglandins

o Fatty alcohols

o Fatty aldehydes

o Fatty  amides

    . Primary amides

=head1 OPTIONS

=over 4

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

Wild card character, *, is also supported in compound abbreviations to generate straight
chain and unsaturated fatty acids.

Examples:

    Specific structures: "18:0" "20:4(5Z,8Z,11Z,14Z)"
                         "28:1(12Z)(2Me,4Me,6Me)"
                         "18:3(6Z,9Z,11E)(13OOH[S])"
                         "18:2(9E,11E)(13OH)"
                         "18:1(10E)(9Ke,10Ep)"
                         "16:1(5Z)(2OMe)" "7:1(2Z)(3Br)"
                         "18:2(9Z,12Z)(10NO2)" "16:2(10E,12Z)(1OH)"
                         "6:0(1CHO)" "12:0(1NH2)"
                         "20:2(5Z,13E)(9OH[S],11OH[R],15OH[S]){8a,12b}"
    All possibilites: *:* or *

With wild card character, +/- can also be used for chain lengths to indicate even and odd lengths;
additionally > and < qualifiers are also allowed to specify length requirements. Examples:

    Odd number chains: "*-:*"
    Even number chains: "*+:*"
    Odd number chains with chain length longer than 18: "*->18:*"
    Even number chains with chain length longer than 14: "*+>14:*"

=item B<-o, --overwrite>

Overwrite existing files

=item B<-r, --root> I<rootname>

New file name is generated using the root: <Root>.sdf. Default for new file names: FAAbbrev.sdf,
<AbbrevFilenName>.sdf, or <FirstAbbrevFileName>1To<Count>.sdf.

=item B<-p, --ProcessMode> I<WriteSDFile|CountOnly>

Specify how abbreviations are processed: generate structures for specified abbreviations along
with generating a SD file or just count the number of structures corresponding to specified
abbreviations without generating any SD file. Possible values: I<WriteSDFile or CountOnly>.
Default: I<WriteSDFile>.

It can take substantial amount of time for generating all the structures and writing out a SD file
for abbreviations containing wild cards. I<CountOnly> value of B<--ProcessMode> option can
be used to get a quick count of number of structures to be generated without writing out any
SD file.

=item B<-w, --workingdir> I<dirname>

Location of working directory. Default: current directory

=back

=head1 EXAMPLES

On some systems, command line scripts may need to be invoked using
I<perl -s FAStrGen.pl>; however, all the examples assume direct invocation
of command line script works.

To generate a FAStructures.sdf file containing a structure specified
by a command line FA abbreviation for straight chain fatty acids, type:

    % FAStrGen.pl -r FAStructures -o "18:0" "9:0"

To generate a FAStructures.sdf file containing a structure specified
by a command line FA abbreviation for methyl branched fatty acids, type:

    % FAStrGen.pl -r FAStructures -o "8:0(6Me)" "18:1(6Z)(17Me)"
    "28:1(12Z)(2Me,4Me,6Me)"

To generate a FAStructures.sdf file containing a structure specified
by a command line FA abbreviation for unsaturated fatty acids, type:

    % FAStrGen.pl -r FAStructures -o "20:4(5Z,8Z,11Z,14Z)" "8:1(5E)"

To generate a FAStructures.sdf file containing a structure specified
by a command line FA abbreviation for hydroperoxy fatty acids, type:

    % FAStrGen.pl -r FAStructures -o "18:2(9E,11E)(13OOH)"
    "18:3(6Z,9Z,11E)(13OOH[S])"

To generate a FAStructures.sdf file containing a structure specified
by a command line FA abbreviation for hydroxy fatty acids, type:

    % FAStrGen.pl -r FAStructures -o "10:0(10OH)" "15:0(2OH,15OH)"
     "18:2(9E,11E)(13OH)" "4:0(3OH[R])"

To generate a FAStructures.sdf file containing a structure specified
by a command line FA abbreviation for oxo fatty acids, type:

    % FAStrGen.pl -r FAStructures -o "10:0(2Ke)" "18:1(10E)(9Ke,10Ep)"

To generate a FAStructures.sdf file containing a structure specified
by a command line FA abbreviation for epoxy fatty acids, type:

    % FAStrGen.pl -r FAStructures -o "18:0(6Ep)"

To generate a FAStructures.sdf file containing a structure specified
by a command line FA abbreviation for methoxy fatty acids, type:

    % FAStrGen.pl -r FAStructures -o "18:1(9E)(12OH,13OH,11OMe)"
    "16:1(5Z)(2OMe)"

To generate a FAStructures.sdf file containing a structure specified
by a command line FA abbreviation for halogenated fatty acids, type:

    % FAStrGen.pl -r FAStructures -o "7:1(2Z)(3Br)" "26:2(5Z,9Z)(2Br)"

To generate a FAStructures.sdf file containing a structure specified
by a command line FA abbreviation for amino fatty acids, type:

    % FAStrGen.pl -r FAStructures -o "13:0(2NH2[S])" "4:0(2NH2,4CN)"

To generate a FAStructures.sdf file containing a structure specified
by a command line FA abbreviation for Cyano fatty acids, type:

    % FAStrGen.pl -r FAStructures -o "4:0(4CN)"

To generate a FAStructures.sdf file containing a structure specified
by a command line FA abbreviation for nitro fatty acids, type:

    % FAStrGen.pl -r FAStructures -o "18:2(9Z,12Z)(10NO2)"

To generate a FAStructures.sdf file containing a structure specified
by a command line FA abbreviation for prostaglanins, type:

    % FAStrGen.pl -r FAStructures -o "20:2(5Z,13E)(9OH[S],11OH[R],
    15OH[S]){8a,12b}"

To generate a FAStructures.sdf file containing a structure specified
by a command line FA abbreviation for fatty alcohols, type:

    % FAStrGen.pl -r FAStructures -o "26:0(1OH)" "16:2(10E,12Z)(1OH)"
    "11:0(1OH,2Me,2Me,9Me,9Me,10OH)"

To generate a FAStructures.sdf file containing a structure specified
by a command line FA abbreviation for fatty aldehydes, type:

    % FAStrGen.pl -r FAStructures -o "6:0(1CHO)" "16:2(2E,4E)(1CHO,6OH)"

To generate a FAStructures.sdf file containing a structure specified
by a command line FA abbreviation for primary amides, type:

    % FAStrGen.pl -r FAStructures -o "12:0(1NH2)"

To enumerate straight chain and unsaturated fatty acids with commonly occuring
chain lengths and generate FAStructures.sdf file, type:

    % FAStrGen.pl -r FAStructures -o "*"

or

    % FAStrGen.pl -r SPStructures -o "*:*"

=head1 AUTHOR

Manish Sud

=head1 CONTRIBUTOR

Eoin Fahy

=head1 SEE ALSO

CLStrGen.pl, GLStrGen.pl, GPStrGen.pl, SPStrGen.pl, STStrGen.pl

=head1 COPYRIGHT

Copyright (C) 2006-2012. The Regents of the University of California. All Rights Reserved.

=head1 LICENSE

Modified BSD License

=cut
