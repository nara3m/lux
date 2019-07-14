#!/usr/bin/perl -w
#
# File: LMStrGen.pl
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
use LMStr;

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
my($SDFileName, $SpecifiedCmpdAbbrevsRef);
ProcessOptions();

my($ExpandedCmdAbbrevsRef, $AbbrevCount);
$ExpandedCmdAbbrevsRef = LMStr::ExpandLMCmpdAbbrevs($SpecifiedCmpdAbbrevsRef);
$AbbrevCount = @$ExpandedCmdAbbrevsRef;

if ($AbbrevCount) {
  print "Generating SD file $SDFileName for $AbbrevCount compound", (($AbbrevCount > 1) ? "s" : "") ,"...\n";
  LMStr::GenerateSDFile($SDFileName, $ExpandedCmdAbbrevsRef);
}
else {
  print "No valid abbreviation found...\n";
}

print "$ScriptName:Done...\n\n";

$EndTime = new Benchmark;
$TotalTime = timediff ($EndTime, $StartTime);
print "Total time: ", timestr($TotalTime), "\n";

###############################################################################

# Process command line options...
sub ProcessOptions {
  # Setup SD file name
  $SDFileName = LMAPSStr::SetupSDFileName('LM', \%Options);
  if (!$Options{overwrite}) {
      if (-e $SDFileName) {
	die "Error: The file $SDFileName already exists\n";
      }
  }
  $SpecifiedCmpdAbbrevsRef = LMAPSStr::SetupCmpdAbbrevs(\%Options);

}

# Setup script usage  and retrieve command line arguments specified using various options...
sub SetupScriptUsage {

  # Retrieve all the options...
  %Options = ();
  $Options{mode} = 'Abbrev';

  if (!GetOptions(\%Options, "help|h", "mode|m=s", "overwrite|o", "root|r=s", "workingdir|w=s")) {
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
}


__END__

=head1 NAME

LMStrGen.pl - Generate arbitrary LIPID MAPS structures

=head1 SYNOPSIS

LMStrGen.pl  LMAbbrev|LMAbbrevFileName ...

LMStrGen.pl [B<-h, --help>] [B<-o, --overwrite>] [B<-r, --root> rootname]
[B<-w, --workingdir> dirname] <arguments>...

=head1 DESCRIPTION

Generate arbitrary LIPID MAPS (LM)  structures using compound abbreviations specified on
a command line or in a CSV/TSV Text file. All the command line arguments represent either
compound abbreviations or file name containing abbreviations. Use mode option to control
the type of command line arguments.

A SD file, containing structures for all LM abbreviations along with ontological information, is
generated as an output.

=head1 SUPPORTED ABBREVIATIONS

Current support for LM structure generation include these headgroups and acyl chains:

   o DIMAXX - sn1, sn2 and sn3
   o DIMA20 - sn1 and sn2
   o DIMA22 - sn1 and sn2
   o DIMA20Me - sn1 and sn2
   o DIMA22Me - sn1 and sn2
   o DIMB20 - sn1 and sn2
   o DIMB22 - sn1 and sn2
   o DIB20Me - sn1 and sn2
   o DIB22Me - sn1 and sn2
   o PIM1 - sn1 and sn2
   o PIM2 - sn1 and sn2
   o PIM3 - sn1 and sn2
   o PIM4 - sn1 and sn2
   o PIM5 - sn1 and sn2
   o PIM6 - sn1 and sn2
   o DAT - sn1 and sn2
   o CoA - sn1

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

Wild card character, *, is not supported in compound abbreviations.

Examples:

    Specific structures: DIMA20(12:0/13:0) DIMB20(17:1(9Z)/18:0)
                         CoA(16:0) PIM1(21:0/22:0)

=item B<-o, --overwrite>

Overwrite existing files

=item B<-r, --root> I<rootname>

New file name is generated using the root: <Root>.sdf. Default for new file names: LMAbbrev.sdf,
<AbbrevFilenName>.sdf, or <FirstAbbrevFileName>1To<Count>.sdf.

=item B<-w, --workingdir> I<dirname>

Location of working directory. Default: current directory

=back

=head1 EXAMPLES

On some systems, command line scripts may need to be invoked using
I<perl -s LMStrGen.pl>; however, all the examples assume direct invocation
of command line script works.

To generate a LMStructures.sdf file containing a structure specified
by a command line LM abbreviation, type:

    % LMStrGen.pl -r LMStructures -o "DIMA20(12:0/13:0)"

To generate a LMStructures.sdf file containing structures specified
by a command line LM abbreviations, type:

    % LMStrGen.pl -r LMStructures -o "DIMA22(16:0/18:0)" "DIMA22Me(18:1(11E)/16:0)"

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
