package LMAPSStr;
#
# File: LMAPSStr.pm
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
use Exporter;
use Text::ParseWords;
use FileUtil;
use TextUtil;

our ($VERSION, @ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);

$VERSION = '2.00';
@ISA = qw(Exporter);
@EXPORT = qw();
@EXPORT_OK = qw(GenerateCmpdAtomLine GenerateCmpdBondLine GenerateCmpdCountsLine GenerateCmpdMiscInfoLine ParseCmpdAtomLine ParseCmpdBondLine ParseCmpdCountsLine RoundToNextInteger SetupCmpdAbbrevs SetupSDFileName StandardizeStereochemistrySpec StandardizeRingStereochemistrySpec);
%EXPORT_TAGS = (
		all  => [@EXPORT, @EXPORT_OK],
		);

# Set up SD file count line...
sub GenerateCmpdCountsLine {
  my($AtomCount, $BondCount, $ChiralFlag, $PropertyCount, $Version, $Line);

  if (@_ == 5) {
    ($AtomCount, $BondCount, $ChiralFlag, $PropertyCount, $Version) = @_;
  }
  elsif (@_ == 3) {
    ($AtomCount, $BondCount, $ChiralFlag) = @_;
    $PropertyCount = 999;
    $Version = "V2000";
  }
  else {
    ($AtomCount, $BondCount) = @_;
    $ChiralFlag = 0;
    $PropertyCount = 999;
    $Version = "V2000";
  }
  $Line = sprintf "%3i%3i%3i%3i%3i%3i%3i%3i%3i%3i%3i%6s", $AtomCount, $BondCount, 0, 0, $ChiralFlag, 0, 0, 0, 0, 0, $PropertyCount, $Version;

  return ($Line);
}

# Parse SD file count line taken from www.mayachemtools.org...
sub ParseCmpdCountsLine {
  my($Line) = @_;
  my($AtomCount, $BondCount, $ChiralFlag, $PropertyCount, $Version);

  ($AtomCount, $BondCount, $ChiralFlag, $PropertyCount, $Version) = unpack("A3A3x3x3A3x3x3x3x3x3A3A6", $Line);
  return ($AtomCount, $BondCount, $ChiralFlag, $PropertyCount, $Version);
}
# Setup SD file atom line taken from www.mayachemtools.org...
sub GenerateCmpdAtomLine {
  my($AtomX, $AtomY, $AtomZ, $AtomSymbol) = @_;
  my($Line);

  $Line = sprintf "%10.4f%10.4f%10.4f %-3s%2i%3i%3i%3i%3i%3i%3i%3i%3i%3i%3i%3i", $AtomX, $AtomY, $AtomZ, $AtomSymbol, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;

  return $Line;
}

# Set up SD file comments line...
sub GenerateCmpdCommentsLine {
  my($CommentsLine);

  $CommentsLine = 'Structure generated using tools available at www.lipidmaps.org';

  return $CommentsLine;
}


# Set up SD file misc line taken from www.mayachemtools.org...
sub GenerateCmpdMiscInfoLine {
  my($Sec, $Min, $Hour, $MDay, $Mon, $Year, $WDay, $YDay, $IsDst) = localtime;

  $Mon += 1;
  $Year += 1900;
  $Year = substr($Year, -2, 2);
  my($TimeStamp) = sprintf "%2i%2i%2i%2i%2i", $Mon, $MDay, $Year, $Hour, $Min;
  $TimeStamp =~ s/ /0/g;
  my($MiscLine) = "  LipdMAPS" . "$TimeStamp" . "2D";

  return $MiscLine;
}

# Parse SD file atom line taken from www.mayachemtools.org...
sub ParseCmpdAtomLine {
  my($Line) = @_;
  my ($AtomX, $AtomY, $AtomZ, $AtomSymbol);

  ($AtomX, $AtomY, $AtomZ, $AtomSymbol) = unpack("A10A10A10xA3", $Line);
  return ($AtomX, $AtomY, $AtomZ, $AtomSymbol);
}

# Setup SD file bond line taken from www.mayachemtools.org...
sub GenerateCmpdBondLine {
  my($FirstAtomNum, $SecondAtomNum, $BondType, $BondStereo, $Line);

  if (@_ == 4) {
    ($FirstAtomNum, $SecondAtomNum, $BondType, $BondStereo) = @_;
  }
  else {
    ($FirstAtomNum, $SecondAtomNum, $BondType) = @_;
    $BondStereo = 0;
  }

  $Line = sprintf "%3i%3i%3i%3i%3i%3i%3i", $FirstAtomNum, $SecondAtomNum, $BondType, $BondStereo, 0, 0, 0;
  return $Line;
}

# Parse SD file bond line taken from www.mayachemtools.org...
sub ParseCmpdBondLine {
  my($Line) = @_;
  my($FirstAtomNum, $SecondAtomNum, $BondType, $BondStereo);

  ($FirstAtomNum, $SecondAtomNum, $BondType, $BondStereo) = map {s/ //g; $_} unpack("A3A3A3A3", $Line);
  return ($FirstAtomNum, $SecondAtomNum, $BondType, $BondStereo);
}

# Round it off...
sub RoundToNextInteger {
  my($Number) = @_;
  return int($Number + .5 * ($Number <=> 0));
}

# Get specified compound abbreviations...
sub SetupCmpdAbbrevs {
  my($OptionsRef) = @_;
  my(@CmpdAbbrevs);

  # Collect compounds specifications...
  @CmpdAbbrevs = ();
  my($Param);
  PARAM: for $Param (@ARGV) {
    if ($OptionsRef->{mode} =~ /^Abbrev$/i) {
      push @CmpdAbbrevs, $Param;
    }
    elsif ($OptionsRef->{mode} =~ /^AbbrevFileName$/i) {
      # Process input files...
      my($Line, $AbbrevFile, $FileDir, $FileName, $FileExt, $InDelim);

      $AbbrevFile = $Param;
      ($FileDir, $FileName, $FileExt) = ParseFileName($AbbrevFile);
      $InDelim = "";
      if (lc($FileExt) eq "csv") {
	$InDelim = "\,";
      }
      elsif (lc($FileExt) eq "tsv") {
	$InDelim = "\t";
      }
      elsif (lc($FileExt) eq "txt") {
	$InDelim = "";
      }
      else {
	warn "Warning: Ignoring file $AbbrevFile: Unknown file extension $FileExt : Allowed extentions: csv or tsv \n";
	next PARAM;
      }
      if (! -e $AbbrevFile) {
	warn "Warning: Ignoring file $AbbrevFile: It doesn't exist\n";
	next PARAM;
      }
      if (!open ABBREVFILE, "$AbbrevFile" ) {
	warn "Warning: Ignoring file $AbbrevFile: Couldn't open it: $! \n";
	next PARAM;
      }
      LINE: while ($Line = GetTextLine(\*ABBREVFILE)) {
           # Ignore comments line...
	if ($Line =~ /^#/ || $Line =~ /^-/) {
	  next LINE;
	}
	if ($InDelim) {
	  push @CmpdAbbrevs, quotewords($InDelim, 0, $Line);
	}
	else {
	  $Line =~ s/ //g;
	  push @CmpdAbbrevs, $Line;
	}
      }
      close ABBREVFILE;
    }
  }
  return \@CmpdAbbrevs;
}

# Setup a new SD file name...
sub SetupSDFileName {
  my($SDFile, $LipidCategory, $OptionsRef);

  ($LipidCategory, $OptionsRef) = @_;

  $SDFile = '';
  if ($OptionsRef->{root}) {
    my ($RootFileDir, $RootFileName, $RootFileExt) = ParseFileName($OptionsRef->{root});
    if ($RootFileName && $RootFileExt) {
      $SDFile = $RootFileName;
    }
    else {
      $SDFile = $OptionsRef->{root};
    }
  }
  else {
    $SDFile = $LipidCategory . "Abbrev";
    if ($OptionsRef->{mode} =~ /^AbbrevFileName$/i) {
      my($FileDir, $FileName, $FileExt) = ParseFileName($ARGV[0]);
      my($FileCount) = @ARGV;
      $SDFile = (@ARGV > 1) ? ("$FileName" . "1To" . "$FileCount") : "$FileName";
    }
  }
  $SDFile .= ".sdf";
  return $SDFile;
}

#
# Map alpha/a to R and b/beta to S...
#
sub StandardizeStereochemistrySpec {
  my($InSpec) = @_;
  my($OutSpec);

  if ($InSpec =~ /^(alpha|R|a)$/i) {
    $OutSpec = 'R';
  }
  elsif ($InSpec =~ /^(beta|S|b)$/i) {
    $OutSpec = 'S';
  }
  else {
    $OutSpec = $InSpec;
  }
  return $OutSpec;
}

#
# Map alpha/a to  alpha b/beta to beta...
#
sub StandardizeRingStereochemistrySpec {
  my($InSpec) = @_;
  my($OutSpec);

  if ($InSpec =~ /^(alpha|a)$/i) {
    $OutSpec = 'alpha';
  }
  elsif ($InSpec =~ /^(beta|b)$/i) {
    $OutSpec = 'beta';
  }
  else {
    $OutSpec = $InSpec;
  }
  return $OutSpec;
}


1;

__END__

=head1 NAME

LMAPSStr - Glycerolipids (GL) structure generation methods

=head1 SYNOPSIS

use LMAPSStr;

use LMAPSStr qw(:all);

=head1 DESCRIPTION

LMAPSStr module provides these methods:

    GenerateCmpdAtomLine - Generate SD file atom data line
    GenerateCmpdBondLine - Generate SD file bond data line
    GenerateCmpdCountsLine - Generate SD file count data line
    GenerateCmpdMiscInfoLine - Generate SD file misc data line
    ParseCmpdAtomLine - Parse SD file atom data line
    ParseCmpdBondLine - Parse SD file bond data line
    ParseCmpdCountsLine - Parse SD file count data line
    RoundToNextInteger - Round up to next integer
    SetupCmpdAbbrevs - Setup lipid abbreviations
    SetupSDFileName - Setup SD file name
    StandardizeStereochemistrySpec - Standardize stereochemistry
    StandardizeStereochemistrySpec - Standardize ring stereochemistry

=head1 METHODS

=over 4

=item B<GenerateCmpdAtomLine>

    $Line = GenerateCmpdAtomLine($AtomX, $AtomY, $AtomZ, $AtomSymbol);

Return a formatted atom data line for SD file.

=item B<GenerateCmpdBondLine>

    $Line = GenerateCmpdBondLine($FirstAtomNum, $SecondAtomNum,
            $BondType, [$BondStereo]);

Return a formatted bond data line for SD file.

=item B<GenerateCmpdCountsLine>

    $Line = GenerateCmpdCountsLine($AtomCount, $BondCount,
            [$ChiralFlag, $PropertyCount, $Version]);

Return a formatted count data line for SD file.

=item B<GenerateCmpdMiscInfoLine>

    $Line = GenerateCmpdMiscInfoLine();

Return a formatted miscellaneous data line for SD file. In addition to a time stamp, LipdMAPS
name is used as the program name.

=item B<ParseCmpdAtomLine>

    ($AtomX, $AtomY, $AtomZ, $AtomSymbol) = ParseCmpdAtomLine($Line);

Parse SD file atom data line and return a list with these values: atom coordinates and
element symbol.

=item B<ParseCmpdBondLine>

    ($FirstAtomNum, $SecondAtomNum, $BondType, $BondStereo) =
        ParseCmpdBondLine($Line);

Parse SD file atom bond data line and return a list containing these values: bond atom numbers and
bond type.

=item B<ParseCmpdCountsLine>

    ($AtomCount, $BondCount, $ChiralFlag, $PropertyCount, $Version) =
         ParseCmpdCountsLine($Line);

Parse SD file count data line and return a list containing these values:  atom/bond count and other
miscellaneous count information.

=item B<RoundToNextInteger>

    $IntegerValue = RoundToNextInteger($Number);

Return an integer by rounding the number off to next integer.

=item B<SetupCmpdAbbrevs>

    $AbbrebArrayRef = SetupCmpdAbbrevs($CmdLineOptionsRef);

Return a reference to an array containing specified compound abbreviations by parsing
command line arguments or processing files containing specified abbreviations.

=item B<SetupSDFileName>

    $SDFilename = SetupSDFileName($LipidCategory, $CmdLineOptionsRef);

Return a SD file name by processing a specified I<-r, --root> option or using default
values.

=item B<StandardizeStereochemistrySpec>

    $StandardizeSpec = StandardizeStereochemistrySpec($StereochemistrySpec);

Return a standardize stereochemistry specification containg R/S instead of a/b or
alpha/beta.

=item B<StandardizeRingStereochemistrySpec>

    $StandardizeSpec = StandardizeRingStereochemistrySpec($StereochemistrySpec);

Return a standardize stereochemistry specification containg alpha/beta instead of
a/b or alpha/beta.

=back

=head1 AUTHOR

Manish Sud

=head1 CONTRIBUTOR

Eoin Fahy

=head1 SEE ALSO

ChainAbbrev.pm, ChainStr.pm

=head1 COPYRIGHT

Copyright (C) 2006-2012. The Regents of the University of California. All Rights Reserved.

=head1 LICENSE

Modified BSD License

=cut
