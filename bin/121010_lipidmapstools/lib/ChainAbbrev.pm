package ChainAbbrev;
#
# File: ChainAbbrev.pm
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

use 5.006;
use strict;
use Exporter;
use Text::ParseWords;
use TextUtil;

our ($VERSION, @ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);

$VERSION = '2.00';
@ISA = qw(Exporter);
@EXPORT = qw();
@EXPORT_OK = qw(ChainAbbrevNameExists ExpandChainAbbrev GetChainAbbrevToNameMap GetChainLenAbbrevSupportedMap GetChainLenAbbrevDbleBondGeometryDataMap GetChainLengthAndMultipleBondCount GetSubstituentsAbbrevToNameMap GetChainLenToNamePrefixMap GetCountToNamePrefixMap GetSubstituentBondOrder GetSupportedChainLenList IsAlkylChainAbbrev IsAlkenylChainAbbrev IsAlkenylChainAbbrevWithImplicitDoubleBond IsChainAbbrevOkay IsChainLengthOkay IsDoubleBondsAbbrevOkay IsRingsAbbrevOkay IsSubstituentsAbbrevOkay IsSubstituentInChainAbbrev IsWildCardInChainAbbrev ParseChainAbbrev ParseSubstituentAbbrev ParseRingAbbrev SetupChainName SetupChainNameWithSubstituents SetupChainSubstituentsName SetupChainsSystematicName SetupChainsAbbrev);
%EXPORT_TAGS = (
		all  => [@EXPORT, @EXPORT_OK],
		);

# Initialize data...
my(%ChainAbbrevToNameMap, %ChainLenAbbrevSupportedMap, %ChainLenAbbrevToDbleBondGeomtryMap, %ChainLenToNamePrefixMap, %CountToNamePrefixMap, %SubstituentAbbrevToNameMap);
_InitializeData();

# Expand chain abbreviation by enumerating all possibilites for wild cards
# and return a reference to expanded list...
sub ExpandChainAbbrev {
  my($Abbrev, @ExpandedAbbrevs);

  ($Abbrev) = @_;
  @ExpandedAbbrevs = ();
  if (!IsWildCardInChainAbbrev($Abbrev)) {
    push @ExpandedAbbrevs, $Abbrev;
    return \@ExpandedAbbrevs;
  }
  # Expand wild cards...
  my($NewAbbrev, $ChainLengthAbbrev, $DoubleBondCountAbbrev, $DoubleBondGeometryAbbrev, $SubstituentsAbbrev, $RingsAbbrev, $ChainLength, $DoubleBondCount, $DoubleBondGeometry, @ChainList, @DoubleBondCountList, @DoubleBondGeometryList);

  ($ChainLengthAbbrev, $DoubleBondCountAbbrev, $DoubleBondGeometryAbbrev, $SubstituentsAbbrev, $RingsAbbrev) = ParseChainAbbrev($Abbrev);

  if (IsEmpty($DoubleBondCountAbbrev)) {
    $DoubleBondCountAbbrev = '0';
    if ($ChainLengthAbbrev =~ /\*/) {
      # Handle */*/* ....
      $DoubleBondCountAbbrev = "*";
    }
  }
  if (IsEmpty($DoubleBondGeometryAbbrev)) {
    $DoubleBondGeometryAbbrev = "*";
  }

  if ($ChainLengthAbbrev =~ /\*/) {
    push @ChainList, GetSupportedChainLenList();
  }

  # Chain length  loop...
  CHAINLEN: for $ChainLength (@ChainList) {
    if (!IsChainLengthOkay($ChainLength, $ChainLengthAbbrev)) {
      next CHAINLEN;
    }
    @DoubleBondCountList = ($DoubleBondCountAbbrev =~ /\*/) ? (sort { $a cmp $b } keys %{$ChainLenAbbrevToDbleBondGeomtryMap{$ChainLength}}) : $DoubleBondCountAbbrev;

    # Double bond count loop...
    DOUBLEBONDCOUNT: for $DoubleBondCount (@DoubleBondCountList) {
      if ($DoubleBondCount && !(exists $ChainLenAbbrevToDbleBondGeomtryMap{$ChainLength}{$DoubleBondCount})) {
	next DOUBLEBONDCOUNT;
      }
      @DoubleBondGeometryList = ($DoubleBondGeometryAbbrev =~ /\*/) ? (sort keys %{$ChainLenAbbrevToDbleBondGeomtryMap{$ChainLength}{$DoubleBondCount}}) :  $DoubleBondGeometryAbbrev;

      # Double bond geomerty loop...
      DOUBLEBONDGEOMERTY: for $DoubleBondGeometry (@DoubleBondGeometryList) {
	if ($DoubleBondGeometry && !(exists $ChainLenAbbrevToDbleBondGeomtryMap{$ChainLength}{$DoubleBondCount}{$DoubleBondGeometry})) {
	  next DOUBLEBONDGEOMERTY;
	}

	if ($ChainLengthAbbrev =~ /\*/ && ($DoubleBondCount =~ /[e]/ || $DoubleBondCount =~ /[p]/)) {
	  # Fix alkyl and alkenyl abbreviations: Use P-16:0, O-16:0 formats instead of 16:0p, P-16:O and 16:0e
	  # for alkyl and alkenyls...
	  $NewAbbrev = FixAlkylAndAlkenylAbbrevFormats($ChainLength, $DoubleBondCount, $DoubleBondGeometry);
	}
	else {
	  $NewAbbrev = ($DoubleBondCount && $DoubleBondGeometry) ? ("$ChainLength:$DoubleBondCount($DoubleBondGeometry)") : ("$ChainLength:$DoubleBondCount");
	}
	if (!IsEmpty($SubstituentsAbbrev)) {
	  $NewAbbrev .= "(${SubstituentsAbbrev})";
	}
	if (!IsEmpty($RingsAbbrev)) {
	  $NewAbbrev .= "(${RingsAbbrev})";
	}
	push @ExpandedAbbrevs, $NewAbbrev;
      }

    }
  }

  return \@ExpandedAbbrevs;
}


# During wild card expansion, this function is used to fix alkyl and alkenyl abbreviations: Use
#  P-16:0 instead 16:0p and P-16:0 to handle alkenyl chains, and O-16:0 instead of 16:0e for alkyl chains...
sub FixAlkylAndAlkenylAbbrevFormats {
  my($ChainLength, $DoubleBondCount, $DoubleBondGeometry) = @_;
  my($Abbrev) = '';

  if ($DoubleBondCount =~ /[p]/) {
    # Use P-16:0 as default format...
    $ChainLength = "P-${ChainLength}";
    $DoubleBondCount =~ s/p//g;
    $DoubleBondCount = 0;
    $DoubleBondGeometry = '';
  }
  elsif ($DoubleBondCount =~ /[e]/) {
    # Use 0-16:0, 0-16:1(1Z) format...
    $ChainLength = "O-${ChainLength}";
    $DoubleBondCount =~ s/e//g;
  }
  $Abbrev = ($DoubleBondCount && $DoubleBondGeometry) ? ("$ChainLength:$DoubleBondCount($DoubleBondGeometry)") : ("$ChainLength:$DoubleBondCount");

  return $Abbrev;
}

# Is this a known chain abbreviation?
#
sub ChainAbbrevNameExists {
  my($ChainAbbrev) = @_;
  my($Abbrev, $ChainLen, $BondCount, $BondGeometry, $Substituents, $Rings);

  ($ChainLen, $BondCount, $BondGeometry, $Substituents, $Rings) = ParseChainAbbrev($ChainAbbrev);
  $Abbrev = $BondGeometry ? "${ChainLen}:${BondCount}(${BondGeometry})" : "${ChainLen}:${BondCount}";
  if (exists $ChainAbbrevToNameMap{$Abbrev}) {
    return 1;
  }
  return 0;
}

# Return chain length and number of double/triple bonds...
sub GetChainLengthAndMultipleBondCount {
  my($ChainAbbrev) = @_;
  my($Abbrev, $ChainLength, $BondCount, $BondGeometry, $DoubleBondCount, $TripleBondCount);

  ($ChainLength, $DoubleBondCount, $TripleBondCount) = (0) x 3;
  ($ChainLength, $BondCount, $BondGeometry) = ParseChainAbbrev($ChainAbbrev);

  if ($BondGeometry =~ /Y/i) {
    my($Geometry, @BondGeometryList);
    @BondGeometryList = ();
    @BondGeometryList = split /,/, $BondGeometry;
    for $Geometry (@BondGeometryList) {
      if ($Geometry =~ /Y/i) {
	$TripleBondCount++;
      }
      else {
	$DoubleBondCount++;
      }
    }
  }
  else {
    $DoubleBondCount = $BondCount;
  }

  return ($ChainLength, $DoubleBondCount, $TripleBondCount);
}

# Return a refernce to %ChainAbbrevToNameMap hash...
sub GetChainAbbrevToNameMap {
  return \%ChainAbbrevToNameMap;
}

# Return a refernce to %ChainLenAbbrevSupportedMap hash..
sub GetChainLenAbbrevSupportedMap {
  return \%ChainLenAbbrevSupportedMap, ;
}

# Return a refernce to %ChainLenAbbrevToDbleBondGeomtryMap hash..
sub GetChainLenAbbrevDbleBondGeometryDataMap {
  return \%ChainLenAbbrevToDbleBondGeomtryMap, ;
}

# Return a refernce to %SubstituentAbbrevToNameMap hash...
sub GetSubstituentsAbbrevToNameMap {
  return \%SubstituentAbbrevToNameMap;
}

# Return a refernce to %ChainLenToNamePrefixMap hash...
sub GetChainLenToNamePrefixMap {
  return \%ChainLenToNamePrefixMap;
}

# Return a refernce to %CountToNamePrefixMap hash...
sub GetCountToNamePrefixMap {
  return \%CountToNamePrefixMap;
}


# Get substituent bond order...
sub GetSubstituentBondOrder {
  my($SubstituentAbbrev) = @_;
  my($BondOrder) = 1;

  if ($SubstituentAbbrev =~ /^Ke$/i) {
    $BondOrder = 2;
  }
  elsif ($SubstituentAbbrev =~ /^My$/i) {
    $BondOrder = 2;
  }

  return $BondOrder;
}

# Does chain abbrev specifies an alkyl chain? e.g. 18:0e/16:0/0:0
sub IsAlkylChainAbbrev {
  my($ChainAbbrev) = @_;
  my($Status, $ChainLength, $DoubleBondCount, $DoubleBondGeometry, $Substituents, $Rings);

  $Status = 0;
  ($ChainLength, $DoubleBondCount, $DoubleBondGeometry, $Substituents, $Rings) = ParseChainAbbrev($ChainAbbrev);
  if ($ChainLength =~ /^O-/) {
    # Format: O-16:0
    $Status = ($DoubleBondCount == 0 ) ? 1 : 0;
  }
  else {
    # Format: 16:0e
    #
    # DoubleBondCount would end have e in its specification and DoubleBondGeometry
    # would be empty...
    #
    $Status = ($DoubleBondCount =~ /[e]/ ) ? 1 : 0;
  }

  return $Status;
}

# Using %ChainLenAbbrevSupportedMap, setup a chain len list without any duplicates
# for alkyl/alkenyl chain len format...
#
sub GetSupportedChainLenList {
  my(@SortedChainLenList, @ChainLenList);

  @SortedChainLenList = ();
  @ChainLenList = ();

  # Collect chain lengths without containing "O-" and "P-"in its values; in other words,
  # ignore alkyl/alkeynl formats...
  my($ChainLen);
  CHAINLEN: for $ChainLen (keys %ChainLenAbbrevSupportedMap) {
    if ($ChainLen =~ /^(O-|P-)/i ) {
      next CHAINLEN;
    }
    push @ChainLenList, $ChainLen;
  }

  # Sort 'em out...
  for $ChainLen (sort {$a <=> $b} @ChainLenList) {
    push @SortedChainLenList, $ChainLen;
  }
  return @SortedChainLenList;
}

# Does chain abbrev specifies an alkenyl chain? e.g. 18:0p/16:0/0:0
sub IsAlkenylChainAbbrev {
  my($ChainAbbrev) = @_;
  my($Status, $ChainLength, $DoubleBondCount, $DoubleBondGeometry, $Substituents, $Rings);

  $Status = 0;
  ($ChainLength, $DoubleBondCount, $DoubleBondGeometry, $Substituents, $Rings) = ParseChainAbbrev($ChainAbbrev);

  if ($ChainLength =~ /^O-/) {
    # Format: O-16:1(1Z), O-16:1(9Z), O-16:2(1Z,9Z)
    $Status = ($DoubleBondCount >= 1 ) ? 1 : 0;
  }
  elsif ($ChainLength =~ /^P-/) {
    # Format: P-16:0. An implict double bond at 1Z...
    $Status = 1;
  }
  else {
    # Format: 16:0p
    # $DoubleBondCount would end up with p in its specification...
    #
    $Status = ($DoubleBondCount =~ /[p]/ ) ? 1 : 0;
  }

  return $Status;
}

# Does chain abbrev specifies an alkenyl chain with implict double bond specification
# at 1Z?
#
# P-18:0, 18:0p implies P-18:0(1Z)
# P-18:1(9Z) implies P-18:2(1Z,9Z)
#
sub IsAlkenylChainAbbrevWithImplicitDoubleBond {
  my($ChainAbbrev) = @_;
  my($Status, $ChainLength, $DoubleBondCount, $DoubleBondGeometry, $Substituents, $Rings);

  $Status = 0;

  ($ChainLength, $DoubleBondCount, $DoubleBondGeometry, $Substituents, $Rings) = ParseChainAbbrev($ChainAbbrev);
  if ($ChainLength =~ /^P-/) {
    # Format: P-16:0. An implict double bond at 1Z is implied.
    $Status = 1;
  }
  else {
    # Format: 16:0p
    # $DoubleBondCount would end uo with p in its specification...
    #
    $Status = ($DoubleBondCount =~ /[p]/ ) ? 1 : 0;
  }

  return $Status;
}

# Check out the chain abrbrev...
#
# . Wild cards are not allowed for substituents and rings
#
sub IsChainAbbrevOkay {
  my($ChainAbbrev, $AllowSubstituents, $AllowRings, $AllowArbitraryChainLenSpec)  = @_;
  my($ChainLength, $DoubleBondCount, $DoubleBondGeometry, $Substituents, $Rings);

  $AllowSubstituents = defined $AllowSubstituents ? $AllowSubstituents : 1;
  $AllowRings = defined $AllowRings ? $AllowRings : 1;
  $AllowArbitraryChainLenSpec = defined $AllowArbitraryChainLenSpec ? $AllowArbitraryChainLenSpec : 0;

  ($ChainLength, $DoubleBondCount, $DoubleBondGeometry, $Substituents, $Rings) = ParseChainAbbrev($ChainAbbrev);

  if (!(length $ChainLength)) {
    warn "Warning: Ignoring chain abbreviation $ChainAbbrev : No chain length specified\n";
    return 0;
  }

  if ($ChainLength =~ /\*/ ) {
    # Assume it's okay and wouble be rechecked during expansion with a valid chain length...
    return 1;
  }

  if ($ChainLength && !$AllowArbitraryChainLenSpec) {
    if (!(exists $ChainLenAbbrevSupportedMap{$ChainLength})) {
      warn "Warning: Ignoring chain abbreviation $ChainAbbrev : Unknown  chain length $ChainLength \n";
      return 0;
    }
  }
  my($NoWildCharInCount) = ( $DoubleBondCount =~ /\*/  ) ? 0 : 1;
  my($NoWildCharInGeometry) = (  $DoubleBondGeometry =~ /\*/ ) ? 0 : 1;

  if ($AllowArbitraryChainLenSpec) {
    # Just check to make sure number of specified double bonds matches double bond geometry specifications...
    if ($NoWildCharInCount && $NoWildCharInGeometry) {
      if ($DoubleBondCount == 0 ) {
	if (length $DoubleBondGeometry) {
	  warn "Warning: Ignoring chain abbreviation $ChainAbbrev : Double bond count must be non zero to specify geometry\n";
	  return 0;
	}
      }
      else {
	my(@DoubleBondGeometry) = split /,/, $DoubleBondGeometry;
	if ($DoubleBondCount != scalar @DoubleBondGeometry) {
	  warn "Warning: Ignoring chain abbreviation $ChainAbbrev : Double bond count doesn't match specify geometry\n";
	  return 0;
	}
      }
    }
  }
  else {
    if ($NoWildCharInCount && $NoWildCharInGeometry) {
      # No wild card in double bond count and geometry...
      if ($ChainAbbrev ne "0:0" && !(ChainAbbrevNameExists($ChainAbbrev))) {
	warn "Warning: Ignoring chain abbreviation $ChainAbbrev : Unknown chain, double bound count, or geometry\n";
	return 0;
      }
    }
    elsif ($NoWildCharInCount) {
      # No wild card in double bond count with possibility of wild cards in geometry...
      if ($DoubleBondCount == 0 ) {
	if (length $DoubleBondGeometry) {
	  warn "Warning: Ignoring chain abbreviation $ChainAbbrev : Double bond count must be non zero to specify geometry\n";
	  return 0;
	}
      }
      else {
	if (!(exists $ChainLenAbbrevToDbleBondGeomtryMap{$ChainLength}{$DoubleBondCount})) {
	  warn "Warning: Ignoring chain abbreviation $ChainAbbrev : Unknown double bond count\n";
	  return 0;
	}
	if ($NoWildCharInGeometry && length $DoubleBondGeometry) {
	  if (!(exists $ChainLenAbbrevToDbleBondGeomtryMap{$ChainLength}{$DoubleBondCount}{$DoubleBondGeometry})) {
	    warn "Warning: Ignoring chain abbreviation $ChainAbbrev : Unknown double bond geometry\n";
	    return 0;
	  }
	}
      }
    }
    else {
      # No wild card in double bond geometry...
      if ($NoWildCharInGeometry && length $DoubleBondGeometry) {
	my(@DoubleBondCountList, $GeometryFound, $BondCount);
	@DoubleBondCountList =  sort keys %{$ChainLenAbbrevToDbleBondGeomtryMap{$ChainLength}};
	$GeometryFound = 0;
	BONDCOUNT: for $BondCount (@DoubleBondCountList) {
	  if (exists $ChainLenAbbrevToDbleBondGeomtryMap{$ChainLength}{$BondCount}{$DoubleBondGeometry}) {
	    $GeometryFound = 1;
	    last BONDCOUNT;
	  }
	}
	if (!$GeometryFound) {
	  warn "Warning: Ignoring chain abbreviation $ChainAbbrev : Unknown double bond geometry\n";
	  return 0;
	}
      }
    }
  }

  if ($Substituents) {
    if (!$AllowSubstituents) {
      warn "Warning: Ignoring chain abbreviation $ChainAbbrev : Substituents are not allowed\n";
      return 0;
    }
    if (!IsSubstituentsAbbrevOkay($ChainAbbrev, $ChainLength, $DoubleBondCount, $DoubleBondGeometry, $Substituents)) {
      return 0;
    }
  }

  if ($Rings) {
    if (!$AllowRings) {
      warn "Warning: Ignoring chain abbreviation $ChainAbbrev : Rings are not allowed\n";
      return 0;
    }
    if (!IsRingsAbbrevOkay($ChainAbbrev, $ChainLength, $Rings)) {
      return 0;
    }
  }

  return 1;
}

# Check to see whether a specific chain length is part of chain abbreviation encompassing
# wild cards.
#
# Notes:
# . Check for  even/odd (+/-) chain lengths in abbreviation with wildcards.
# . Check for  >ChanLen and <ChainLen in abbreviation with wildcards.
# . Check combinations of +>ChainLen and ->ChainLen<ChainLen in abbreviation with wildcards.
#
sub IsChainLengthOkay {
  my($ChainLength, $ChainLengthAbbrev) = @_;

  if ($ChainLengthAbbrev !~ /\*/) {
    return 1;
  }

  # Check even or odd chain length specifications...
  if ($ChainLengthAbbrev =~ /\+/) {
    # Even chain length specified...
    if ($ChainLength % 2) {
	return 0;
      }
  }
  elsif ($ChainLengthAbbrev =~ /\-/) {
    # Even chain length specified...
    if (!($ChainLength % 2)) {
      return 0;
    }
  }

  # Check explicit chain length specifications...
  my($Length) = 0;
  if ($ChainLengthAbbrev =~ /\>/) {
    ($Length) = ($ChainLengthAbbrev =~ /\>([0-9]+)/);
    if ($ChainLength < $Length) {
      return 0;
    }
  }
  if ($ChainLengthAbbrev =~ /\</) {
    ($Length) = ($ChainLengthAbbrev =~ /\<([0-9]+)/);
    if ($ChainLength > $Length) {
      return 0;
    }
  }
  return 1;
}

#
# Make sure it's a valid rings abbreviation...
#
#  . Only 5 and 6 membered rings are supported.
#  . Stereochemistry: R/S, a/b, alpha/beta or none
#  . Ring positions must be less than chain length...
#
#  e.g.
#   8a, 12b
#   8R, 12S
#   8alpha, 12beta
#   8, 12
#
sub IsRingsAbbrevOkay {
  my($ChainAbbrev, $ChainLength, $Rings) = @_;
  if ($Rings =~ /\*/) {
    warn "Warning: Ignoring chain abbreviation $ChainAbbrev : Wild card found in rings abbreviation\n";
    return 0;
  }
  my(@RingWords);
  @RingWords = quotewords(',', 0, $Rings);
  if (@RingWords <= 1 || @RingWords > 2) {
    warn "Warning: Ignoring chain abbreviation $ChainAbbrev : Rings specificaton, $Rings, must contain only two comma delimited values\n";
    return 0;
  }
  my($RingAbbrev, $Pos, $StereoChemistry, @RingPositions, @RingStreroChemistry);
  @RingPositions = (); @RingStreroChemistry = ();
  for $RingAbbrev (@RingWords) {
    ($Pos, $StereoChemistry) = ParseRingAbbrev($RingAbbrev);
    if (IsEmpty($Pos)) {
      warn "Warning: Ignoring chain abbreviation $ChainAbbrev : Ring position specificaton is empty\n";
      return 0;
    }
    if ($Pos > $ChainLength) {
      warn "Warning: Ignoring chain abbreviation $ChainAbbrev : Ring position specificaton, $Pos, must be  less than chain length\n";
      return 0;
    }
    if (!IsEmpty($StereoChemistry)) {
      if ($StereoChemistry !~ /^(a|b|alpha|beta|U)$/) {
	warn "Warning: Ignoring chain abbreviation $ChainAbbrev : Uknown ring position stereochemistry specificaton $StereoChemistry. Supported values: a, b, alpha or beta.\n";
	return 0;
      }
    }
    push @RingPositions, $Pos;
    push @RingStreroChemistry, $StereoChemistry;
  }
  if ($RingPositions[0] >= $RingPositions[1]) {
      warn "Warning: Ignoring chain abbreviation $ChainAbbrev : First ring position specificaton, $RingPositions[0], must be  less than second ring position specification $RingPositions[1]\n";
      return 0;
  }
  my($RingSize) = $RingPositions[1] - $RingPositions[0] + 1;
  if ($RingSize < 5 || $RingSize > 6) {
    warn "Warning: Ignoring chain abbreviation $ChainAbbrev : Unknown ring size, $RingSize\n";
    return 0;
  }
  # Make sure ring specification lies with in the chain length...
  if ($RingPositions[0] >= $ChainLength ||  $RingPositions[1] > $ChainLength) {
      warn "Warning: Ignoring chain abbreviation $ChainAbbrev : First and second ring positions, $RingPositions[0] and $RingPositions[1], must be less than the chain length, $ChainLength...\n";
      return 0;
  }

  return 1;
}

#
# Make sure it's a valid double bond abbreviation for any arbitrary chain length and bonds...
#
#   . Number of double bonds and corresponding specification match
#   . Only double and triple bonds are allowed (Y and Z)
#   . No multiple bond specification at position 1
#   . No double and triple bond specification next to each other
#   . Bond geometry specification not allowed at terminal chain atom
#   . All position specification are less than chain length
#
sub IsDoubleBondsAbbrevOkay {
  my($ChainAbbrev, $ChainLength, $DoubleBondCount, $DoubleBondGeometry) = @_;
  my($MultipleBondGeometry, $Position, $Geometry, @DoubleBondGeometryList, %DoubleBondChainAtomPosToGeometry);

  %DoubleBondChainAtomPosToGeometry = ();
  @DoubleBondGeometryList = ();
  @DoubleBondGeometryList = split /,/, $DoubleBondGeometry;

  if ($DoubleBondCount != @DoubleBondGeometryList) {
    warn "Warning: Ignoring chain abbreviation $ChainAbbrev : Mutiple bond geometry specification, $DoubleBondGeometry, doesn't match the Total number of multiple bonds, $DoubleBondCount, specified.\n";
    return 0;
  }

  # Check geometry specification...
  for $MultipleBondGeometry (@DoubleBondGeometryList) {
    ($Position, $Geometry) = ($MultipleBondGeometry =~ /^([0-9]+)([a-zA-Z]+)$/);
    if (!defined($Position) || !defined($Geometry)) {
      warn "Warning: Ignoring chain abbreviation $ChainAbbrev : Invalid bond geometry specification: $MultipleBondGeometry\n";
      return 0;
    }
    $DoubleBondChainAtomPosToGeometry{$Position} = $Geometry;
    if ($Geometry && $Geometry !~ /^[EZY]$/i) {
      warn "Warning: Ignoring chain abbreviation $ChainAbbrev : Unknown bond geometry, $Geometry, specified at position, $Position, in $DoubleBondGeometry: Allowed values: E, Z or Y.\n";
      return 0;
    }
    if ($Position == 1) {
      warn "Warning: Ignoring chain abbreviation $ChainAbbrev : Bond position and geometry specification, ${Position}${Geometry}, not supported at chain atom position 1. \n";
      return 0;
    }
    if ($Position >= $ChainLength) {
      warn "Warning: Ignoring chain abbreviation $ChainAbbrev : Bond position in , ${Position}${Geometry}, must be less than the chain length, $ChainLength\n";
      return 0;
    }
  }

  # Check positions...
  my($NextGeometry, $NextPosition);
  for $Position (sort keys %DoubleBondChainAtomPosToGeometry) {
    $Geometry = $DoubleBondChainAtomPosToGeometry{$Position};
    $NextPosition = $Position + 1;
    if (exists $DoubleBondChainAtomPosToGeometry{$NextPosition}) {
      $NextGeometry = $DoubleBondChainAtomPosToGeometry{$NextPosition};
      warn "Warning: Ignoring chain abbreviation $ChainAbbrev : Invalid positions in bond specifications,  ${Position}${Geometry},${NextPosition}${NextGeometry}: Position specification not supported at two sequential chain atoms.\n";
      return 0;
    }
  }

  return 1;
}

#
# Make sure it's a valid substituent abbreviation...
#
#   . Stereo chemistry for substituents doesn't contain any numbers
#   . Substituent position must be > 1 and <= chain length (< chain length for epoxy).
#     By default substituent at position 1 is not allowed: This position is usually part of
#     the template and can cause problems in structure generation.
#   . Allow only upto two sustituents 2, 1 and no substitents for sp3, sp2 and sp carbons except
#     for terminal chain atom which could have one extra substituent at each position.
#   . Stereochemistry specification is only allowed for one substituent.
#
#  e.g.
#   9Ke,15OH[S]
#   9Ke,12OH[S],14Ep[R S])
#
sub IsSubstituentsAbbrevOkay {
  my($ChainAbbrev, $ChainLength, $DoubleBondCount, $DoubleBondGeometry, $Substituents) = @_;

  if ($Substituents =~ /\*/) {
    warn "Warning: Ignoring chain abbreviation $ChainAbbrev : Wild card found in substituents abbreviation\n";
    return 0;
  }
  if ($ChainLength =~ /^O-/) {
    $ChainLength =~ s/^O-//g;
  }
  elsif ($ChainLength =~ /^P-/) {
    $ChainLength =~ s/^P-//g;
  }
  # Make sure it's a known substituent and the position specified is less than the chain length...
  my($Substituent, $SubstituentPos, $SubstituentAbbrev, $SubstituentStereoChemistry, @SubstituentsList, %SubstituentsDataMap);
  %SubstituentsDataMap = ();
  %{$SubstituentsDataMap{SubstituentsPos}} = ();
  %{$SubstituentsDataMap{SubstituentsAtPosCount}} = ();
  %{$SubstituentsDataMap{SubstituentsAtPos}} = ();
  %{$SubstituentsDataMap{SubstituentsAtPosStereoSpecCount}} = ();
  %{$SubstituentsDataMap{SubstituentsAtPosStereoSpec}} = ();
  %{$SubstituentsDataMap{SubstituentsAtPosBondOrderCount}} = ();

  (@SubstituentsList) = split /\,/, $Substituents;
 SUBSTITUENT: for $Substituent (@SubstituentsList) {
    ($SubstituentPos, $SubstituentAbbrev, $SubstituentStereoChemistry) = ParseSubstituentAbbrev($Substituent);

    # Validate position value...
    if (IsEmpty($SubstituentPos)) {
      warn "Warning: Ignoring chain abbreviation $ChainAbbrev : Missing substituent position value\n";
      return 0;
    }
    if ($SubstituentPos =~ /[^0-9]+/) {
      warn "Warning: Ignoring chain abbreviation $ChainAbbrev : Substituent position value must be an integer\n";
      return 0;
    }
    if ($SubstituentPos <= 1 && $SubstituentAbbrev !~ /^(OH|NH2|CHO|COOH)$/i) {
      # Allowed for OH, CHO, NH2, to support Fatty Alcohols, aldehydes, primary amides using different templates...
      warn "Warning: Ignoring chain abbreviation $ChainAbbrev : Substituent position value, $SubstituentPos, for substituent, $SubstituentAbbrev, must be greater than 1.\n";
      return 0;
    }
    if ($SubstituentPos > 1 && $SubstituentAbbrev =~ /^(COOH)$/i) {
      # COOH only allowed at position 1...
      warn "Warning: Ignoring chain abbreviation $ChainAbbrev : Substituent position value, $SubstituentPos, for substituent, $SubstituentAbbrev, must be 1.\n";
      return 0;
    }

    if ($SubstituentPos > $ChainLength) {
      warn "Warning: Ignoring chain abbreviation $ChainAbbrev : Substituent position value, $SubstituentPos, must be less than chain length\n";
      return 0;
    }
    if ($SubstituentAbbrev =~ /[CE]p/) {
      if ($SubstituentPos >= $ChainLength) {
	warn "Warning: Ignoring chain abbreviation $ChainAbbrev : Epoxy substituent position value, $SubstituentPos, must be less than chain length minus one\n";
	return 0;
      }
    }
    if (!exists $SubstituentAbbrevToNameMap{$SubstituentAbbrev}) {
      warn "Warning: Ignoring chain abbreviation $ChainAbbrev : Unknown substituent $SubstituentAbbrev\n";
      return 0;
    }

    if (exists $SubstituentsDataMap{SubstituentsPos}{$SubstituentPos}) {
      $SubstituentsDataMap{SubstituentsAtPos}{$SubstituentPos} .= ",${SubstituentAbbrev}";
      $SubstituentsDataMap{SubstituentsAtPosCount}{$SubstituentPos} += 1;
      $SubstituentsDataMap{SubstituentsAtPosStereoSpec}{$SubstituentPos} .= ",${SubstituentStereoChemistry}";
      if (!IsEmpty($SubstituentStereoChemistry)) {
	$SubstituentsDataMap{SubstituentsAtPosStereoSpecCount}{$SubstituentPos} += 1;
      }
      $SubstituentsDataMap{SubstituentsAtPosBondOrderCount}{$SubstituentPos} += GetSubstituentBondOrder($SubstituentAbbrev);
    }
    else {
      $SubstituentsDataMap{SubstituentsPos}{$SubstituentPos} = $SubstituentPos;
      $SubstituentsDataMap{SubstituentsAtPos}{$SubstituentPos} = $SubstituentAbbrev;
      $SubstituentsDataMap{SubstituentsAtPosCount}{$SubstituentPos} = 1;
      $SubstituentsDataMap{SubstituentsAtPosStereoSpec}{$SubstituentPos} = $SubstituentStereoChemistry;
      $SubstituentsDataMap{SubstituentsAtPosStereoSpecCount}{$SubstituentPos} = IsEmpty($SubstituentStereoChemistry) ? 0 : 1;
      $SubstituentsDataMap{SubstituentsAtPosBondOrderCount}{$SubstituentPos} = GetSubstituentBondOrder($SubstituentAbbrev);
    }

    # Validate SubstituentStereoChemistry value...
    # Format: $StereoChemistry and for Ep: $Number$StereoChemistry
    if (IsEmpty($SubstituentStereoChemistry)) {
      next SUBSTITUENT;
    }
    if ($SubstituentStereoChemistry =~ /[^URSab ]/) {
      warn "Warning: Ignoring chain abbreviation $ChainAbbrev : Unknown substituent stereo chemistry $SubstituentStereoChemistry\n";
      return 0;
    }
    if ($SubstituentAbbrev =~ /^[CE]p$/i && $SubstituentStereoChemistry) {
      my(@StereoChemistryWords);
      @StereoChemistryWords = quotewords(' ', 0, $SubstituentStereoChemistry);
      if (@StereoChemistryWords > 2) {
	warn "Warning: Ignoring chain abbreviation $ChainAbbrev :  Stereochemistry specified for more than two epoxy positions\n";
	return 0;
      }
      # Not supported at for now...
      warn "Warning: Ignoring chain abbreviation $ChainAbbrev :  Stereochemistry specification for epoxy is not supported.\n";
      return 0;
    }
  }
  # Make sure substitunts at position 1 are valid...
  if (exists $SubstituentsDataMap{SubstituentsPos}{1}) {
    if ($SubstituentsDataMap{SubstituentsAtPosCount}{1} > 1) {
      warn "Warning: Ignoring chain abbreviation $ChainAbbrev :  More than one substituent specified at position 1.\n";
      return 0;
    }
    if ($SubstituentsDataMap{SubstituentsAtPosStereoSpecCount}{1} >= 1) {
      warn "Warning: Ignoring chain abbreviation $ChainAbbrev :  Stereochemistry specification at position 1 is nore supported.\n";
      return 0;
    }
  }

  # Do a check for thsese conditions:
  #   . Allow only upto two sustituents 2, 1 and no substitents for sp3, sp2 and sp carbons except
  #     for terminal chain atom which could have one extra substituent at each position.
  #   . Stereochemistry specification is only allowed for one substituent.
  #   . Check for multiple bonds in substituents and chiral center to figure out maximum number of
  #     substituents.

  my($SubstituentsAtPos, $SubstituentsAtPosCount, $SubstituentsAtPosStereoSpec, $SubstituentsAtPosStereoSpecCount, $SubstituentsAtPosBondOrderCount, $ChainAtomPosBondCount, $MultipleBondGeometry, $BeforeSubstituentPos, $Position, $Geometry, @DoubleBondGeometryList, %DoubleBondChainAtomPosToGeometry);

  for $SubstituentPos (sort {$a <=> $b} keys %{$SubstituentsDataMap{SubstituentsPos}}) {
    $SubstituentsAtPos = $SubstituentsDataMap{SubstituentsAtPos}{$SubstituentPos};
    $SubstituentsAtPosCount = $SubstituentsDataMap{SubstituentsAtPosCount}{$SubstituentPos};
    $SubstituentsAtPosStereoSpec = $SubstituentsDataMap{SubstituentsAtPosStereoSpec}{$SubstituentPos};
    $SubstituentsAtPosStereoSpecCount = $SubstituentsDataMap{SubstituentsAtPosStereoSpecCount}{$SubstituentPos};
    $SubstituentsAtPosBondOrderCount = $SubstituentsDataMap{SubstituentsAtPosBondOrderCount}{$SubstituentPos};

    if ($SubstituentsAtPosStereoSpecCount > 1) {
      warn "Warning: Ignoring chain abbreviation $ChainAbbrev : Stereochemistry specified for more than one substituent at position $SubstituentPos.\n";
      return 0;
    }

    # Is substituent chain atom position involved in multiple bonds...
    @DoubleBondGeometryList = ();
    %DoubleBondChainAtomPosToGeometry = ();

   @DoubleBondGeometryList = split /,/, $DoubleBondGeometry;
    for $MultipleBondGeometry (@DoubleBondGeometryList) {
      ($Position, $Geometry) = ($MultipleBondGeometry =~ /^([0-9]+)([a-zA-Z]+)$/);
      $DoubleBondChainAtomPosToGeometry{$Position} = $Geometry;
    }
    $BeforeSubstituentPos = $SubstituentPos - 1;
    if ($SubstituentPos == $ChainLength) {
      $ChainAtomPosBondCount = 1;
      if (exists $DoubleBondChainAtomPosToGeometry{$BeforeSubstituentPos}) {
	if ($DoubleBondChainAtomPosToGeometry{$BeforeSubstituentPos} =~ /^Y$/i) {
	  $ChainAtomPosBondCount = 3;
	}
	else {
	  $ChainAtomPosBondCount = 2;
	}
      }
    }
    else {
      $ChainAtomPosBondCount = 2;
      if (exists $DoubleBondChainAtomPosToGeometry{$SubstituentPos}) {
	if ($DoubleBondChainAtomPosToGeometry{$SubstituentPos} =~ /^Y$/i) {
	  $ChainAtomPosBondCount = 4;
	}
	else {
	  $ChainAtomPosBondCount = (exists $DoubleBondChainAtomPosToGeometry{$BeforeSubstituentPos}) ? 4 : 3;
	}
      }
    }
    # Make sure total bond count including substituents and chain atoms is less than 4...
    if (($ChainAtomPosBondCount + $SubstituentsAtPosBondOrderCount) > 4) {
      warn "Warning: Ignoring chain abbreviation $ChainAbbrev : Total number of substituents, $SubstituentsAtPosCount, specified at position $SubstituentPos  is not valid.\n";
      return 0;
    }
  }

  return 1;
}

# Check for substituents in chain abbreviation...
#
sub IsSubstituentInChainAbbrev {
  my($Abbrev) = @_;
  my($ChainLengthAbbrev, $DoubleBondCountAbbrev, $DoubleBondGeometryAbbrev, $SubstituentsAbbrev, $RingsAbbrev, $Status);

  ($ChainLengthAbbrev, $DoubleBondCountAbbrev, $DoubleBondGeometryAbbrev, $SubstituentsAbbrev) = ParseChainAbbrev($Abbrev);
  $Status =  IsEmpty($SubstituentsAbbrev) ? 0 : 1;

  return $Status;
}

# Check for wild card in chain length, double bond count or double bond
# geometry specifications...
#
sub IsWildCardInChainAbbrev {
  my($Abbrev) = @_;

  my($ChainLengthAbbrev, $DoubleBondCountAbbrev, $DoubleBondGeometryAbbrev, $SubstituentsAbbrev, $RingsAbbrev, $Status);

  ($ChainLengthAbbrev, $DoubleBondCountAbbrev, $DoubleBondGeometryAbbrev, $SubstituentsAbbrev, $RingsAbbrev) = ParseChainAbbrev($Abbrev);
  $Status = (($ChainLengthAbbrev =~ /\*/) || ($DoubleBondCountAbbrev =~ /\*/) || ($DoubleBondGeometryAbbrev =~ /\*/)) ? 1 : 0;
  return $Status;
}

# Parse chain abbreviation to retrieve chain length, double bonds count and their geometry...
# Format:
#   ChainLen:BondCount(BondGeometry)(Substituents){Ring}
#
#   . Optional fields: Substituents and Rings; BondGeometry when BondCount = 0
#   . Stereo chemistry for substituents doesn't contain any numbers
#   . Substituent position must be > 0 and <= chain length (< chain length for epoxy)
#
#  e.g.
#   18:1(9Z)(6Ke,15OH[S]){8,9}
#   18:1(9Z)(6Ke,12OH[S],14Ep[R S]){8,9}
#   18:1(9Z)
#   18:0(9Ke,15OH)
#   and so on...
#
sub ParseChainAbbrev {
  my($ChainAbbrev) = @_;
  my($ChainLength, $DoubleBondCount, $DoubleBondGeometry, $Substituents, $Rings, $ChainSpec, $LeftOverChainSpec);

  $ChainLength = ''; $DoubleBondCount = ''; $DoubleBondGeometry = ''; $Substituents = ''; $Rings = '';
  $ChainSpec = $ChainAbbrev; $LeftOverChainSpec = '';
  ($ChainLength, $DoubleBondCount, $LeftOverChainSpec) = $ChainSpec =~ /^(.+?)\:(.+?)(?=\(|\{|$)(.*?)$/;

  $ChainLength = defined $ChainLength ? $ChainLength : '';
  $DoubleBondCount = defined $DoubleBondCount ? $DoubleBondCount : '';

  if (IsEmpty($ChainLength)) {
    #In all likelyhood, it's a wild card abbreviation. Just assign the whole abbrev to chain
    # length and return: the caller would figure out what to do...
    $ChainLength = $ChainAbbrev;
    return ($ChainLength, $DoubleBondCount, $DoubleBondGeometry, $Substituents, $Rings);
  }

  $ChainLength = defined $ChainLength ? $ChainLength : '';
  $DoubleBondCount = defined $DoubleBondCount ? $DoubleBondCount : '';
  $ChainSpec = defined $LeftOverChainSpec ? $LeftOverChainSpec : ''; $LeftOverChainSpec = '';
  if ($ChainSpec && $DoubleBondCount) {
    ($DoubleBondGeometry, $LeftOverChainSpec) = $ChainSpec =~ /^\((.+?)\)(.*?)$/;
    $DoubleBondGeometry = defined $DoubleBondGeometry ? $DoubleBondGeometry : '';
    $ChainSpec = defined $LeftOverChainSpec ? $LeftOverChainSpec : ''; $LeftOverChainSpec = '';
  }
  if ($ChainSpec =~ /^\(/) {
    # Substituents...
    ($Substituents, $LeftOverChainSpec) = $ChainSpec =~ /^\((.+?)\)(.*?)$/;
    $Substituents = defined $Substituents ? $Substituents : '';
    $ChainSpec = defined $LeftOverChainSpec ? $LeftOverChainSpec : ''; $LeftOverChainSpec = '';
  }
  if ($ChainSpec && ($ChainSpec =~ /^\{/)) {
    # Rings...
    ($Rings, $LeftOverChainSpec) = $ChainSpec =~ /^\{(.+?)\}(.*?)$/;
    $Rings = defined $Rings ? $Rings : '';
    $ChainSpec = defined $LeftOverChainSpec ? $LeftOverChainSpec : ''; $LeftOverChainSpec = '';
  }
  return ($ChainLength, $DoubleBondCount, $DoubleBondGeometry, $Substituents, $Rings);
}

# Process substituent abbreviation...
# Format:
#   SubstituentPosSubstituentName[SubstituentStereoChemistry]
#
#   Optional field: SubstituentStereoChemistry
#
#  e.g.
#   9Ke
#   15OH[S]
#   14Ep[14R 15S]
sub ParseSubstituentAbbrev {
  my($Abbrev) = @_;
  my($Pos, $Name, $StereoChemistry, $AbbrevSpec, $LeftOverSpec);

  $Pos = ''; $Name = ''; $StereoChemistry = '';
  $AbbrevSpec = $Abbrev; $LeftOverSpec = '';
  ($Pos, $LeftOverSpec) = $AbbrevSpec =~ /^([0-9]+)(.*?)$/;

  $Pos = defined $Pos ? $Pos : '';
  $AbbrevSpec = defined $LeftOverSpec ? $LeftOverSpec : $Abbrev; $LeftOverSpec = '';

  ($Name, $LeftOverSpec) = $AbbrevSpec =~ /^(.+?)(?=\[|$)(.*?)$/;
  $Name = defined $Name ? $Name : '';
  $AbbrevSpec = defined $LeftOverSpec ? $LeftOverSpec : ''; $LeftOverSpec = '';

  if ($AbbrevSpec && ($AbbrevSpec =~ /^\[/)) {
    ($StereoChemistry, $LeftOverSpec) = $AbbrevSpec =~ /^\[(.+?)\](.*?)$/;
    $StereoChemistry = defined $StereoChemistry ? $StereoChemistry : '';
    $AbbrevSpec = defined $LeftOverSpec ? $LeftOverSpec : ''; $LeftOverSpec = '';
  }
  return ($Pos, $Name, $StereoChemistry);
}

# Process ring abbreviation...
# Format:
#   $Pos$StereoChemistry
#
#   Optional field: StereoChemistry
#
#  e.g.
#   8a,12b
#   8,12
#   8R,12S
sub ParseRingAbbrev {
  my($Abbrev) = @_;
  my($Pos, $StereoChemistry, $AbbrevSpec, $LeftOverSpec);

  $Pos = ''; $StereoChemistry = '';
  $AbbrevSpec = $Abbrev; $LeftOverSpec = '';
  ($Pos, $LeftOverSpec) = $AbbrevSpec =~ /^([0-9]+)(.*?)$/;

  $Pos = defined $Pos ? $Pos : '';
  $StereoChemistry = defined $LeftOverSpec ? $LeftOverSpec : '';

  return ($Pos, $StereoChemistry);
}

# Set up chain name from chain abbreviation...
#
sub SetupChainNameWithSubstituents {
  my($CmpdAbbrevTemplateDataMapRef, $ChainIndex) = @_;
  my($ChainName, $SubstituentsName, $ChainNameWithSubstituents);

  $ChainName = SetupChainName($CmpdAbbrevTemplateDataMapRef, $ChainIndex);
  $SubstituentsName = SetupChainSubstituentsName($CmpdAbbrevTemplateDataMapRef, $ChainIndex);

  $ChainNameWithSubstituents = $ChainName;
  if (!IsEmpty($SubstituentsName)) {
    $ChainNameWithSubstituents =~ s/^\(//g;
    $ChainNameWithSubstituents =~ s/\)$//g;
    $ChainNameWithSubstituents = "($SubstituentsName-$ChainNameWithSubstituents)";
  }

  return $ChainNameWithSubstituents;
}


# Set up chain name from chain abbreviation...
#
sub SetupChainName {
  my($CmpdAbbrevTemplateDataMapRef, $ChainIndex) = @_;
  my($ChainName, $ChainAbbrev);

  $ChainAbbrev = $CmpdAbbrevTemplateDataMapRef->{SnChainAdd}[$ChainIndex] ? $CmpdAbbrevTemplateDataMapRef->{SnAbbrev}[$ChainIndex] : '0:0';

  $ChainName = '';
  if (IsEmpty($ChainAbbrev) || $ChainAbbrev eq "0:0") {
    return $ChainName;
  }

  # Is it a most likely chain abbreviation?
  if (exists $ChainAbbrevToNameMap{$ChainAbbrev}) {
    $ChainName = $ChainAbbrevToNameMap{$ChainAbbrev};
    return $ChainName;
  }

  # Use chain length and double bond geometry specifications to generate chain name...

  my($ChainNamePrefix, $ChainNameMiddle, $ChainNameSuffix, $ChainLength, $DoubleBondCount, $DoubleBondGeometry, $BondPos, $BondGeometry);

  $ChainLength = $CmpdAbbrevTemplateDataMapRef->{SnChainLength}[$ChainIndex];
  $DoubleBondCount = $CmpdAbbrevTemplateDataMapRef->{SnDblBondCount}[$ChainIndex];

  $DoubleBondGeometry = '';
  for $BondPos (sort {$a <=> $b} keys %{$CmpdAbbrevTemplateDataMapRef->{SnDbleBondPosToGeometry}[$ChainIndex]}) {
    $BondGeometry = $CmpdAbbrevTemplateDataMapRef->{SnDbleBondPosToGeometry}[$ChainIndex]{$BondPos};
    if ($DoubleBondGeometry) {
      $DoubleBondGeometry .= ",${BondPos}${BondGeometry}";
    }
    else {
      $DoubleBondGeometry = "${BondPos}${BondGeometry}";
    }
  }

  $ChainNamePrefix = exists $ChainLenToNamePrefixMap{$ChainLength} ? $ChainLenToNamePrefixMap{$ChainLength} : '';
  if (!$ChainNamePrefix) {
    warn "Warning: Problems in setting chain name: Chain name prefix is not known for chain length $ChainLength...\n";
  }

  $ChainNameSuffix = '';
  if (IsAlkylChainAbbrev($ChainAbbrev)) {
    $ChainNameSuffix = 'yl';
  }
  elsif (IsAlkenylChainAbbrev($ChainAbbrev)) {
    $ChainNameSuffix = 'enyl';
  }
  else {
    $ChainNameSuffix = $DoubleBondCount ? 'enoyl' : 'anoyl';
  }

  $ChainNameMiddle = '';
  if ($DoubleBondCount > 1) {
    if (exists $CountToNamePrefixMap{$DoubleBondCount}) {
      $ChainNameMiddle = "a$CountToNamePrefixMap{$DoubleBondCount}";
    }
    else {
      warn "Warning: Problems in setting chain name: Double bond count prefix is not known for $DoubleBondCount...\n";
    }
  }

  $ChainName = "${ChainNamePrefix}${ChainNameMiddle}${ChainNameSuffix}";
  if ($DoubleBondCount && $DoubleBondGeometry) {
    $ChainName = "($DoubleBondGeometry-${ChainName})";
  }

  return $ChainName;
}

# Setup substituents name...
sub SetupChainSubstituentsName {
  my($CmpdAbbrevTemplateDataMapRef, $ChainIndex) = @_;
  my($SubstituentsName, $SubstituentPosNum, $SubstituentsCount, $SubstituentIndex, $SubstituentAbbrev, $SubstituentStereoChemistry, $SubstituentName, $SubtituentCount, $CountPrefix, $SubstituentPosPrefix, %SubstituentsDataMap);

  $SubstituentsName = '';
  if (!$CmpdAbbrevTemplateDataMapRef->{SnSubstituentsCount}[$ChainIndex]) {
    return $SubstituentsName;
  }
  # To allow sorting and counting at the current position, setup a hash using substituent abbreviation...
  %SubstituentsDataMap = ();
  @{$SubstituentsDataMap{Name}} = ();
  %{$SubstituentsDataMap{Abbrev}} = ();
  %{$SubstituentsDataMap{Count}} = ();
  %{$SubstituentsDataMap{NameWithCount}} = ();
  %{$SubstituentsDataMap{PosNumAndStereoPrefix}} = ();
  %{$SubstituentsDataMap{FirstPosNum}} = ();

  for $SubstituentPosNum (sort {$a <=> $b} keys %{$CmpdAbbrevTemplateDataMapRef->{SnSubstituentsPosInfo}[$ChainIndex]}) {
    $SubstituentsCount = @{$CmpdAbbrevTemplateDataMapRef->{SnSubstituentsPosInfo}[$ChainIndex]{$SubstituentPosNum}{Spec}};

    for $SubstituentIndex (0 .. ($SubstituentsCount - 1)) {
      $SubstituentAbbrev = $CmpdAbbrevTemplateDataMapRef->{SnSubstituentsPosInfo}[$ChainIndex]{$SubstituentPosNum}{Abbrev}[$SubstituentIndex];
      $SubstituentStereoChemistry = $CmpdAbbrevTemplateDataMapRef->{SnSubstituentsPosInfo}[$ChainIndex]{$SubstituentPosNum}{StereoChemistry}[$SubstituentIndex];
      $SubstituentName = $SubstituentAbbrevToNameMap{$SubstituentAbbrev};

      if (exists $SubstituentsDataMap{Abbrev}{$SubstituentName}) {
	$SubstituentsDataMap{Count}{$SubstituentName} += 1;
	$SubstituentsDataMap{PosNumAndStereoPrefix}{$SubstituentName} .= ",${SubstituentPosNum}${SubstituentStereoChemistry}";
      }
      else {
	push @{$SubstituentsDataMap{Name}}, $SubstituentName;
	$SubstituentsDataMap{Abbrev}{$SubstituentName} = $SubstituentAbbrev;
	$SubstituentsDataMap{Count}{$SubstituentName} = 1;
	$SubstituentsDataMap{PosNumAndStereoPrefix}{$SubstituentName} = "${SubstituentPosNum}${SubstituentStereoChemistry}";
	$SubstituentsDataMap{FirstPosNum}{$SubstituentName} = $SubstituentPosNum;
      }
    }
  }
  # Setup name with count prefix...
  my($SubstituentNameWithCount);
  for $SubstituentName (@{$SubstituentsDataMap{Name}}) {
    $SubtituentCount = $SubstituentsDataMap{Count}{$SubstituentName};
    $CountPrefix = $CountToNamePrefixMap{$SubtituentCount};
    $SubstituentNameWithCount = "${CountPrefix}${SubstituentName}";
    $SubstituentsDataMap{NameWithCount}{$SubstituentName} = $SubstituentNameWithCount;
  }

  # Setup a data map based on first position of the substituents and names along with count
  # prefix at that position for further sorting by names...
  my($FirstPosNum, $SubstituentPosNumAndStereoPrefix, %FirstPosNumDataMap);
  %FirstPosNumDataMap = ();
  %{$FirstPosNumDataMap{FirstPosNum}} = ();
  %{$FirstPosNumDataMap{NameWithCount}} = ();
  %{$FirstPosNumDataMap{PosNumAndStereoPrefix}} = ();

  for $SubstituentName (@{$SubstituentsDataMap{Name}}) {
    $FirstPosNum = $SubstituentsDataMap{FirstPosNum}{$SubstituentName};
    $SubstituentNameWithCount = $SubstituentsDataMap{NameWithCount}{$SubstituentName};
    $SubstituentPosNumAndStereoPrefix = $SubstituentsDataMap{PosNumAndStereoPrefix}{$SubstituentName};

    if (exists $FirstPosNumDataMap{FirstPosNum}{$FirstPosNum}) {
      $FirstPosNumDataMap{NameWithCount}{$FirstPosNum}{$SubstituentNameWithCount} = $SubstituentNameWithCount;
      $FirstPosNumDataMap{PosNumAndStereoPrefix}{$FirstPosNum}{$SubstituentNameWithCount} = $SubstituentPosNumAndStereoPrefix;
    }
    else {
      $FirstPosNumDataMap{FirstPosNum}{$FirstPosNum} = $FirstPosNum;
      %{$FirstPosNumDataMap{NameWithCount}{$FirstPosNum}} = ();
      $FirstPosNumDataMap{NameWithCount}{$FirstPosNum}{$SubstituentNameWithCount} = $SubstituentNameWithCount;
      %{$FirstPosNumDataMap{PosNumAndStereoPrefix}{$FirstPosNum}} = ();
      $FirstPosNumDataMap{PosNumAndStereoPrefix}{$FirstPosNum}{$SubstituentNameWithCount} = $SubstituentPosNumAndStereoPrefix;
    }
  }

  # Setup the substituents name...
  for $FirstPosNum (sort {$a <=> $b} keys %{$FirstPosNumDataMap{FirstPosNum}}) {
    for $SubstituentNameWithCount (sort keys %{$FirstPosNumDataMap{NameWithCount}{$FirstPosNum}}) {
      $SubstituentPosNumAndStereoPrefix = $FirstPosNumDataMap{PosNumAndStereoPrefix}{$FirstPosNum}{$SubstituentNameWithCount};
      if ($SubstituentsName) {
	$SubstituentsName .= "-${SubstituentPosNumAndStereoPrefix}-${SubstituentNameWithCount}";
      }
      else {
	$SubstituentsName = "${SubstituentPosNumAndStereoPrefix}-${SubstituentNameWithCount}";
      }
    }
  }
  return $SubstituentsName;
}

# Setup systematic names corresponding to all specified chains...
sub SetupChainsSystematicName {
  my($CmpdAbbrevTemplateDataMapRef, $ChainsIndicesRef, $ChainsPositionsRef) = @_;
  my($ChainIndex, $ChainCount, $Abbrev, $ChainAbbrev, $ChainPosition, $ChainName, $ChainLength, $DoubleBondCount, $DoubleBondGeometry, $Substituents, $SubstituentsAbbrev, $SubstituentsName, $CountName, $NameWithCount, $CountNamePrefix, $ChainSysName, $SystematicName, $Count, $ChainsSystematicName, @Abbrevs, @ChainAbbrevs, @SubstituentsAbbrevs, @ChainNames, @SubstituentsNames, %ChainSysNamesDataMap);

  @Abbrevs = ();
  @ChainAbbrevs = ();
  @SubstituentsAbbrevs = ();

  @ChainNames = ();
  @SubstituentsNames = ();

  for $ChainIndex (@{$ChainsIndicesRef}) {
    $Abbrev = $CmpdAbbrevTemplateDataMapRef->{SnChainAdd}[$ChainIndex] ? $CmpdAbbrevTemplateDataMapRef->{SnAbbrev}[$ChainIndex] : '0:0';
    push @Abbrevs, $Abbrev;

    # Setup chain and substituent abbreviations...
    $ChainAbbrev = $Abbrev;
    $SubstituentsAbbrev = '';

    if ($Abbrev !~ /^0:0$/) {
      ($ChainLength, $DoubleBondCount, $DoubleBondGeometry, $Substituents) = ParseChainAbbrev($Abbrev);
      if (!IsEmpty($Substituents)) {
	$ChainAbbrev = IsEmpty($DoubleBondGeometry) ? "${ChainLength}:${DoubleBondCount}" : "${ChainLength}:${DoubleBondCount}(${DoubleBondGeometry})";
	$SubstituentsAbbrev = $Substituents;
      }
    }
    push @ChainAbbrevs, $ChainAbbrev;
    push @SubstituentsAbbrevs, $SubstituentsAbbrev;

    # Setup chain and substituents names...
    $ChainName = SetupChainName($CmpdAbbrevTemplateDataMapRef, $ChainIndex);
    $SubstituentsName = SetupChainSubstituentsName($CmpdAbbrevTemplateDataMapRef, $ChainIndex);

    push @ChainNames, $ChainName;
    push @SubstituentsNames, $SubstituentsName;
  }


  %ChainSysNamesDataMap = ();
  @{$ChainSysNamesDataMap{Name}} = ();
  %{$ChainSysNamesDataMap{Count}} = ();
  %{$ChainSysNamesDataMap{CountPrefix}} = ();
  %{$ChainSysNamesDataMap{NameWithCount}} = ();

  CHAININDEX: for $ChainIndex (@{$ChainsIndicesRef}) {
    $ChainPosition = (defined $ChainsPositionsRef) ? $ChainsPositionsRef->[$ChainIndex]: ($ChainIndex + 1);

    $ChainName = $ChainNames[$ChainIndex];
    if (IsEmpty($ChainName)) {
      next CHAININDEX;
    }

    $ChainSysName = $ChainName;

    $SubstituentsName = $SubstituentsNames[$ChainIndex];
    if (!IsEmpty($SubstituentsName)) {
      $ChainSysName =~ s/^\(//g;
      $ChainSysName =~ s/\)$//g;
      $ChainSysName = "($SubstituentsName-$ChainSysName)";
    }

    if (exists $ChainSysNamesDataMap{Count}{$ChainSysName}) {
      $ChainSysNamesDataMap{Count}{$ChainSysName} += 1;
      $ChainSysNamesDataMap{CountPrefix}{$ChainSysName} .= ",$ChainPosition";
    }
    else {
      push @{$ChainSysNamesDataMap{Name}}, $ChainSysName;
      $ChainSysNamesDataMap{Count}{$ChainSysName} = 1;
      $ChainSysNamesDataMap{CountPrefix}{$ChainSysName} = $ChainPosition;
    }
  }

  #  Setup name with count...
  for $ChainSysName (@{$ChainSysNamesDataMap{Name}}) {
    $Count = $ChainSysNamesDataMap{Count}{$ChainSysName};

    if ($Count > 1) {
      $CountNamePrefix = exists $CountToNamePrefixMap{$Count} ? $CountToNamePrefixMap{$Count} : '';
      $CountName = "$ChainSysNamesDataMap{CountPrefix}{$ChainSysName}-${CountNamePrefix}";

      $ChainSysNamesDataMap{NameWithCount}{$ChainSysName} = ($ChainSysName =~ /^\(/) ? "${CountName}-${ChainSysName}" : "${CountName}${ChainSysName}";

    }
    else {
      # Just the chain position...
      $CountName = "$ChainSysNamesDataMap{CountPrefix}{$ChainSysName}";

      $ChainSysNamesDataMap{NameWithCount}{$ChainSysName} = "${CountName}-${ChainSysName}";
    }
  }

  # Setup systematic name...
  $ChainsSystematicName = '';
  for $ChainSysName (@{$ChainSysNamesDataMap{Name}}) {
    $NameWithCount = $ChainSysNamesDataMap{NameWithCount}{$ChainSysName};

    if ($ChainsSystematicName) {
      $ChainsSystematicName .= "-$NameWithCount";
    }
    else {
      $ChainsSystematicName = "$NameWithCount";
    }
  }

  return $ChainsSystematicName;
}

# Setup abbreviations corresponding to all chains...
sub SetupChainsAbbrev {
  my($CmpdAbbrevTemplateDataMapRef, $ChainsIndicesRef) = @_;
  my($ChainIndex, $Abbrev, $ChainsAbbrev, @Abbrevs);

  $ChainsAbbrev = '';
  for $ChainIndex (@{$ChainsIndicesRef}) {
    $Abbrev = $CmpdAbbrevTemplateDataMapRef->{SnChainAdd}[$ChainIndex] ? $CmpdAbbrevTemplateDataMapRef->{SnAbbrev}[$ChainIndex] : '0:0';
    push @Abbrevs, $Abbrev;
  }
  $ChainsAbbrev = join("/", @Abbrevs);

  return $ChainsAbbrev;
}


# Initialize data...
sub _InitializeData {
  _InitializeChainAbbrevData();
  _InitializeChainLenAbbrevDbleBondGeometryData();
}

# Setup chain abbreviation to name map...
sub _InitializeChainAbbrevData {
  %ChainAbbrevToNameMap = ();

  # In order to facilitate chain abbreviation validation, ChainAbbrevToNameMap contains multiple
  # enrites for Alkyl and Alkenyl chains using different formats:
  #
  # Alkyl ether chain: 16:0e, O-16:0
  # Alkenyl ether chain: 16:0p, P-16:0, O-16:1(1Z), O-16:1(9Z)
  #
  # Internally, 16:0e is mapped to O-16:0; both 16:0p, P-16:0 get mapped to O-16:1(1Z) for structure
  # generation based on chain lengths and explicit bond specifications.
  #
  # During wild card expansion, GetSupportedChainLenList() function takes out the duplicates by
  # ignoring chain lengths starting with P- and O- and only using the ones containing "e" and "p"
  # in their double bond geomtry specifications.
  #
  %ChainAbbrevToNameMap = (
			   '2:0' => 'acetyl',
			   '3:0' => 'propionyl',
			   '4:0' => 'butyryl',
			   '5:0' => 'valeryl',
			   '6:0' => 'hexanoyl',
			   '7:0' => 'heptanoyl',
			   '8:0' => 'octanoyl',
			   '9:0' => 'nonanoyl',
			   '10:0' => 'decanoyl',
			   'O-10:0' => 'tetradecyl', # EDIT: 27 May 2015 : Added by Chakravarhy Marella
			   '10:1(7Z)' => '(7Z-decanoyl1)', # EDIT : 27 Nov 2012 : Added by Chakravarthy Marella
			   '11:1(7Z)' => '(7Z-decanoyl2)', # EDIT : 27 Nov 2012 : Added by Chakravarthy Marella
			   '12:1(7Z)' => '(7Z-decanoyl3)', # EDIT : 27 Nov 2012 : Added by Chakravarthy Marella
			   '13:1(7Z)' => '(7Z-decanoyl4)', # EDIT : 27 Nov 2012 : Added by Chakravarthy Marella
			   '15:1(7Z)' => '(7Z-decanoyl5)', # EDIT : 27 Nov 2012 : Added by Chakravarthy Marella
			   '17:1(7Z)' => '(7Z-decanoyl6)', # EDIT : 27 Nov 2012 : Added by Chakravarthy Marella
			   '19:1(7Z)' => '(7Z-decanoyl7)', # EDIT : 27 Nov 2012 : Added by Chakravarthy Marella
			   '20:1(7Z)' => '(7Z-decanoyl8)', # EDIT : 27 Nov 2012 : Added by Chakravarthy Marella
			   '11:0' => 'undecanoyl',
			   '12:0' => 'dodecanoyl',
			   'O-12:0' => 'dodecanoyl', # EDIT: 27 May 2015 : Added by Chakravarhy Marella
			   '12:1(9Z)' => '(9Z-dodecanoyl)', # EDIT : 10 Oct 2012 : Added by Chakravarthy Marella
			   'O-12:1(9Z)' => '(9Z-dodecanoyl)', # EDIT : 10 Oct 2012 : Added by Chakravarthy Marella
			   '13:0' => 'tridecanoyl',
			   'O-13:0' => 'tridecanoyl', # EDIT: 27 May 2015 : Added by Chakravarhy Marella
			   '13:1(9Z)' => '(9Z-tridecanoyl)', # EDIT: 27 May 2015 : Added by Chakravarhy Marella
			   'O-13:1(9Z)' => '(9Z-tridecanoyl)', # EDIT: 27 May 2015 : Added by Chakravarhy Marella
			   '14:0' => 'tetradecanoyl',
			   '14:0e' => 'tetradecyl',
			   'O-14:0' => 'tetradecyl',
			   '14:1(9Z)' => '(9Z-tetradecenoyl)',
			   'O-14:1(9Z)' => '(9Z-tetradecenoyl)',# Added on 27.05.215, CM
			   '15:0' => 'pentadecanoyl',
			   'O-15:0' => 'pentadecanoyl', # Added on 27.05.215, CM
			   '15:1(9Z)' => '(9Z-pentadecenoyl)',
			   'O-15:1(9Z)' => '(9Z-pentadecenoyl)', # Added on 27.05.215, CM
			   '16:0' => 'hexadecanoyl',
			   '16:0e' => 'hexadecyl',
			   'O-16:0' => 'hexadecyl',
			   '16:0p' => '(1Z-hexadecenyl)',
			   'P-16:0' => '(1Z-hexadecenyl)',
			   'O-16:1(1Z)' => '(1Z-hexadecenyl)',
			   '16:1e(9Z)' => '(9Z-hexadecenyl)',
			   'O-16:1(9Z)' => '(9Z-hexadecenyl)',
			   '16:1e(11Z)' => '(11Z-hexadecenyl)',
			   'O-16:1(11Z)' => '(11Z-hexadecenyl)',
			   '16:1(7Z)' => '(7Z-hexadecenoyl)',
			   '16:1(9Z)' => '(9Z-hexadecenoyl)',
			   '16:2e(1Z,9Z)' => '(1Z,9Z-hexadecadienyl)',
			   '16:2(9Z,12Z)' => '(9Z,12Z-hexadecadienoyl)', # EDIT : 10 Oct 2012 : Added by Chakravarthy Marella
			   'O-16:2(9Z,12Z)' => '(9Z,12Z-hexadecadienoyl)', # EDIT : 10 Oct 2012 : Added by Chakravarthy Marella
			   'O-16:2(1Z,9Z)' => '(1Z,9Z-hexadecadienyl)',
			   '16:2e(1Z,11Z)' => '(1Z,11Z-hexadecadienyl)',
			   'O-16:2(1Z,11Z)' => '(1Z,11Z-hexadecadienyl)',
			   'P-16:1(9Z)' => '(1Z,9Z-hexadecadienyl)',
			   'P-16:1(11Z)' => '(1Z,11Z-hexadecadienyl)',
			   '17:0' => 'heptadecanoyl',
			   'O-17:0' => 'heptadecanoyl', # Added on 27.05.2015, CM
			   '17:1(9Z)' => '(9Z-heptadecenoyl)',
			   'O-17:1(9Z)' => '(9Z-heptadecenoyl)',# Added on 27.05.2015, CM
			   '17:2(9Z,12Z)' => '(9Z,12Z-heptadecadienoyl)',
			   'O-17:2(9Z,12Z)' => '(9Z,12Z-heptadecadienoyl)',# Added on 27.05.2015, CM
			   '18:0' => 'octadecanoyl',
			   '18:0e' => 'octadecyl',
			   'O-18:0' => 'octadecyl',
			   '18:0p' => '(1Z-octadecenyl)',
			   'P-18:0' => '(1Z-octadecenyl)',
			   'O-18:1(1Z)' => '(1Z-octadecenyl)',
			   '18:1e(9Z)' => '(9Z-octadecenyl)',
			   'O-18:1(9Z)' => '(9Z-octadecenyl)',
			   '18:1e(11Z)' => '(11Z-octadecenyl)',
			   'O-18:1(11Z)' => '(11Z-octadecenyl)',
			   '18:1(4E)' => '(4E-octadecenoyl)',
			   '18:1(6Z)' => '(6Z-octadecenoyl)',
			   '18:1(7Z)' => '(7Z-octadecenoyl)',
			   '18:1(9Z)' => '(9Z-octadecenoyl)',
			   'O-18:1(9Z)' => '(9Z-octadecenoyl)', # Added on 27.05.2015, CM
			   '18:1(9E)' => '(9E-octadecenoyl)',
			   '18:1(11E)' => '(11E-octadecenoyl)',
			   '18:1(11Z)' => '(11Z-octadecenoyl)',
			   '18:1(13Z)' => '(13Z-octadecenoyl)',
			   '18:1(17Z)' => '(13Z-octadecenoyl)',
			   '18:2e(1Z,9Z)' => '(1Z,9Z-octadecadienyl)',
			   'O-18:2(1Z,9Z)' => '(1Z,9Z-octadecadienyl)',
			   '18:2e(1Z,11Z)' => '(1Z,11Z-octadecadienyl)',
			   'O-18:2(1Z,11Z)' => '(1Z,11Z-octadecadienyl)',
			   'P-18:1(9Z)' => '(1Z,9Z-octadecadienyl)',
			   'P-18:1(11Z)' => '(1Z,11Z-octadecadienyl)',
			   '18:2(2E,4E)' => '(2E,4E-octadecadienoyl)',
			   '18:2(6Z,9Z)' => '(6Z,9Z-octadecadienoyl)',
			   '18:2(9E,11E)' => '(9E,11E-octadecadienoyl)',
			   '18:2(9Z,11Z)' => '(9Z,11Z-octadecadienoyl)',
			   '18:2(9Z,12Z)' => '(9Z,12Z-octadecadienoyl)',
			   'O-18:2(9Z,12Z)' => '(9Z,12Z-octadecadienoyl)', # Added on 27 May 2015, CM
			   '18:2(9E,12E)' => '(9E,12E-octadecadienoyl)',
			   '18:3(6Z,9Z,12Z)' => '(6Z,9Z,12Z-octadecatrienoyl)',
			   '18:3(9Z,12Z,15Z)' => '(9Z,12Z,15Z-octadecatrienoyl)',
			   'O-18:3(9Z,12Z,15Z)' => '(9Z,12Z,15Z-octadecatrienoyl)', # Added on 27 May 2015, CM
			   '18:3(9Z,11Z,13Z)' => '(9Z,11Z,13Z-octadecatrienoyl)', # EDIT : 13 May 2015 : Added by Chakravarthy Marella
			   '18:4(6Z,9Z,12Z,15Z)' => '(6Z,9Z,12Z,15Z-octadecatetraenoyl)',
			   '18:4(9E,11E,13E,15E)' => '(9E,11E,13E,15E-octadecatetraenoyl)',
			   '19:0' => 'nonadecanoyl',
			   'O-19:0' => 'nonadecanoyl', # Added on 27.05.2015, CM
			   '19:1(9Z)' => '(9Z-nonadecanoyl)',	# EDIT : 18 May 2015 : Added by Chakravarthy Marella for Drosophila Lipidome
			   'O-19:1(9Z)' => '(9Z-nonadecanoyl)',	# EDIT : 18 May 2015 : Added by Chakravarthy Marella for Drosophila Lipidome
			   '20:0' => 'eicosanoyl',
			   '20:0e' => 'eicosyl',
			   'O-20:0' => 'eicosyl',
			   '20:0p' => '(1Z-eicosenyl)',
			   'P-20:0' => '(1Z-eicosenyl)',
			   'O-20:1(1Z)' => '(1Z-eicosenyl)',
			   '20:1e(9Z)' => '(9Z-eicosenyl)',
			   'O-20:1(9Z)' => '(9Z-eicosenyl)',
			   '20:1e(11Z)' => '(11Z-eicosenyl)',
			   'O-20:1(11Z)' => '(11Z-eicosenyl)',
			   '20:1(2Z)' => '(2Z-eicosenyl)', # EDIT : 12 May 2012 : Added by Chakravarthy Marella for Drosophila Lipidome
			   '20:1(9Z)' => '(9Z-eicosenoyl)', # EDIT : 23 Apr 2015 : Added by Chakravarthy Marella
			   'O-20:1(9Z)' => '(9Z-eicosenoyl)', # EDIT : 23 Apr 2015 : Added by Chakravarthy Marella
			   '20:1(11E)' => '(11E-eicosenoyl)',
			   '20:1(13Z)' => '(13Z-eicosenoyl)',
			   '20:1(13E)' => '(13E-eicosenoyl)',
			   '20:2e(1Z,9Z)' => '(1Z,9Z-eicosadienyl)',
			   'O-20:2(1Z,9Z)' => '(1Z,9Z-eicosadienyl)',
			   '20:2e(1Z,11Z)' => '(1Z,11Z-eicosadienyl)',
			   'O-20:2(1Z,11Z)' => '(1Z,11Z-eicosadienyl)',
			   'P-20:1(9Z)' => '(1Z,9Z-eicosadienyl)',
			   'P-20:1(11Z)' => '(1Z,11Z-eicosadienyl)',
			   '20:2(9Z,12Z)' => '(9Z,11Z-eicosadienoyl)', # EDIT : 18 May 2015 : Added by Chakravarthy Marella
			   'O-20:2(9Z,12Z)' => '(9Z,11Z-eicosadienoyl)', # EDIT : 18 May 2015 : Added by Chakravarthy Marella
			   '20:2(11Z,14Z)' => '(11Z,14Z-eicosadienoyl)',
			   '20:2(5Z,8Z)' => '(5Z,8Z-eicosadienoyl)',
			   '20:3(9Z,12Z,15Z)' => '(9Z,11Z,13Z-eicosatrienoyl)', # EDIT : 18 May 2015 : Added by Chakravarthy Marella
			   'O-20:3(9Z,12Z,15Z)' => '(9Z,11Z,13Z-eicosatrienoyl)', # EDIT : 18 May 2015 : Added by Chakravarthy Marella
			   '20:3(8Z,11Z,14Z)' => '(8Z,11Z,14Z-eicosatrienoyl)',
			   '20:3(5Z,8Z,11Z)' => '(5Z,8Z,11Z-eicosatrienoyl)',
			   '20:4(5Z,8Z,11Z,13E)' => '(5Z,8Z,11Z,13E-eicosatetraenoyl)',
			   '20:4(5Z,8Z,11Z,14Z)' => '(5Z,8Z,11Z,14Z-eicosatetraenoyl)',
			   '20:4(5Z,8Z,10E,14Z)' => '(5Z,8Z,10E,14Z-eicosatetraenoyl)',
			   '20:4(5E,8E,11E,14E)' => '(5E,8E,11E,14E-eicosatetraenoyl)',
			   '20:4(6E,8Z,11Z,14Z)' => '(6E,8Z,11Z,14Z-eicosatetraenoyl)',
			   '20:4(7E,10E,13E,16E)' => '(7E,10E,13E,16E-eicosatetraenoyl)',
			   '20:5(5Z,8Z,11Z,14Z,17Z)' => '(5Z,8Z,11Z,14Z,17Z-eicosapentaenoyl)',
			   '21:0' => 'heneicosanoyl',
			   '22:0' => 'docosenyl',
			   '22:0e' => 'docosenyl',
			   'O-22:0' => 'docosenyl',
			   '22:0p' => '(1Z-docosenyl)',
			   'P-22:0' => '(1Z-docosenyl)',
			   'O-22:1(1Z)' => '(1Z-docosenyl)',
			   '22:1e(9Z)' => '(9Z-docosenyl)',
			   'O-22:1(9Z)' => '(9Z-docosenyl)',
			   '22:1e(11Z)' => '(11Z-docosenyl)',
			   'O-22:1(11Z)' => '(11Z-docosenyl)',
			   '22:1(9Z)' => '(9Z-docosenoyl)', # Added on 18 May 2015 by Chakravarthy Marella
			   '22:1(13Z)' => '(13Z-docosenoyl)',
			   '22:2e(1Z,9Z)' => '(1Z,9Z-docosenyl)',
			   'O-22:2(1Z,9Z)' => '(1Z,9Z-docosenyl)',
			   '22:2e(1Z,11Z)' => '(1Z,11Z-docosenyl)',
			   'O-22:2(1Z,11Z)' => '(1Z,11Z-docosenyl)',
			   'P-22:1(9Z)' => '(1Z,9Z-docosenyl)',
			   'P-22:1(11Z)' => '(1Z,11Z-docosenyl)',
			   '22:2(9Z,12Z)' => '(9Z,12Z-docosadienoyl)', # EDIT : 18 May 2015 : Added by Chakravarthy Marella
			   'O-22:2(9Z,12Z)' => '(9Z,12Z-docosadienoyl)', # EDIT : 18 May 2015 : Added by Chakravarthy Marella
			   '22:2(13Z,16Z)' => '(13Z,16Z-docosadienoyl)',
			   '22:3(9Z,12Z,15Z)' => '(9Z,12Z,15Z-docosadienoyl)',# EDIT : 18 May 2015 : Added by Chakravarthy Marella
			   'O-22:3(9Z,12Z,15Z)' => '(9Z,12Z,15Z-docosadienoyl)',# EDIT : 18 May 2015 : Added by Chakravarthy Marella
			   '22:4(7Z,10Z,13Z,16Z)' => '(7Z,10Z,13Z,16Z-docosatetraenoyl)',
			   '22:5(4Z,7Z,10Z,13Z,16Z)' => '(4Z,7Z,10Z,13Z,16Z-docosapentaenoyl)',
			   '22:5(7Z,10Z,13Z,16Z,19Z)' => '(7Z,10Z,13Z,16Z,19Z-docosapentaenoyl)',
			   '22:6(4Z,7Z,10Z,12E,16Z,19Z)' => '(4Z,7Z,10Z,12E,16Z,19Z-docosahexaenoyl)',
			   '22:6(4Z,7Z,10Z,13Z,16Z,19Z)' => '(4Z,7Z,10Z,13Z,16Z,19Z-docosahexaenoyl)',
			   '23:0' => 'tricosanoyl',
			   '24:0' => 'tetracosanoyl',
			   '24:1(9Z)' => '(9Z-tetracosenoyl)', # Added on 18 May 2015 by Chakravarthy Marella
			   'O-24:1(9Z)' => '(9Z-tetracosenoyl)', # Added on 18 May 2015 by Chakravarthy Marella
			   '24:1(15Z)' => '(15Z-tetracosenoyl)',
			   '24:2(9Z,12Z)' => '(9Z,12Z-tetracosenoyl)', # Added on 18 May 2015 by Chakravarthy Marella
			   'O-24:2(9Z,12Z)' => '(9Z,12Z-tetracosenoyl)', # Added on 18 May 2015 by Chakravarthy Marella
			   '24:3(9Z,12Z,15Z)' => '(9Z,12Z,15Z-tetracosenoyl)', # Added on 18 May 2015 by Chakravarthy Marella
			   'O-24:3(9Z,12Z,15Z)' => '(9Z,12Z,15Z-tetracosenoyl)', # Added on 18 May 2015 by Chakravarthy Marella
			   '24:4(5Z,8Z,11Z,14Z)' => '(5Z,8Z,11Z,14Z-tetracosatetraenoyl)',
			   '25:0' => 'pentacosanoyl',
			   '26:0' => 'hexacosanoyl',
			   '26:1(5Z)' => '(5Z-hexacosenoyl)',
			   '26:1(9Z)' => '(5Z-hexacosenoyl)', # Added on 23 Apr 2015 by Chakravarthy Marella
			   '26:2(5Z,9Z)' => '(5Z,9Z-hexacosadienoyl)',
			   '26:2(5Z,9E)' => '(5Z,9E-hexacosadienoyl)',
			   '26:2(5E,9Z)' => '(5Z,9E-hexacosadienoyl)',
			   '27:0' => 'heptacosanoyl',
			   '28:0' => 'octacosanoyl',
			   '29:0' => 'nonacosanoyl',
			   '30:0' => 'triacontanoyl',
			   '31:0' => 'hentriacontanoyl',
			   '32:0' => 'dotriacontanoyl',
			   '33:0' => 'tritriacontanoyl',
			   '34:0' => 'tetratriacontanoyl',
			   '35:0' => 'pentatriacontanoyl',
			   '36:0' => 'hexatriacontanoyl',
			   '37:0' => 'heptatriacontanoyl',
			   '38:0' => 'octatriacontanoyl',
			   '39:0' => 'nonatriacontanoyl'
			  );

  # Setup substituents map...
  #
  %SubstituentAbbrevToNameMap = ();
  %SubstituentAbbrevToNameMap = (
				 'OH' => 'hydroxy',
				 'NH2' => 'amino',
				 'SH' => 'thio',
				 'Me' => 'methyl',
				 'Et' => 'ethyl',
				 'Pr' => 'propyl',
				 'OMe' => 'methoxy',
				 'OAc' => 'acetoxy',
				 'Ke' => 'oxo',
				 'Ep' => 'epoxy',
				 'Cp' => 'cyclopropyl',
				 'My' => 'methylene',
				 'OOH' => 'hydroperoxy',
				 'Br' => 'bromo',
				 'Cl' => 'chloro',
				 'F' => 'fluro',
				 'I' => 'iodo',
				 'CN' => 'cyno',
				 'NO2' => 'nitro',
				 'COOH' => 'carboxy',
				 'CHO' => 'aldehyde'
				);

  # Setup  chain length to name prefix map for dynamic generation of systematic names...
  #
  %ChainLenToNamePrefixMap = ();
  %ChainLenToNamePrefixMap = (
			      '2' => 'eth',
			      '3' => 'prop',
			      '4' => 'but',
			      '5' => 'pent',
			      '6' => 'hex',
			      '7' => 'hept',
			      '8' => 'oct',
			      '9' => 'non',
			      '10' => 'dec',
			      '11' => 'undec',
			      '12' => 'dodec',
			      '13' => 'tridec',
			      '14' => 'tetradec',
			      '15' => 'pentadec',
			      '16' => 'hexadec',
			      '17' => 'heptadec',
			      '18' => 'octadec',
			      '19' => 'nonadec',
			      '20' => 'eicos',
			      '21' => 'heneicos',
			      '22' => 'docos',
			      '23' => 'tricos',
			      '24' => 'tetracos',
			      '25' => 'pentacos',
			      '26' => 'hexacos',
			      '27' => 'heptacos',
			      '28' => 'octacos',
			      '29' => 'nonacos',
			      '30' => 'triacont',
			      '31' => 'hentriacont',
			      '32' => 'dotriacont',
			      '33' => 'tritriacont',
			      '34' => 'tetratriacont',
			      '35' => 'pentatriacont',
			      '36' => 'hexatriacont',
			      '37' => 'heptatriacont',
			      '38' => 'octatriacont',
			      '39' => 'nonatriacont',
			      '40' => 'tetracont',
			      '41' => 'hentriacont',
			      '42' => 'dotetracont',
			      '43' => 'tritetracont',
			      '44' => 'tetratetracont',
			      '45' => 'pentatetracont',
			      '46' => 'hexatetracont',
			      '47' => 'heptatetracont',
			      '48' => 'octatetracont',
			      '49' => 'nonatetracont',
			      '50' => 'pentacont'
			     );

  # Setup  to name prefix map for dynamic generation of systematic names for compounds
  # involving multiple bonds and substituents...
  %CountToNamePrefixMap = ();
  %CountToNamePrefixMap = (
			   '1' => '',
			   '2' => 'di',
			   '3' => 'tri',
			   '4' => 'tetra',
			   '5' => 'penta',
			   '6' => 'hexa',
			   '7' => 'hepta',
			   '8' => 'octa',
			   '9' => 'nona',
			   '10' => 'deca',
			   '11' => 'undeca',
			   '12' => 'dodeca',
			   '13' => 'trideca',
			   '14' => 'tetradeca',
			   '15' => 'pentadeca',
			   '16' => 'hexadeca',
			   '17' => 'heptadeca',
			   '18' => 'octadeca',
			   '19' => 'nonadeca',
			   '20' => 'eicosa'
			  );
}

# To facilitate expansion of wild card specs - 18:1(*), 18:* and so on, set up
#  ChainLenToDbleBondGeomtryMap
sub _InitializeChainLenAbbrevDbleBondGeometryData {
  my($ChainAbbrev, $ChainLength, $DoubleBondSpec, $DoubleBondCount, $DoubleBondGeometry);

  %ChainLenAbbrevSupportedMap = ();
  %ChainLenAbbrevToDbleBondGeomtryMap = ();

  for $ChainAbbrev (sort keys %ChainAbbrevToNameMap) {
    ($ChainLength, $DoubleBondCount, $DoubleBondGeometry) = ParseChainAbbrev($ChainAbbrev);
    if (exists $ChainLenAbbrevSupportedMap{$ChainLength}) {
      if (exists $ChainLenAbbrevToDbleBondGeomtryMap{$ChainLength}{$DoubleBondCount}) {
	$ChainLenAbbrevToDbleBondGeomtryMap{$ChainLength}{$DoubleBondCount}{$DoubleBondGeometry} = $DoubleBondGeometry;
      }
      else {
	$ChainLenAbbrevToDbleBondGeomtryMap{$ChainLength}{$DoubleBondCount}{$DoubleBondGeometry} = $DoubleBondGeometry;
      }
    }
    else {
      $ChainLenAbbrevSupportedMap{$ChainLength} = $ChainLength;
      $ChainLenAbbrevToDbleBondGeomtryMap{$ChainLength}{$DoubleBondCount}{$DoubleBondGeometry} = $DoubleBondGeometry;
    }
  }

}

1;

__END__

=head1 NAME

ChainAbbrev - Methods for processing chain abbreviations

=head1 SYNOPSIS

use ChainAbbrev;

use ChainAbbrev qw(:all);

=head1 DESCRIPTION

ChainAbbrev module provides these methods:

    ChainAbbrevNameExists - Is it a supported chain abbreviation
    ExpandChainAbbrev - Expand wild cards in chain abbreviation
    GetChainAbbrevToNameMap - Get chain name
    GetChainLenAbbrevSupportedMap - Get reference to supported chain
                                    abbreviations data
    GetChainLenAbbrevDbleBondGeometryDataMap - Get reference to supported
                                               double bond geometry data
    GetChainLengthAndMultipleBondCount - Get chain length and number of
                                         double and triple bonds
    GetChainLenToNamePrefixMap - Get chain name prefix
    GetCountToNamePrefixMap - Get count prefix
    GetSubstituentsAbbrevToNameMap - Get substituents name
    GetSubstituentBondOrder - Get substituent bond order
    GetSupportedChainLenList - Get supported chain lengths
    IsAlkylChainAbbrev - Is it a alkyl chain abbreviation
    IsAlkenylChainAbbrev - Is it a alkenyl chain abbreviation
    IsChainAbbrevOkay - Is it a valid chain abbreviation
    IsDoubleBondsAbbrevOkay - Is it a valid double bond abbreviation
    IsRingsAbbrevOkay - Is it a valid ring abbreviation
    IsSubstituentsAbbrevOkay - Is it a valid substituent abbreviation
    IsWildCardInChainAbbrev - Does chain abbreviation contains a wild card
    ParseChainAbbrev - Parse chain abbreviation
    ParseRingAbbrev - Parse ring abbreviation
    ParseSubstituentAbbrev - Parse substituent abbreviation
    SetupChainSubstituentsName - Set up substituent name

=head1 METHODS

=over 4

=item B<ChainAbbrevNameExists>

    $Status = ChainAbbrevNameExists($ChainAbbrev);

Return 1 or 0 based on whether it's a supported chain name.

=item B<ExpandChainAbbrev>

    $AbbrevArrayRef = ExpandChainAbbrev($Abbrev);

Return a reference to an array containing complete chain abbreviations. Wild card
characters in chain abbreviation name are expanded to generate fully qualified
chain abbreviations.

=item B<GetChainAbbrevToNameMap>

    $AbbrevNameHashRef = GetChainAbbrevToNameMap();

Return a reference to hash with chain abbreviation/name as key/value pair.

=item B<GetChainLenAbbrevSupportedMap>

    $ChainLenHashRef = GetChainLenAbbrevSupportedMap();

Return a reference to hash with supported chain length as hash key.

=item B<GetChainLenAbbrevDbleBondGeometryDataMap>

    $ChainLenDblBondHashRef =
        GetChainLenAbbrevDbleBondGeometryDataMap();

Return a reference to hash containing information about chain length, number of
double bonds and geometry of double bonds.

=item B<GetChainLengthAndMultipleBondCount>

    ($ChainLength, $DoubleBondCount, $TripleBondCount) =
        GetChainLengthAndMultipleBondCount($ChainAbbrev);

Parse chain abbreviation and return these values: chain length; number of
double and triple bonds.

=item B<GetChainLenToNamePrefixMap>

    $ChainNameHashRef = GetChainLenToNamePrefixMap();

Return a reference to hash with chain length/name prefix as key/value pair.

=item B<GetCountToNamePrefixMap>

    $CountHashRef = GetCountToNamePrefixMap();

Return a reference to hash with count/name prefix as key/value pair.

=item B<GetSubstituentsAbbrevToNameMap>

    $AbbrevNameHashRef = GetSubstituentsAbbrevToNameMap();

Return a reference to hash with substituents abbreviation/name as key/value pair.

=item B<GetSubstituentBondOrder>

    $BondOrder = GetSubstituentBondOrder($SubstituentAbbrev);

Return bond order for a sustituent.

=item B<GetSupportedChainLenList>

    $ChainLengthListRef = GetSupportedChainLenList();

Return a reference to a sorted list containing supported chain lengths.

=item B<IsAlkylChainAbbrev>

    $Status = IsAlkylChainAbbrev($ChainAbbrev);

Return 1 or 0 based on whether it's a alkyl chain abbreviation.

=item B<IsAlkenylChainAbbrev>

    $Status = IsAlkenylChainAbbrev($ChainAbbrev);

Return 1 or 0 based on whether it's a alkenyl chain abbreviation.

=item B<IsChainAbbrevOkay>

    $Status = IsChainAbbrevOkay($ChainAbbrev);

Return 1 or 0 based on whether chain abbreviation is valid.

=item B<IsDoubleBondsAbbrevOkay>

    $Status = IsDoubleBondsAbbrevOkay($ChainAbbrev, $ChainLength,
        $DoubleBondCount, $DoubleBondGeometry);

Return 1 or 0 based on whether chain abbreviation contains a valid multiple bond specification.

=item B<IsRingsAbbrevOkay>

    $Status = IsRingsAbbrevOkay($ChainAbbrev, $ChainLength, $Rings);

Return 1 or 0 based on whether chain abbreviation contains a valid ring specification.

=item B<IsSubstituentsAbbrevOkay>

    $Status = IsSubstituentsAbbrevOkay($ChainAbbrev, $ChainLength,
        $DoubleBondCount, $DoubleBondGeometry, $Substituents);

Return 1 or 0 based on whether chain abbreviation contains a valid substituents specification.

=item B<IsWildCardInChainAbbrev>

    $Status = IsWildCardInChainAbbrev($ChainAbbrev);

Return 1 or 0 based on whether chain abbreviation contains any wild card character.

=item B<ParseChainAbbrev>

    ($ChainLength, $DoubleBondCount, $DoubleBondGeometry) =
        ParseChainAbbrev($ChainAbbrev);

Parse chain abbreviation and return these values: chain length, number of double bonds,
and geometry of double bonds.

=item B<ParseRingAbbrev>

    ($Pos, $StereoChemistry) = ParseRingAbbrev($ChainAbbrev);

Parse chain abbreviation and return these values: ring position and stereochemistry
specificaton at the ring.

=item B<ParseSubstituentAbbrev>

    ($Pos, $Name, $StereoChemistry) =
        ParseSubstituentAbbrev($SubstituentAbbrev);

Parse substituent abbreviation and return these values: position of the substituent
on the chain, name, and stereochemistry of the substituent.

=item B<SetupChainSubstituentsName>

    $SubstituentsName = SetupChainSubstituentsName(
        $CmpdAbbrevTemplateDataMapRef, $ChainIndex);

Return systematic name for substituents after ordering and grouping substituents by their
position.

=back

=head1 AUTHOR

Manish Sud

=head1 CONTRIBUTOR

Eoin Fahy

=head1 SEE ALSO

LMAPSStr.pm, ChainStr.pm

=head1 COPYRIGHT

Copyright (C) 2006-2012. The Regents of the University of California. All Rights Reserved.

=head1 LICENSE

Modified BSD License

=cut
