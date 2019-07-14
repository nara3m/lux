package SPChainAbbrev;
#
# File: SPChainAbbrev.pm
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
use TextUtil;
use ChainAbbrev;

our ($VERSION, @ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);

$VERSION = '2.00';
@ISA = qw(Exporter);
@EXPORT = qw();
@EXPORT_OK = qw(ChainAbbrevNameExists ExpandChainAbbrev GetChainAbbrevToNameMap GetChainLenAbbrevSupportedMap GetChainLenAbbrevDbleBondGeometyDataMap GetSupportedChainLenList IsChainAbbrevOkay IsSphingosineChainAbbrev IsSphingosineC18ChainAbbrev IsSphinganineC18ChainAbbrev);
%EXPORT_TAGS = (
		all  => [@EXPORT, @EXPORT_OK],
		);

# Initialize data...
my(%Sn1ChainAbbrevToNameMap, %Sn1ChainLenAbbrevSupportedMap, %Sn1ChainLenAbbrevToDbleBondGeomtryMap, %Sn2ChainAbbrevToNameMap, %Sn2ChainLenAbbrevSupportedMap, %Sn2ChainLenAbbrevToDbleBondGeomtryMap);
_InitializeData();

# Expand chain abbreviation by enumerating all possibilites for wild cards
# and return a reference to expanded list...
sub ExpandChainAbbrev {
  my($Abbrev, $ChainType, @ExpandedAbbrevs);

  ($Abbrev, $ChainType) = @_;
  @ExpandedAbbrevs = ();
  if (!ChainAbbrev::IsWildCardInChainAbbrev($Abbrev)) {
    push @ExpandedAbbrevs, $Abbrev;
    return \@ExpandedAbbrevs;
  }
  # Expand wild cards...
  my($NewAbbrev, $ChainLengthAbbrev, $DoubleBondCountAbbrev, $DoubleBondGeometryAbbrev, $SubstituentsAbbrev, $ChainLength, $DoubleBondCount, $DoubleBondGeometry, @ChainList, @DoubleBondCountList, @DoubleBondGeometryList);

  ($ChainLengthAbbrev, $DoubleBondCountAbbrev, $DoubleBondGeometryAbbrev, $SubstituentsAbbrev) = ChainAbbrev::ParseChainAbbrev($Abbrev);

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
    push @ChainList, GetSupportedChainLenList($ChainType);
  }
  my($ChainLenAbbrevToDbleBondGeomtryRef);
  $ChainLenAbbrevToDbleBondGeomtryRef = GetChainLenAbbrevDbleBondGeometyDataMap($ChainType);

  # Chain length  loop...
  CHAINLEN: for $ChainLength (@ChainList) {
    if (!ChainAbbrev::IsChainLengthOkay($ChainLength, $ChainLengthAbbrev)) {
      next CHAINLEN;
    }
    @DoubleBondCountList = ($DoubleBondCountAbbrev =~ /\*/) ? (sort { $a cmp $b } keys %{$ChainLenAbbrevToDbleBondGeomtryRef->{$ChainLength}}) : $DoubleBondCountAbbrev;

    # Double bond count loop...
    DOUBLEBONDCOUNT: for $DoubleBondCount (@DoubleBondCountList) {
      if ($DoubleBondCount && !(exists $ChainLenAbbrevToDbleBondGeomtryRef->{$ChainLength}{$DoubleBondCount})) {
	next DOUBLEBONDCOUNT;
      }
      @DoubleBondGeometryList = ($DoubleBondGeometryAbbrev =~ /\*/) ? (sort keys %{$ChainLenAbbrevToDbleBondGeomtryRef->{$ChainLength}{$DoubleBondCount}}) :  $DoubleBondGeometryAbbrev;

      # Double bond geomerty loop...
      DOUBLEBONDGEOMERTY: for $DoubleBondGeometry (@DoubleBondGeometryList) {
	if ($DoubleBondGeometry && !(exists $ChainLenAbbrevToDbleBondGeomtryRef->{$ChainLength}{$DoubleBondCount}{$DoubleBondGeometry})) {
	  next DOUBLEBONDGEOMERTY;
	}

	$NewAbbrev = ($DoubleBondCount && $DoubleBondGeometry) ? ("$ChainLength:$DoubleBondCount($DoubleBondGeometry)") : ("$ChainLength:$DoubleBondCount");
	if (!IsEmpty($SubstituentsAbbrev)) {
	  $NewAbbrev .= "(${SubstituentsAbbrev})";
	}
	push @ExpandedAbbrevs, $NewAbbrev;
      }

    }
  }

  return \@ExpandedAbbrevs;
}

# Is this a known chain abbreviation?
#
sub ChainAbbrevNameExists {
  my($ChainAbbrev, $ChainType) = @_;
  my($Abbrev, $ChainLen, $BondCount, $BondGeomety, $Substituents, $ChainAbbrevToNameRef);

  $ChainAbbrevToNameRef = GetChainAbbrevToNameMap($ChainType);

  ($ChainLen, $BondCount, $BondGeomety, $Substituents) = ChainAbbrev::ParseChainAbbrev($ChainAbbrev);
  $Abbrev = $BondGeomety ? "${ChainLen}:${BondCount}(${BondGeomety})" : "${ChainLen}:${BondCount}";
  if (exists $ChainAbbrevToNameRef->{$Abbrev}) {
    return 1;
  }
  return 0;
}

# Return a refernce to abbreviation data hash for sn1 or sn2..
sub GetChainAbbrevToNameMap {
  my($ChainType) = @_;

  if ($ChainType =~ /^Sn1$/i) {
    return \%Sn1ChainAbbrevToNameMap;
  }
  elsif ($ChainType =~ /^Sn2$/i) {
    return \%Sn2ChainAbbrevToNameMap;
  }
  else {
    return undef;
  }
}

# Return a refernce to chain length hash for sn1 or sn2..
sub GetChainLenAbbrevSupportedMap {
  my($ChainType) = @_;

  if ($ChainType =~ /^Sn1$/i) {
    return \%Sn1ChainLenAbbrevSupportedMap;
  }
  elsif ($ChainType =~ /^Sn2$/i) {
    return \%Sn2ChainLenAbbrevSupportedMap;
  }
  else {
    return undef;
  }
}

# Return a refernce to double bond geometry hash for sn1 or sn2..
sub GetChainLenAbbrevDbleBondGeometyDataMap {
  my($ChainType) = @_;

  if ($ChainType =~ /^Sn1$/i) {
    return \%Sn1ChainLenAbbrevToDbleBondGeomtryMap;
  }
  elsif ($ChainType =~ /^Sn2$/i) {
    return \%Sn2ChainLenAbbrevToDbleBondGeomtryMap;
  }
  else {
    return undef;
  }
}

# Setup a supported chain length for sn1 or sn2 positions...
#
sub GetSupportedChainLenList {
  my($ChainType) = @_;
  my(@SortedChainLenList, @ChainLenList, $ChainLenAbbrevSupportedRef);

  $ChainLenAbbrevSupportedRef = GetChainLenAbbrevSupportedMap($ChainType);

  @SortedChainLenList = ();
  @ChainLenList = ();

  my($ChainLen);
  CHAINLEN: for $ChainLen (keys %{$ChainLenAbbrevSupportedRef}) {
    push @ChainLenList, $ChainLen;
  }

  # Sort 'em out...
  for $ChainLen (sort {$a <=> $b} @ChainLenList) {
    push @SortedChainLenList, $ChainLen;
  }
  return @SortedChainLenList;
}


# Check out the chain abrbrev...
#
# . Wild cards are not allowed for substituents and rings
#
sub IsChainAbbrevOkay {
  my($ChainAbbrev, $ChainType, $AllowSubstituents, $AllowRings, $AllowArbitraryChainLenSpec) = @_;
  my($ChainLength, $DoubleBondCount, $DoubleBondGeometry, $Substituents, $RingsAbbrev);

  $AllowSubstituents = defined $AllowSubstituents ? $AllowSubstituents : 1;
  $AllowRings = defined $AllowRings ? $AllowRings : 1;
  $AllowArbitraryChainLenSpec = defined $AllowArbitraryChainLenSpec ? $AllowArbitraryChainLenSpec : 0;

  ($ChainLength, $DoubleBondCount, $DoubleBondGeometry, $Substituents, $RingsAbbrev) = ChainAbbrev::ParseChainAbbrev($ChainAbbrev);

  if (!(length $ChainLength)) {
    warn "Warning: Ignoring chain abbreviation $ChainAbbrev : No chain length specified.\n";
    return 0;
  }

  if ($ChainLength =~ /\*/ ) {
    # Assume it's okay and wouble be rechecked during expansion with a valid chain length...
    return 1;
  }

  my($ChainLenAbbrevSupportedRef, $ChainLenAbbrevToDbleBondGeomtryRef);

  $ChainLenAbbrevSupportedRef = GetChainLenAbbrevSupportedMap($ChainType);
  $ChainLenAbbrevToDbleBondGeomtryRef = GetChainLenAbbrevDbleBondGeometyDataMap($ChainType);

  if ($ChainLength && !$AllowArbitraryChainLenSpec) {
    if (!(exists $ChainLenAbbrevSupportedRef->{$ChainLength})) {
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
      if ($ChainAbbrev ne "0:0" && !(ChainAbbrevNameExists($ChainAbbrev, $ChainType))) {
	warn "Warning: Ignoring chain abbreviation $ChainAbbrev : Unknown chain, double bound count, or geometry\n";
	return 0;
      }
    }
    elsif ($NoWildCharInCount) {
      if ($DoubleBondCount == 0 ) {
	if (length $DoubleBondGeometry) {
	  warn "Warning: Ignoring chain abbreviation $ChainAbbrev : Double bond count must be non zero to specify geometry.\n";
	  return 0;
	}
      }
      else {
	if (!(exists $ChainLenAbbrevToDbleBondGeomtryRef->{$ChainLength}{$DoubleBondCount})) {
	  warn "Warning: Ignoring chain abbreviation $ChainAbbrev : Unknown double bond count.\n";
	  return 0;
	}
	if ($NoWildCharInGeometry && length $DoubleBondGeometry) {
	  if (!(exists $ChainLenAbbrevToDbleBondGeomtryRef->{$ChainLength}{$DoubleBondCount}{$DoubleBondGeometry})) {
	    warn "Warning: Ignoring chain abbreviation $ChainAbbrev : Unknown double bond geometry.\n";
	    return 0;
	  }
	}
      }
    }
    else {
      if ($NoWildCharInGeometry && length $DoubleBondGeometry) {
	my(@DoubleBondCountList, $GeometryFound, $BondCount);
	@DoubleBondCountList =  sort keys %{$ChainLenAbbrevToDbleBondGeomtryRef->{$ChainLength}};
	$GeometryFound = 0;
           BONDCOUNT: for $BondCount (@DoubleBondCountList) {
	  if (exists $ChainLenAbbrevToDbleBondGeomtryRef->{$ChainLength}{$BondCount}{$DoubleBondGeometry}) {
	    $GeometryFound = 1;
	    last BONDCOUNT;
	  }
	}
	if (!$GeometryFound) {
	  warn "Warning: Ignoring chain abbreviation $ChainAbbrev : Unknown double bond geometry.\n";
	  return 0;
	}
      }
    }
  }

  if ($Substituents) {
    if (!$AllowSubstituents) {
      warn "Warning: Ignoring chain abbreviation $ChainAbbrev : Substituents are not allowed.\n";
      return 0;
    }
    if (!ChainAbbrev::IsSubstituentsAbbrevOkay($ChainAbbrev, $ChainLength, $DoubleBondCount, $DoubleBondGeometry, $Substituents)) {
      return 0;
    }
  }
  if ($RingsAbbrev) {
    if (!$AllowRings) {
      warn "Warning: Ignoring chain abbreviation $ChainAbbrev : Rings are not allowed\n";
      return 0;
    }
    if (!ChainAbbrev::IsRingsAbbrevOkay($ChainAbbrev, $ChainLength, $RingsAbbrev)) {
      return 0;
    }
  }

  return 1;
}

# Initialize data...
sub _InitializeData {
  _InitializeChainAbbrevData();
  _InitializeChainLenAbbrevDbleBondGeometryData('Sn1');
  _InitializeChainLenAbbrevDbleBondGeometryData('Sn2');
}


# Setup chain abbreviation to name map...
sub _InitializeChainAbbrevData {
  # sn1 chain corresponds to sphingoid base and the chain length includes all the three
  # back bone carbon atoms as well. sphingnine, shinga-4-enine, and 4,8-sphingadiene
  # are represented by <ChainLen>:0, <ChainLen>:1(4E) and <ChainLen>:2(4E,8E)
  # respectively.
  #
  # For chains containing
  #
  %Sn1ChainAbbrevToNameMap = (
			   '14:0' => 'tetradecanoyl',
			   '14:1(4E)' => '4E-tetradecenoyl',
			   '14:1(8E)' => '8E-tetradecenoyl',
			   '14:2(4E,6E)' => '4E,6E-tetradecadienoyl',
			   '14:2(4E,8E)' => '4E,8E-tetradecadienoyl',
			   '15:0' => 'pentadecanoyl',
			   '15:1(4E)' => '4E-pentadecenoyl',
			   '15:1(8E)' => '8E-pentadecenoyl',
			   '15:2(4E,6E)' => '4E,6E-pentadecadienoyl',
			   '15:2(4E,8E)' => '4E,8E-pentadecadienoyl',
			   '16:0' => 'hexadecanoyl',
			   '16:1(4E)' => '4E-hexadecenoyl',
			   '16:1(8E)' => '8E-hexadecenoyl',
			   '16:2(4E,6E)' => '4E,6E-hexadecadienoyl',
			   '16:2(4E,8E)' => '4E,8E-hexadecadienoyl',
			   '17:0' => 'heptadecanoyl',
			   '17:1(4E)' => '4E-heptadecenoyl',
			   '17:1(8E)' => '8E-heptadecenoyl',
			   '17:2(4E,8E)' => '4E,8E-heptadecadienoyl',
			   '18:0' => 'octadecanoyl',
			   '18:1(4E)' => '4E-octadecenoyl',
			   '18:1(8E)' => '8E-octadecenoyl',
			   '18:2(4E,8E)' => '4E,8E-octadecadienoyl',
			   '18:2(4E,14Z)' => '4E,14Z-octadecadienoyl',
			   '19:0' => 'nonadecanoyl',
			   '19:1(4E)' => '4E-nonadecenoyl',
			   '19:1(8E)' => '8E-nonadecenoyl',
			   '19:2(4E,8E)' => '4E,8E-nonadecadienoyl',
			   '20:0' => 'eicosanoyl',
			   '20:1(4E)' => '4E-eicosenoyl',
			   '20:1(8E)' => '8E-eicosenoyl',
			   '20:2(4E,8E)' => '4E,8E-eicosadienoyl',
			   '21:0' => 'heneicosanoyl',
			   '21:1(4E)' => '4E-heneicosenoyl',
			   '21:1(8E)' => '8E-heneicosenoyl',
			   '21:2(4E,8E)' => '4E,8E-heneicosadienoyl',
			   '22:0' => 'docosanoyl',
			   '22:1(4E)' => '4E-docosenoyl',
			   '22:1(8E)' => '8E-docosenoyl',
			   '22:2(4E,8E)' => '4E,8E-docosadienoyl',
			  );

  # sn2 chain corresponds to N-acyl group...
  #
  %Sn2ChainAbbrevToNameMap = (
			   '2:0' => 'acetyl',
			   '10:0' => 'decanoyl',
			   '12:0' => 'dodecanoyl',
			   '13:0' => 'tridecanoyl',
			   '14:0' => 'tetradecanoyl',
			   '15:0' => 'pentadecanoyl',
			   '16:0' => 'hexadecanoyl',
			   '16:1(9Z)' => '9Z-hexadecenoyl',
			   '17:0' => 'heptadecanoyl',
			   '17:1(9Z)' => '9Z-heptadecenoyl',
			   '18:0' => 'octadecanoyl',
			   '18:1(9Z)' => '9Z-octadecenoyl',
			   '18:2(9Z,12Z)' => '9Z,12Z-octadecadienoyl',
			   '19:0' => 'nonadecanoyl',
			   '20:0' => 'eicosanoyl',
			   '20:1(11Z)' => '11Z-eicosenoyl',
			   '21:0' => 'heneicosanoyl',
			   '22:0' => 'docosanoyl',
			   '22:1(13Z)' => '13Z-docosenoyl',
			   '23:0' => 'tricosanoyl',
			   '24:0' => 'tetracosanoyl',
			   '24:1(15Z)' => '15Z-tetracosenoyl',
			   '24:4(5Z,8Z,11Z,14Z)' => '5Z,8Z,11Z,14Z-tetracosatetraenoyl',
			   '25:0' => 'pentacosanoyl',
			   '26:0' => 'hexacosanoyl',
			   '26:1(17Z)' => '17Z-hexacosenoyl',
			   '27:0' => 'heptacosanoyl'
			  );
}

# To facilitate expansion of wild card specs - 18:1(*), 18:* and so on, set up
#  ChainLenToDbleBondGeomtryMap
sub _InitializeChainLenAbbrevDbleBondGeometryData {
  my($ChainType) = @_;

  my($ChainAbbrev, $ChainLength, $DoubleBondSpec, $DoubleBondCount, $DoubleBondGeometry, $ChainAbbrevToNameRef, $ChainLenAbbrevSupportedMapRef, $ChainLenAbbrevToDbleBondGeomtryMapRef);

  if ($ChainType =~ /^Sn1$/i) {
    %Sn1ChainLenAbbrevSupportedMap = ();
    %Sn1ChainLenAbbrevToDbleBondGeomtryMap = ();
    $ChainAbbrevToNameRef = \%Sn1ChainAbbrevToNameMap;
    $ChainLenAbbrevSupportedMapRef = \%Sn1ChainLenAbbrevSupportedMap;
    $ChainLenAbbrevToDbleBondGeomtryMapRef = \%Sn1ChainLenAbbrevToDbleBondGeomtryMap;
  }
  else {
    %Sn2ChainLenAbbrevSupportedMap = ();
    %Sn2ChainLenAbbrevToDbleBondGeomtryMap = ();
    $ChainAbbrevToNameRef = \%Sn2ChainAbbrevToNameMap;
    $ChainLenAbbrevSupportedMapRef = \%Sn2ChainLenAbbrevSupportedMap;
    $ChainLenAbbrevToDbleBondGeomtryMapRef = \%Sn2ChainLenAbbrevToDbleBondGeomtryMap;
  }

  for $ChainAbbrev (sort keys %{$ChainAbbrevToNameRef}) {
    ($ChainLength, $DoubleBondCount, $DoubleBondGeometry) = ChainAbbrev::ParseChainAbbrev($ChainAbbrev);
    if (exists $ChainLenAbbrevSupportedMapRef->{$ChainLength}) {
      if (exists $ChainLenAbbrevToDbleBondGeomtryMapRef->{$ChainLength}{$DoubleBondCount}) {
	$ChainLenAbbrevToDbleBondGeomtryMapRef->{$ChainLength}{$DoubleBondCount}{$DoubleBondGeometry} = $DoubleBondGeometry;
      }
      else {
	$ChainLenAbbrevToDbleBondGeomtryMapRef->{$ChainLength}{$DoubleBondCount}{$DoubleBondGeometry} = $DoubleBondGeometry;
      }
    }
    else {
      $ChainLenAbbrevSupportedMapRef->{$ChainLength} = $ChainLength;
      $ChainLenAbbrevToDbleBondGeomtryMapRef->{$ChainLength}{$DoubleBondCount}{$DoubleBondGeometry} = $DoubleBondGeometry;
    }
  }
}

# Does chain abbrev specifies a shing-4-enine shingoid base or its analogs? e.g. 17:1(4E)
sub IsSphingosineChainAbbrev {
  my($ChainAbbrev) = @_;
  my($Status, $ChainLength, $DoubleBondCount, $DoubleBondGeometry);

  $Status = 0;
  ($ChainLength, $DoubleBondCount, $DoubleBondGeometry) = ChainAbbrev::ParseChainAbbrev($ChainAbbrev);

  # $Status = ($DoubleBondCount == 1 && $DoubleBondGeometry eq "4E") ? 1 : 0;
  $Status = (($DoubleBondCount == 1 && $DoubleBondGeometry eq "4E") || ($DoubleBondCount == 2 && $DoubleBondGeometry =~ /^4E/)) ? 1 : 0;

  return $Status;
}

# Does chain abbrev specifies a shing-4-enine shingoid base? e.g. 18:1(4E)
sub IsSphingosineC18ChainAbbrev {
  my($ChainAbbrev) = @_;
  my($Status, $ChainLength, $DoubleBondCount, $DoubleBondGeometry);

  $Status = 0;
  ($ChainLength, $DoubleBondCount, $DoubleBondGeometry) = ChainAbbrev::ParseChainAbbrev($ChainAbbrev);

  # $Status = ($ChainLength == 18 && $DoubleBondCount == 1 && $DoubleBondGeometry eq "4E") ? 1 : 0;
  $Status = ($ChainLength == 18 && (($DoubleBondCount == 1 && $DoubleBondGeometry eq "4E") || ($DoubleBondCount == 2 && $DoubleBondGeometry =~ /^4E/)) ) ? 1 : 0;

  return $Status;
}

# Does chain abbrev specifies a shinganine shingoid base? e.g. 18:0
sub IsSphinganineC18ChainAbbrev {
  my($ChainAbbrev) = @_;
  my($Status, $ChainLength, $DoubleBondCount, $DoubleBondGeometry);

  $Status = 0;
  ($ChainLength, $DoubleBondCount, $DoubleBondGeometry) = ChainAbbrev::ParseChainAbbrev($ChainAbbrev);
  $Status = ($ChainLength == 18 && $DoubleBondCount == 0) ? 1 : 0;

  return $Status;
}

1;

__END__

=head1 NAME

SPChainAbbrev - Methods for processing SP chain abbreviations

=head1 SYNOPSIS

use SPChainAbbrev;

use SPChainAbbrev qw(:all);

=head1 DESCRIPTION

SPChainAbbrev module provides these methods:

    ChainAbbrevNameExists - Is it a supported chain abbreviation
    ExpandChainAbbrev - Expand wild cards in chain abbreviation
    GetChainAbbrevToNameMap - Get chain name
    GetChainLenAbbrevSupportedMap - Get reference to supported chain
                                    abbreviations data
    GetChainLenAbbrevDbleBondGeometyDataMap - Get reference to supported chain
                                              double bond geometry data
    GetSupportedChainLenList - Get supported chain lengths
    IsSphingosineChainAbbrev - Is it a sphingosine chain abbreviation
    IsSphingosineC18ChainAbbrev - Is it a sphingosine C18 abbreviation
    IsSphinganineC18ChainAbbrev - Is it a sphinganine C18 abbreviation

=head1 METHODS

=over 4

=item B<ChainAbbrevNameExists>

    $Status = ChainAbbrevNameExists($ChainAbbrev, $ChainType);

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

=item B<GetChainLenAbbrevDbleBondGeometyDataMap>

    $ChainLenDblBondHashRef = GetChainLenAbbrevDbleBondGeometyDataMap();

Return a reference to hash containing information about chain length, number of
double bonds and geometry of double bonds.

=item B<GetSupportedChainLenList>

    $ChainLengthListRef = GetSupportedChainLenList();

Return a reference to a sorted list containing supported chain lengths.

=item B<IsChainAbbrevOkay>

    $Status = IsChainAbbrevOkay($ChainAbbrev);

Return 1 or 0 based on whether chain abbreviation is valid.

=item B<IsSphingosineChainAbbrev>

    $Status = IsSphingosineChainAbbrev($ChainAbbrev);

Return 1 or 0 based on whether it's a sphingosine chain abbreviation.

=item B<IsSphingosineC18ChainAbbrev>

    $Status = IsSphingosineC18ChainAbbrev($ChainAbbrev);

Return 1 or 0 based on whether it's a sphingosine abbreviation with chain length
of 18.

=item B<IsSphinganineC18ChainAbbrev>

    $Status = IsSphinganineC18ChainAbbrev($ChainAbbrev);

Return 1 or 0 based on whether it's a sphinganine abbreviation with chain length
of 18.

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
