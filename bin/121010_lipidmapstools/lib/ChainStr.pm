package ChainStr;
#
# File: ChainStr.pm
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
use LMAPSStr;
use ChainAbbrev;

our ($VERSION, @ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);

$VERSION = '2.00';
@ISA = qw(Exporter);
@EXPORT = qw();
@EXPORT_OK = qw(AssignSubstituentStereoChemistry GenerateAtomBlockLines GenerateBondBlockLines GenerateChainStrData GenerateCmpdCountsLine IsAnySubstituentSpecifiedWithStereoChemistry SetupTemplateDataMap);
%EXPORT_TAGS = (
		all  => [@EXPORT, @EXPORT_OK],
		);

# Assign stereochemistry for substituents...
sub AssignSubstituentStereoChemistry {
  my($CmpdAbbrevTemplateDataMapRef, $CmdDataLinesRef) = @_;
  my($ChainIndex, $SubstituentPosNum, $SubstituentsCount, $SubstituentIndex, $SubstituentSpec, $SubstituentAbbrev, $SubstituentStereoChemistry, $ChiralCenterAtomNum, $SubstituentAtomNum);

  CHAIN: for $ChainIndex (0 .. $#{$CmpdAbbrevTemplateDataMapRef->{SnChainAdd}}) {
    if (!$CmpdAbbrevTemplateDataMapRef->{SnSubstituentsCount}[$ChainIndex]) {
      next CHAIN;
    }
    for $SubstituentPosNum (sort {$a <=> $b} keys %{$CmpdAbbrevTemplateDataMapRef->{SnSubstituentsPosInfo}[$ChainIndex]}) {
      $SubstituentsCount = @{$CmpdAbbrevTemplateDataMapRef->{SnSubstituentsPosInfo}[$ChainIndex]{$SubstituentPosNum}{Spec}};
      SUBSTITUENT: for $SubstituentIndex (0 .. ($SubstituentsCount - 1)) {
	$SubstituentSpec = $CmpdAbbrevTemplateDataMapRef->{SnSubstituentsPosInfo}[$ChainIndex]{$SubstituentPosNum}{Spec}[$SubstituentIndex];
	$SubstituentAbbrev = $CmpdAbbrevTemplateDataMapRef->{SnSubstituentsPosInfo}[$ChainIndex]{$SubstituentPosNum}{Abbrev}[$SubstituentIndex];
	$SubstituentStereoChemistry = $CmpdAbbrevTemplateDataMapRef->{SnSubstituentsPosInfo}[$ChainIndex]{$SubstituentPosNum}{StereoChemistry}[$SubstituentIndex];
	if ($SubstituentStereoChemistry) {
	  $ChiralCenterAtomNum = $CmpdAbbrevTemplateDataMapRef->{SnSubstituentsPosInfo}[$ChainIndex]{$SubstituentPosNum}{ChiralCenterAtomNum}[$SubstituentIndex];
	  $SubstituentAtomNum = $CmpdAbbrevTemplateDataMapRef->{SnSubstituentsPosInfo}[$ChainIndex]{$SubstituentPosNum}{SubstituentAtomNum}[$SubstituentIndex];
	  _AssignStereoChemistry($CmpdAbbrevTemplateDataMapRef, $CmdDataLinesRef, $ChiralCenterAtomNum, $SubstituentAtomNum, $SubstituentStereoChemistry)
	}
      }
    }
  }
}

# Atom block lines...
sub GenerateAtomBlockLines {
    my($CmpdAbbrevTemplateDataMapRef, $Sn1AtomBlockLinesRef, $Sn2AtomBlockLinesRef, $Sn3AtomBlockLinesRef, $Sn4AtomBlockLinesRef, $StrDataLines);

    ($CmpdAbbrevTemplateDataMapRef, $Sn1AtomBlockLinesRef, $Sn2AtomBlockLinesRef, $Sn3AtomBlockLinesRef, $Sn3AtomBlockLinesRef) = (undef) x 5;
    $StrDataLines = '';
    if (@_ == 5) {
      ($CmpdAbbrevTemplateDataMapRef, $Sn1AtomBlockLinesRef, $Sn2AtomBlockLinesRef, $Sn3AtomBlockLinesRef, $Sn4AtomBlockLinesRef) = @_;
    }
    elsif (@_ == 4) {
      ($CmpdAbbrevTemplateDataMapRef, $Sn1AtomBlockLinesRef, $Sn2AtomBlockLinesRef, $Sn3AtomBlockLinesRef) = @_;
    }
    elsif (@_ == 3) {
      ($CmpdAbbrevTemplateDataMapRef, $Sn1AtomBlockLinesRef, $Sn2AtomBlockLinesRef) = @_;
    }
    elsif (@_ == 2) {
      ($CmpdAbbrevTemplateDataMapRef, $Sn1AtomBlockLinesRef) = @_;
    }
    else {
      return $StrDataLines;
    }

    # Template atom block...
    $StrDataLines .= JoinWords(\@{$CmpdAbbrevTemplateDataMapRef->{AtomBlockLines}}, "\n", 0);
    # Sn2 atom block...
    if ($CmpdAbbrevTemplateDataMapRef->{SnChainAdd}[1] && defined($Sn2AtomBlockLinesRef) && @{$Sn2AtomBlockLinesRef} ) { $StrDataLines .= "\n" . JoinWords($Sn2AtomBlockLinesRef, "\n", 0); }
    # Sn1 atom block...
    if ($CmpdAbbrevTemplateDataMapRef->{SnChainAdd}[0] && defined($Sn1AtomBlockLinesRef) && @{$Sn1AtomBlockLinesRef}) { $StrDataLines .= "\n" . JoinWords($Sn1AtomBlockLinesRef, "\n", 0); }
    # Sn3 atom block...
    if ($CmpdAbbrevTemplateDataMapRef->{SnChainAdd}[2] && defined($Sn3AtomBlockLinesRef) && @{$Sn3AtomBlockLinesRef}) { $StrDataLines .= "\n" . JoinWords($Sn3AtomBlockLinesRef, "\n", 0); }
    # Sn4 atom block...
    if ($CmpdAbbrevTemplateDataMapRef->{SnChainAdd}[3] && defined($Sn4AtomBlockLinesRef) && @{$Sn4AtomBlockLinesRef}) { $StrDataLines .= "\n" . JoinWords($Sn4AtomBlockLinesRef, "\n", 0); }

    return $StrDataLines;

}

# Bond block lines...
sub GenerateBondBlockLines {
    my($CmpdAbbrevTemplateDataMapRef, $Sn1BondBlockLinesRef, $Sn2BondBlockLinesRef, $Sn3BondBlockLinesRef, $Sn4BondBlockLinesRef, $StrDataLines);

    ($CmpdAbbrevTemplateDataMapRef, $Sn1BondBlockLinesRef, $Sn2BondBlockLinesRef, $Sn3BondBlockLinesRef, $Sn4BondBlockLinesRef) = (undef) x 5;
    $StrDataLines = '';

    if (@_ == 5) {
      ($CmpdAbbrevTemplateDataMapRef, $Sn1BondBlockLinesRef, $Sn2BondBlockLinesRef, $Sn3BondBlockLinesRef, $Sn4BondBlockLinesRef) = @_;
    }
    elsif (@_ == 4) {
      ($CmpdAbbrevTemplateDataMapRef, $Sn1BondBlockLinesRef, $Sn2BondBlockLinesRef, $Sn3BondBlockLinesRef) = @_;
    }
    elsif (@_ == 3) {
      ($CmpdAbbrevTemplateDataMapRef, $Sn1BondBlockLinesRef, $Sn2BondBlockLinesRef) = @_;
    }
    elsif (@_ == 2) {
      ($CmpdAbbrevTemplateDataMapRef, $Sn1BondBlockLinesRef) = @_;
    }
    else {
      return $StrDataLines;
    }

    # Template bond block..
    $StrDataLines .= JoinWords(\@{$CmpdAbbrevTemplateDataMapRef->{BondBlockLines}}, "\n", 0);
    # Sn2 bond block...
    if ($CmpdAbbrevTemplateDataMapRef->{SnChainAdd}[1] && defined($Sn2BondBlockLinesRef) && @{$Sn2BondBlockLinesRef}) { $StrDataLines .= "\n" . JoinWords($Sn2BondBlockLinesRef, "\n", 0); }
    # Sn1 bond block...
    if ($CmpdAbbrevTemplateDataMapRef->{SnChainAdd}[0] && defined($Sn1BondBlockLinesRef) && @{$Sn1BondBlockLinesRef} ) { $StrDataLines .= "\n" . JoinWords($Sn1BondBlockLinesRef, "\n", 0); }
    # Sn3 bond block...
    if ($CmpdAbbrevTemplateDataMapRef->{SnChainAdd}[2] && defined($Sn3BondBlockLinesRef) && @{$Sn3BondBlockLinesRef}) { $StrDataLines .= "\n". JoinWords($Sn3BondBlockLinesRef, "\n", 0); }
    # Sn4 bond block...
    if ($CmpdAbbrevTemplateDataMapRef->{SnChainAdd}[3] && defined($Sn4BondBlockLinesRef) && @{$Sn4BondBlockLinesRef}) { $StrDataLines .= "\n". JoinWords($Sn4BondBlockLinesRef, "\n", 0); }

    return $StrDataLines;

}

# Generate atom and bond block lines for a chain
sub GenerateChainStrData {
  my($ChainType, $CmpdAbbrevTemplateDataMapRef, $ChainIndex, $SubstituentAtomNum, @AtomBlockLines, @BondBlockLines, @SubstituentAtomBlockLines, @SubstituentBondBlockLines);

  ($ChainType, $CmpdAbbrevTemplateDataMapRef) = @_;

  @AtomBlockLines = ();
  @BondBlockLines = ();

  @SubstituentAtomBlockLines = ();
  @SubstituentBondBlockLines = ();
  $SubstituentAtomNum = 0;

  $ChainIndex = 0;
  SWITCH: {
      if ($ChainType =~ /^Sn1$/i) {$ChainIndex = 0; last SWITCH;}
      if ($ChainType =~ /^Sn2$/i) {$ChainIndex = 1; last SWITCH;}
      if ($ChainType =~ /^Sn3$/i) {$ChainIndex = 2; last SWITCH;}
      if ($ChainType =~ /^Sn4$/i) {$ChainIndex = 3; last SWITCH;}
      return (\@AtomBlockLines, \@BondBlockLines);
  }

  if (!$CmpdAbbrevTemplateDataMapRef->{SnChainAdd}[$ChainIndex]) {
    return (\@AtomBlockLines, \@BondBlockLines);
  }

  if ($CmpdAbbrevTemplateDataMapRef->{SnSubstituentsCount}[$ChainIndex]) {
    $SubstituentAtomNum = _SetupSubstituentStartAtomNum($ChainIndex, $CmpdAbbrevTemplateDataMapRef);
  }

  my($XOffset, $CurSnX, $CurSnY, $ChainAtomNum, $NewSnX, $NewSnY, $AlkyneOffset, $AlkyneYOffset);

  $XOffset = 0.7200;
  $AlkyneYOffset = 0.825;
  $AlkyneOffset = 0;

  ($CurSnX, $CurSnY) = LMAPSStr::ParseCmpdAtomLine($CmpdAbbrevTemplateDataMapRef->{SnAtomLines}[$ChainIndex]);

  # Does chain growth needs to be reversed to go from right to let instead of left to right?
  if (exists($CmpdAbbrevTemplateDataMapRef->{SnReverseChainGrowth}) && $CmpdAbbrevTemplateDataMapRef->{SnReverseChainGrowth}[$ChainIndex]) {
    $XOffset = - $XOffset;
  }

  # Go over atoms in chain...

  for $ChainAtomNum (1 .. ($CmpdAbbrevTemplateDataMapRef->{SnChainLength}[$ChainIndex] - $CmpdAbbrevTemplateDataMapRef->{SnCarbonCount}[$ChainIndex])) {

    $NewSnX = $CurSnX - $XOffset;
    $NewSnX = LMAPSStr::RoundToNextInteger($NewSnX * 1000);
    $NewSnX = 0.0001 + $NewSnX/1000;

    $NewSnY = ($CurSnY > ($CmpdAbbrevTemplateDataMapRef->{SnYMidPoints}[$ChainIndex] + $AlkyneOffset)) ? ($CmpdAbbrevTemplateDataMapRef->{SnY1}[$ChainIndex] + $AlkyneOffset) : ($CmpdAbbrevTemplateDataMapRef->{SnY2}[$ChainIndex] + $AlkyneOffset);

    # Handle end of  double bond...
    #
    # Bond order would be set for (ChainAtomNum -1) to ChainAtomNum later on...
    #
    my($DoubleBondEndAtomNum) = _GetDoubleBondEndAtomNum($ChainIndex, $ChainAtomNum, $CmpdAbbrevTemplateDataMapRef);
    if (exists $CmpdAbbrevTemplateDataMapRef->{SnDbleBondPosToGeometry}[$ChainIndex]{$DoubleBondEndAtomNum}) {
      if ($CmpdAbbrevTemplateDataMapRef->{SnDbleBondPosToGeometry}[$ChainIndex]{$DoubleBondEndAtomNum} eq "Z") {
	$NewSnY = ($CurSnY > ($CmpdAbbrevTemplateDataMapRef->{SnYMidPoints}[$ChainIndex] + $AlkyneOffset)) ? ($CmpdAbbrevTemplateDataMapRef->{SnY2}[$ChainIndex] + $AlkyneOffset) : ($CmpdAbbrevTemplateDataMapRef->{SnY1}[$ChainIndex] + $AlkyneOffset);
      }
      elsif ($CmpdAbbrevTemplateDataMapRef->{SnDbleBondPosToGeometry}[$ChainIndex]{$DoubleBondEndAtomNum} eq "Y") {
	if ($CurSnY > ($CmpdAbbrevTemplateDataMapRef->{SnYMidPoints}[$ChainIndex] + $AlkyneOffset)) {
	  $NewSnY = $CmpdAbbrevTemplateDataMapRef->{SnY2}[$ChainIndex] + ($AlkyneYOffset/2) + $AlkyneOffset;
	  $AlkyneOffset += $AlkyneYOffset;
	}
	else {
	  $NewSnY = $CmpdAbbrevTemplateDataMapRef->{SnY1}[$ChainIndex] - ($AlkyneYOffset/2) + $AlkyneOffset;
	  $AlkyneOffset -= $AlkyneYOffset;
	}
      }
    }

    # Correct for any ring positions...
    my($RingPosNum) = $ChainAtomNum + $CmpdAbbrevTemplateDataMapRef->{SnCarbonCount}[$ChainIndex];
    if ($CmpdAbbrevTemplateDataMapRef->{SnRing}[$ChainIndex] && $CmpdAbbrevTemplateDataMapRef->{SnRingPosInfo}[$ChainIndex]{RingPosToModify}{$RingPosNum}) {
      my($ModifiedSnX, $ModifiedSnY);
      ($ModifiedSnX, $ModifiedSnY) = _ModifyRingChainAtomPosition($ChainIndex, $ChainAtomNum, $CurSnX, $CurSnY, $CmpdAbbrevTemplateDataMapRef);
      $NewSnX = $ModifiedSnX; $NewSnY = $ModifiedSnY;
    }

    # Store current atom line...
    my($SnLine) = LMAPSStr::GenerateCmpdAtomLine($NewSnX, $NewSnY, "0.0000", "C");
    push @AtomBlockLines, $SnLine;

    # Setup first and second atom number for bonds...
    #
    my($PreviousAtomNum, $CurrentAtomNum, $BondOrder, $BondLine);
    ($PreviousAtomNum, $CurrentAtomNum, $BondOrder) = _SetupAtomNumsAndBondOrder($ChainIndex, $ChainAtomNum, $CmpdAbbrevTemplateDataMapRef);

    if ($CmpdAbbrevTemplateDataMapRef->{SnRing}[$ChainIndex] && ($RingPosNum == $CmpdAbbrevTemplateDataMapRef->{SnRingPosInfo}[$ChainIndex]{StartRingPos}) && $CmpdAbbrevTemplateDataMapRef->{SnRingPosInfo}[$ChainIndex]{StartRingStereoChemistry}) {
      my($BondStereo) = ($CmpdAbbrevTemplateDataMapRef->{SnRingPosInfo}[$ChainIndex]{StartRingStereoChemistry} =~ /^alpha$/i) ? 6 : 1;
      $BondLine = LMAPSStr::GenerateCmpdBondLine($CurrentAtomNum, $PreviousAtomNum, $BondOrder, $BondStereo);
    }
    elsif ($CmpdAbbrevTemplateDataMapRef->{SnRing}[$ChainIndex] && ($RingPosNum == ($CmpdAbbrevTemplateDataMapRef->{SnRingPosInfo}[$ChainIndex]{EndRingPos} + 1)) && $CmpdAbbrevTemplateDataMapRef->{SnRingPosInfo}[$ChainIndex]{EndRingStereoChemistry}) {
      my($BondStereo) = ($CmpdAbbrevTemplateDataMapRef->{SnRingPosInfo}[$ChainIndex]{EndRingStereoChemistry} =~ /^beta$/i) ? 1 : 6;
      $BondLine = LMAPSStr::GenerateCmpdBondLine($PreviousAtomNum, $CurrentAtomNum, $BondOrder, $BondStereo);
    }
    else {
      $BondLine = LMAPSStr::GenerateCmpdBondLine($PreviousAtomNum, $CurrentAtomNum, $BondOrder);
    }

    push @BondBlockLines, $BondLine;

    # An additional bond from starting to ending ring atom...
    if ($CmpdAbbrevTemplateDataMapRef->{SnRing}[$ChainIndex] && ($RingPosNum == $CmpdAbbrevTemplateDataMapRef->{SnRingPosInfo}[$ChainIndex]{StartRingPos})) {
      my($StartRingAtomNum, $EndRingAtomNum);
      $StartRingAtomNum = $CurrentAtomNum;
      $EndRingAtomNum = $StartRingAtomNum + $CmpdAbbrevTemplateDataMapRef->{SnRingPosInfo}[$ChainIndex]{RingSize} - 1;

      my($RingBondLine) = LMAPSStr::GenerateCmpdBondLine($StartRingAtomNum, $EndRingAtomNum, 1);
      push @BondBlockLines, $RingBondLine;
    }

    # Handle substituents...
    my($SubstituentPosNum) = $ChainAtomNum + $CmpdAbbrevTemplateDataMapRef->{SnCarbonCount}[$ChainIndex];
    if (exists $CmpdAbbrevTemplateDataMapRef->{SnSubstituentsPosInfo}[$ChainIndex]{$SubstituentPosNum}) {
      # Substituent(s) found at ChainAtomNum: It's two off from substituent position...
      my($SubstituentAtomBlockLinesRef, $SubstituentBondBlockLinesRef);
      ($SubstituentAtomBlockLinesRef, $SubstituentBondBlockLinesRef) = _SetupSubstituentsAtomAndBondBlockLines($ChainIndex, $ChainAtomNum, $SubstituentAtomNum, $CurrentAtomNum, $NewSnX, $NewSnY, $AlkyneOffset, $CmpdAbbrevTemplateDataMapRef);
      $SubstituentAtomNum += $CmpdAbbrevTemplateDataMapRef->{SnSubstituentsPosInfo}[$ChainIndex]{$SubstituentPosNum}{Count};
      push @SubstituentAtomBlockLines, @{$SubstituentAtomBlockLinesRef};
      push @SubstituentBondBlockLines, @{$SubstituentBondBlockLinesRef};
    }

    # Reversing strcture growth for ring in the chain...
    if ($CmpdAbbrevTemplateDataMapRef->{SnRing}[$ChainIndex] && $CmpdAbbrevTemplateDataMapRef->{SnRingPosInfo}[$ChainIndex]{ReverseStructureGrowth}{$RingPosNum}) {
      # Just reverse XOffset...
      $XOffset = - $XOffset;

      # Setup Y1, Y2 and mid point positions...
      my($SnYDiff);
      $SnYDiff = $CmpdAbbrevTemplateDataMapRef->{SnY2}[$ChainIndex] - $CmpdAbbrevTemplateDataMapRef->{SnY1}[$ChainIndex];
      $SnYDiff = abs $SnYDiff;

      # SnY1 is always less than SnY2 to allow for correct generation of Y positions using mid points...
      $CmpdAbbrevTemplateDataMapRef->{SnY1}[$ChainIndex] = $NewSnY;
      $CmpdAbbrevTemplateDataMapRef->{SnY2}[$ChainIndex] = $NewSnY + $SnYDiff;

      $CmpdAbbrevTemplateDataMapRef->{SnYMidPoints}[$ChainIndex] = ($CmpdAbbrevTemplateDataMapRef->{SnY1}[$ChainIndex] + $CmpdAbbrevTemplateDataMapRef->{SnY2}[$ChainIndex])/2;
    }

    # Update chain atom coordinates...
    $CurSnX = $NewSnX;
    $CurSnY = $NewSnY;

  }
  # Add substituent data to atom and bond block lines...
  if (@SubstituentAtomBlockLines) {
    push @AtomBlockLines, @SubstituentAtomBlockLines;
    push @BondBlockLines, @SubstituentBondBlockLines;
  }

  return (\@AtomBlockLines, \@BondBlockLines);
}

# Chain atom loop goes from 1 to chain length minus Sn carbon count and the numbering is offset
# by 2 and the bond order is set for (ChainAtomNum -1) to ChainAtomNum. So, the specified bond
# gemetry position ends up the end of double bond in this scheme...
sub _GetDoubleBondEndAtomNum {
  my($ChainIndex, $ChainAtomNum, $CmpdAbbrevTemplateDataMapRef) = @_;

  my($DoubleBondEndAtomNum) = $ChainAtomNum  + $CmpdAbbrevTemplateDataMapRef->{SnCarbonCount}[$ChainIndex] - 1;

  return $DoubleBondEndAtomNum;
}

# Setup atom numbers and bond order...
sub _SetupAtomNumsAndBondOrder {
  my($ChainIndex, $ChainAtomNum, $CmpdAbbrevTemplateDataMapRef) = @_;

  my($PreviousAtomNum, $CurrentAtomNum, $BondOrder);

  my($Sn1ChainLength, $Sn2ChainLength) = ($CmpdAbbrevTemplateDataMapRef->{SnChainLength}[0], $CmpdAbbrevTemplateDataMapRef->{SnChainLength}[1]);
  my($Sn1AcylCarbons, $Sn2AcylCarbons) = ($CmpdAbbrevTemplateDataMapRef->{SnCarbonCount}[0], $CmpdAbbrevTemplateDataMapRef->{SnCarbonCount}[1]);
  my($Sn1SubstituentCount, $Sn2SubstituentCount) = ($CmpdAbbrevTemplateDataMapRef->{SnSubstituentsCount}[0], $CmpdAbbrevTemplateDataMapRef->{SnSubstituentsCount}[1]);

  if ($ChainIndex == 1) {
    # It's the first block in MolFile data...
    # sn2 acyl chain...
    #
    # Position corresponding to ChainAtomNum - 1
    $PreviousAtomNum = ($ChainAtomNum == 1) ? $CmpdAbbrevTemplateDataMapRef->{SnAtomNums}[$ChainIndex] : ($CmpdAbbrevTemplateDataMapRef->{AtomCount} - 1 + $ChainAtomNum);
    # Position corresponding to ChainAtomNum
    $CurrentAtomNum = $CmpdAbbrevTemplateDataMapRef->{AtomCount} + $ChainAtomNum;
  }
  elsif ($ChainIndex == 0) {
    # sn1 acyl chain..
    if ($Sn2ChainLength >= 2) {
      # Add Sn2 counts as it's present before Sn1 blocks in MolFile data..
      #
      # Position corresponding to ChainAtomNum - 1
      $PreviousAtomNum = ($ChainAtomNum == 1) ? $CmpdAbbrevTemplateDataMapRef->{SnAtomNums}[$ChainIndex] : ($CmpdAbbrevTemplateDataMapRef->{AtomCount} - 1 + $ChainAtomNum + $Sn2ChainLength - $Sn2AcylCarbons + $Sn2SubstituentCount) ;
      # Position corresponding to ChainAtomNum
      $CurrentAtomNum = $CmpdAbbrevTemplateDataMapRef->{AtomCount} + $ChainAtomNum + $Sn2ChainLength - $Sn2AcylCarbons + $Sn2SubstituentCount;
    }
    else {
      # Position corresponding to ChainAtomNum - 1
      $PreviousAtomNum = ($ChainAtomNum == 1) ? $CmpdAbbrevTemplateDataMapRef->{SnAtomNums}[$ChainIndex] : ($CmpdAbbrevTemplateDataMapRef->{AtomCount} - 1 + $ChainAtomNum);
      # Position corresponding to ChainAtomNum
      $CurrentAtomNum = $CmpdAbbrevTemplateDataMapRef->{AtomCount} + $ChainAtomNum;
    }
  }
  elsif ($ChainIndex == 2) {
    # sn3 acyl chain..
    if ($Sn2ChainLength >= 2) {
      # Add Sn2 and Sn1 counts as it's present before Sn3 blocks in MolFile data..
      #
      # Position corresponding to ChainAtomNum - 1
      $PreviousAtomNum = ($ChainAtomNum == 1) ? $CmpdAbbrevTemplateDataMapRef->{SnAtomNums}[$ChainIndex] : ($CmpdAbbrevTemplateDataMapRef->{AtomCount} - 1 + $ChainAtomNum + $Sn1ChainLength + $Sn1SubstituentCount + $Sn2ChainLength + $Sn2SubstituentCount - $Sn1AcylCarbons - $Sn2AcylCarbons);
      # Position corresponding to ChainAtomNum
      $CurrentAtomNum = $CmpdAbbrevTemplateDataMapRef->{AtomCount} + $ChainAtomNum + $Sn1ChainLength + $Sn1SubstituentCount + $Sn2ChainLength + $Sn2SubstituentCount - $Sn1AcylCarbons - $Sn2AcylCarbons;
    }
    else {
      # Add Sn1 counts as it's present before Sn3 blocks in MolFile data..
      # Position corresponding to ChainAtomNum - 1
      $PreviousAtomNum = ($ChainAtomNum == 1) ? $CmpdAbbrevTemplateDataMapRef->{SnAtomNums}[$ChainIndex] : ($CmpdAbbrevTemplateDataMapRef->{AtomCount} - 1 + $ChainAtomNum + $Sn1ChainLength + $Sn1SubstituentCount - $Sn1AcylCarbons);
      # Position corresponding to ChainAtomNum
      $CurrentAtomNum = $CmpdAbbrevTemplateDataMapRef->{AtomCount} + $ChainAtomNum + $Sn1ChainLength + $Sn1SubstituentCount - $Sn1AcylCarbons;
    }
  }
  elsif ($ChainIndex == 3) {
    # sn4 chain during structure generation for Cardiolipins...

    my($Sn3ChainLength) = $CmpdAbbrevTemplateDataMapRef->{SnChainLength}[2];
    my($Sn3AcylCarbons) = $CmpdAbbrevTemplateDataMapRef->{SnCarbonCount}[2];
    my($Sn3SubstituentCount) = $CmpdAbbrevTemplateDataMapRef->{SnSubstituentsCount}[2];

    if ($Sn2ChainLength >= 2 && $Sn3ChainLength >= 2) {
      # Position corresponding to ChainAtomNum - 1
      $PreviousAtomNum = ($ChainAtomNum == 1) ? $CmpdAbbrevTemplateDataMapRef->{SnAtomNums}[$ChainIndex] : ($CmpdAbbrevTemplateDataMapRef->{AtomCount} - 1 + $ChainAtomNum + $Sn1ChainLength + $Sn1SubstituentCount + $Sn2ChainLength + $Sn2SubstituentCount + $Sn3ChainLength + $Sn3SubstituentCount - $Sn1AcylCarbons - $Sn2AcylCarbons - $Sn3AcylCarbons);

      # Position corresponding to ChainAtomNum
      $CurrentAtomNum = $CmpdAbbrevTemplateDataMapRef->{AtomCount} + $ChainAtomNum + $Sn1ChainLength + $Sn1SubstituentCount + $Sn2ChainLength + $Sn2SubstituentCount + $Sn3ChainLength + $Sn3SubstituentCount - $Sn1AcylCarbons - $Sn2AcylCarbons - $Sn3AcylCarbons;
    }
    elsif ($Sn2ChainLength >= 2) {
      # Position corresponding to ChainAtomNum - 1
      $PreviousAtomNum = ($ChainAtomNum == 1) ? $CmpdAbbrevTemplateDataMapRef->{SnAtomNums}[$ChainIndex] : ($CmpdAbbrevTemplateDataMapRef->{AtomCount} - 1 + $ChainAtomNum + $Sn1ChainLength + $Sn1SubstituentCount + $Sn2ChainLength + $Sn2SubstituentCount - $Sn1AcylCarbons - $Sn2AcylCarbons);
      # Position corresponding to ChainAtomNum
      $CurrentAtomNum = $CmpdAbbrevTemplateDataMapRef->{AtomCount} + $ChainAtomNum + $Sn1ChainLength + $Sn1SubstituentCount + $Sn2ChainLength + $Sn2SubstituentCount - $Sn1AcylCarbons - $Sn2AcylCarbons;
    }
    elsif ($Sn3ChainLength >= 2) {
      # Position corresponding to ChainAtomNum - 1
      $PreviousAtomNum = ($ChainAtomNum == 1) ? $CmpdAbbrevTemplateDataMapRef->{SnAtomNums}[$ChainIndex] : ($CmpdAbbrevTemplateDataMapRef->{AtomCount} - 1 + $ChainAtomNum + $Sn1ChainLength + $Sn1SubstituentCount + $Sn3ChainLength + $Sn3SubstituentCount - $Sn1AcylCarbons - $Sn3AcylCarbons);
      # Position corresponding to ChainAtomNum
      $CurrentAtomNum = $CmpdAbbrevTemplateDataMapRef->{AtomCount} + $ChainAtomNum + $Sn1ChainLength + $Sn1SubstituentCount + $Sn3ChainLength + $Sn3SubstituentCount - $Sn1AcylCarbons - $Sn3AcylCarbons;
    }
    else {
      # Position corresponding to ChainAtomNum - 1
      $PreviousAtomNum = ($ChainAtomNum == 1) ? $CmpdAbbrevTemplateDataMapRef->{SnAtomNums}[$ChainIndex] : ($CmpdAbbrevTemplateDataMapRef->{AtomCount} - 1 + $ChainAtomNum + $Sn1ChainLength + $Sn1SubstituentCount - $Sn1AcylCarbons);
      # Position corresponding to ChainAtomNum
      $CurrentAtomNum = $CmpdAbbrevTemplateDataMapRef->{AtomCount} + $ChainAtomNum + $Sn1ChainLength + $Sn1SubstituentCount - $Sn1AcylCarbons;
    }
  }
  # Set the bond order...
  #my($DoubleBondEndAtomNum) = $ChainAtomNum + 1;
  my($DoubleBondEndAtomNum) = _GetDoubleBondEndAtomNum($ChainIndex, $ChainAtomNum, $CmpdAbbrevTemplateDataMapRef);
  $BondOrder = 1;
  if (exists $CmpdAbbrevTemplateDataMapRef->{SnDbleBondPosToGeometry}[$ChainIndex]{$DoubleBondEndAtomNum}) {
    # It's being set for (ChainAtomNum -1) to ChainAtomNum as NextAtom, which is ChainAtomNum + 1, corresponds to
    # double bond position being defined in the specification; in other words, this position takes into
    # account the offset caused by the number of Sn carbons. So, double bond from PreviousAtomNum to CurrentAtomNum
    # corresponds to double bond from (ChainAtomNum - 1) to ChainAtomNum; consquently, ChainAtomNum - 1 + 2 corresponds
    # to actually specified bond position...
    $BondOrder = ($CmpdAbbrevTemplateDataMapRef->{SnDbleBondPosToGeometry}[$ChainIndex]{$DoubleBondEndAtomNum} eq "Y") ? 3 : 2;

  }

  return ($PreviousAtomNum, $CurrentAtomNum, $BondOrder);
}

# Setup atom number for first substituent....
sub _SetupSubstituentStartAtomNum {
  my($ChainIndex, $CmpdAbbrevTemplateDataMapRef) = @_;
  my($SubstituentAtomNum);

  my($Sn1ChainLength, $Sn2ChainLength, $Sn3ChainLength) = ($CmpdAbbrevTemplateDataMapRef->{SnChainLength}[0], $CmpdAbbrevTemplateDataMapRef->{SnChainLength}[1], $CmpdAbbrevTemplateDataMapRef->{SnChainLength}[2]);
  my($Sn1AcylCarbons, $Sn2AcylCarbons, $Sn3AcylCarbons) = ($CmpdAbbrevTemplateDataMapRef->{SnCarbonCount}[0], $CmpdAbbrevTemplateDataMapRef->{SnCarbonCount}[1], $CmpdAbbrevTemplateDataMapRef->{SnCarbonCount}[2]);
  my($Sn1SubstituentCount, $Sn2SubstituentCount, $Sn3SubstituentCount) = ($CmpdAbbrevTemplateDataMapRef->{SnSubstituentsCount}[0], $CmpdAbbrevTemplateDataMapRef->{SnSubstituentsCount}[1], $CmpdAbbrevTemplateDataMapRef->{SnSubstituentsCount}[2],);

  if ($ChainIndex == 1) {
    # It's the first block in MolFile data...
    # sn2 acyl chain...
    $SubstituentAtomNum = $CmpdAbbrevTemplateDataMapRef->{AtomCount} + $Sn2ChainLength - $Sn2AcylCarbons;
  }
  elsif ($ChainIndex == 0) {
    # sn1 acyl chain..
    $SubstituentAtomNum = $CmpdAbbrevTemplateDataMapRef->{AtomCount} + $Sn2ChainLength - $Sn2AcylCarbons + $Sn2SubstituentCount + $Sn1ChainLength - $Sn1AcylCarbons;
  }
  elsif ($ChainIndex == 2) {
    # sn3 acyl chain..
    $SubstituentAtomNum = $CmpdAbbrevTemplateDataMapRef->{AtomCount} + $Sn2ChainLength - $Sn2AcylCarbons + $Sn2SubstituentCount + $Sn1ChainLength - $Sn1AcylCarbons + $Sn1SubstituentCount + $Sn3ChainLength - $Sn3AcylCarbons;
  }
  elsif ($ChainIndex == 3) {
    my($Sn4ChainLength) = $CmpdAbbrevTemplateDataMapRef->{SnChainLength}[3];
    my($Sn4AcylCarbons) = $CmpdAbbrevTemplateDataMapRef->{SnCarbonCount}[3];
    my($Sn4SubstituentCount) = $CmpdAbbrevTemplateDataMapRef->{SnSubstituentsCount}[3];
    # sn4 acyl chain..
    $SubstituentAtomNum = $CmpdAbbrevTemplateDataMapRef->{AtomCount} + $Sn3ChainLength - $Sn3AcylCarbons + $Sn3SubstituentCount  + $Sn2ChainLength - $Sn2AcylCarbons + $Sn2SubstituentCount + $Sn1ChainLength - $Sn1AcylCarbons + $Sn1SubstituentCount + $Sn4ChainLength - $Sn4AcylCarbons;
  }
  return $SubstituentAtomNum;
}

#
# Tasks:
#  . Handle substituents attached to atoms involved with double bonds (Ke)...
#  . Handle substituents attached to atoms involved with muliple bonds (Ep)...
#
sub _SetupSubstituentsAtomAndBondBlockLines {
  my($ChainIndex, $ChainAtomNum, $SubstituentAtomNum, $CurrentAtomNum, $ChainAtomX, $ChainAtomY, $AlkyneOffset, $CmpdAbbrevTemplateDataMapRef) = @_;
  my($SubstituentPosNum, $SubstituentAtomX, $SubstituentAtomY, @SubstituentAtomBlockLines, @SubstituentBondBlockLines);

  @SubstituentAtomBlockLines = (); @SubstituentBondBlockLines = ();

  $SubstituentPosNum = $ChainAtomNum + $CmpdAbbrevTemplateDataMapRef->{SnCarbonCount}[$ChainIndex];
  if (!exists $CmpdAbbrevTemplateDataMapRef->{SnSubstituentsPosInfo}[$ChainIndex]{$SubstituentPosNum}) {
    return (\@SubstituentAtomBlockLines, \@SubstituentBondBlockLines);
  }

  # Is substituent position involved in cis double bonds...
  my($SubstituentPosInvolvedInMultipleBonds, $SubstituentAtStartOfCisDoubleBond, $SubstituentAtEndOfCisDoubleBond, $SubstituentAtEndOfChain);
  $SubstituentPosInvolvedInMultipleBonds = 0; $SubstituentAtStartOfCisDoubleBond = 0; $SubstituentAtEndOfCisDoubleBond = 0;
  $SubstituentAtEndOfChain = 0;

  # Is substituents position involved in a ring...
  my($SubstituentPosInvolvedInRing);
  $SubstituentPosInvolvedInRing = 0;
  if ($CmpdAbbrevTemplateDataMapRef->{SnRing}[$ChainIndex] && ($SubstituentPosNum >= $CmpdAbbrevTemplateDataMapRef->{SnRingPosInfo}[$ChainIndex]{StartRingPos} && $SubstituentPosNum <= $CmpdAbbrevTemplateDataMapRef->{SnRingPosInfo}[$ChainIndex]{EndRingPos})) {
    $SubstituentPosInvolvedInRing = 1;
  }

  my($DoubleBondEndAtomNum) = _GetDoubleBondEndAtomNum($ChainIndex, $ChainAtomNum, $CmpdAbbrevTemplateDataMapRef);
  if (exists $CmpdAbbrevTemplateDataMapRef->{SnDbleBondPosToGeometry}[$ChainIndex]{$SubstituentPosNum}) {
    $SubstituentPosInvolvedInMultipleBonds = 1;
    if ($CmpdAbbrevTemplateDataMapRef->{SnDbleBondPosToGeometry}[$ChainIndex]{$SubstituentPosNum} eq "Z") {
      # Start of a cis double bond...
      $SubstituentAtStartOfCisDoubleBond = 1;
    }
  }
  elsif (exists $CmpdAbbrevTemplateDataMapRef->{SnDbleBondPosToGeometry}[$ChainIndex]{$DoubleBondEndAtomNum}) {
    $SubstituentPosInvolvedInMultipleBonds = 1;
    if ($CmpdAbbrevTemplateDataMapRef->{SnDbleBondPosToGeometry}[$ChainIndex]{$DoubleBondEndAtomNum} eq "Z") {
      # End of a cis double bond...
      $SubstituentAtEndOfCisDoubleBond = 1;
    }
  }

  # Is ChainAtomY above or below mid point?
  my($ChainAtomYAboveMidPoint, $ChainAtomYBelowMidPoint);
  $ChainAtomYAboveMidPoint = ($ChainAtomY > ($CmpdAbbrevTemplateDataMapRef->{SnYMidPoints}[$ChainIndex] + $AlkyneOffset)) ? 1 : 0;
  $ChainAtomYBelowMidPoint = $ChainAtomYAboveMidPoint ? 0 : 1;

  # Any correction for chain proximity...
  my($CorrectForChainProximity) = 0;
  if ($ChainIndex == 0) {
    # Sn1 acyl chain...
    if ($CmpdAbbrevTemplateDataMapRef->{SnChainAdd}[1] && $ChainAtomYBelowMidPoint) {
      # Sn2 acyl chain is present and the substituent is down...
      $CorrectForChainProximity = 1;
    }
  }
  elsif ($ChainIndex == 1) {
    # Sn2 acyl chain...
    if ($CmpdAbbrevTemplateDataMapRef->{SnChainAdd}[0] && $ChainAtomYAboveMidPoint) {
      # Sn1 acyl chain is present and the substituent is up...
      $CorrectForChainProximity = 1;
    }
  }
  # Is substituent at the end of chain?
  my($ChainLength) = $CmpdAbbrevTemplateDataMapRef->{SnChainLength}[$ChainIndex];
  $SubstituentAtEndOfChain = ($ChainLength == $SubstituentPosNum) ? 1 : 0;

  # Go over the substituents...
  my($Index, $Spec, $Abbrev, $StereoChemistry, $XOffset, $YOffset, $XShift, $YShift, $AtomX, $AtomY, $SubStituentsCount);

  $XOffset = 0.72;
  $YOffset = 0.825;

  $SubStituentsCount = @{$CmpdAbbrevTemplateDataMapRef->{SnSubstituentsPosInfo}[$ChainIndex]{$SubstituentPosNum}{Spec}};

  for $Index (0 .. ($SubStituentsCount - 1)) {
    $Spec = $CmpdAbbrevTemplateDataMapRef->{SnSubstituentsPosInfo}[$ChainIndex]{$SubstituentPosNum}{Spec}[$Index];
    $Abbrev = $CmpdAbbrevTemplateDataMapRef->{SnSubstituentsPosInfo}[$ChainIndex]{$SubstituentPosNum}{Abbrev}[$Index];
    $StereoChemistry = $CmpdAbbrevTemplateDataMapRef->{SnSubstituentsPosInfo}[$ChainIndex]{$SubstituentPosNum}{StereoChemistry}[$Index];

    $SubstituentAtomNum++;

    # Store $CurrentAtomNum and $SubstituentAtomNum atom numbers to assign specified chirality for substituents...
    $CmpdAbbrevTemplateDataMapRef->{SnSubstituentsPosInfo}[$ChainIndex]{$SubstituentPosNum}{ChiralCenterAtomNum}[$Index] = $CurrentAtomNum;
    $CmpdAbbrevTemplateDataMapRef->{SnSubstituentsPosInfo}[$ChainIndex]{$SubstituentPosNum}{SubstituentAtomNum}[$Index] = $SubstituentAtomNum;

    $SubstituentAtomX = 0; $SubstituentAtomY = 0;
    $AtomX = 0; $AtomY = 0;

    if ($SubstituentPosInvolvedInRing) {
      ($AtomX, $AtomY) = _SetupRingSubstituentAtomPositions($ChainIndex, $ChainAtomNum, $Index, $ChainAtomX, $ChainAtomY, $CmpdAbbrevTemplateDataMapRef);
    }
    elsif ($Abbrev =~ /^Ep$/i) {
      # Handle substituents on cis double bonds...
      my($EpYOffset) = 0.2;
      if ($SubstituentAtStartOfCisDoubleBond) {
	# Move to left...
	$XShift = -$XOffset/2;
	$YShift = $YOffset - $EpYOffset;
	if ($CorrectForChainProximity) {
	  $YShift = $YShift * 0.7;
	}
      }
      else {
	$XShift = -$XOffset;
	$YShift = $YOffset - $EpYOffset*2;
      }
      $AtomX = $ChainAtomX + $XShift;
      $AtomY = $ChainAtomYAboveMidPoint ? ($ChainAtomY + $YShift) : ($ChainAtomY - $YShift);
    }
    elsif ($Abbrev =~ /^Cp$/i) {
      # Handle substituents on cis double bonds...
      my($CpYOffset) = 0.2;
      if ($SubstituentAtStartOfCisDoubleBond) {
	# Move to left...
	$XShift = -$XOffset/2;
	$YShift = $YOffset - $CpYOffset;
	if ($CorrectForChainProximity) {
	  $YShift = $YShift * 0.7;
	}
      }
      else {
	$XShift = -$XOffset;
	$YShift = $YOffset - $CpYOffset*2;
      }
    
      $AtomX = $ChainAtomX + $XShift;
      $AtomY = $ChainAtomYAboveMidPoint ? ($ChainAtomY + $YShift) : ($ChainAtomY - $YShift);
    }
    elsif ($SubstituentPosInvolvedInMultipleBonds) {
      if ($SubstituentAtStartOfCisDoubleBond || $SubstituentAtEndOfCisDoubleBond) {
	# Handle substituents on cis double bonds...
	if ($SubstituentAtStartOfCisDoubleBond) {
	  # Move to right...
	  $XShift = $XOffset;
	}
	elsif ($SubstituentAtEndOfCisDoubleBond) {
	  # Move to left...
	  $XShift = -$XOffset;
	}
	$YShift = $YOffset;
	if ($CorrectForChainProximity) {
	  $YShift = $YOffset/2;
	}
	
	$AtomX = $ChainAtomX + $XShift;
	$AtomY = $ChainAtomYAboveMidPoint ? ($ChainAtomY + $YShift) : ($ChainAtomY - $YShift);

	# Move X half way close to the middle...
	$AtomX = ($ChainAtomX + $AtomX)/2;
      }
      else {
	# Straight up or down...
	$YShift = $YOffset;
	if ($CorrectForChainProximity) {
	  $YShift = $YOffset/2;
	}
	$AtomX = $ChainAtomX;
	$AtomY = $ChainAtomYAboveMidPoint ? ($ChainAtomY + $YShift) : ($ChainAtomY - $YShift);
      }
    }
    else {
      if ($SubStituentsCount == 3 && $SubstituentAtEndOfChain) {
	# All other positions allow maximum of two sustituents...
	if ($Index == 0) {
	  # First substituent left down or up...
	  $AtomX = $ChainAtomX - $XOffset;
	  $AtomY = ($ChainAtomY > ($CmpdAbbrevTemplateDataMapRef->{SnYMidPoints}[$ChainIndex] + $AlkyneOffset)) ? ($CmpdAbbrevTemplateDataMapRef->{SnY1}[$ChainIndex] + $AlkyneOffset) : ($CmpdAbbrevTemplateDataMapRef->{SnY2}[$ChainIndex] + $AlkyneOffset);
	}
	else {
	  if ($Index == 1) {
	    # First substituent on the left...
	    $XShift =  -$XOffset;
	    $YShift = $YOffset;
	  }
	  else {
	    # Second substituent on the right...
	    $XShift =  $XOffset;
	    $YShift = $YOffset;
	  }
	  if ($CorrectForChainProximity) {
	    $YShift = $YOffset/2;
	  }
	  $AtomX = $ChainAtomX + $XShift;
	  $AtomY = $ChainAtomYAboveMidPoint ? ($ChainAtomY + $YShift) : ($ChainAtomY - $YShift);
	  # Move X half way close to the middle...
	  $AtomX = ($ChainAtomX + $AtomX)/2;
	}
      }
      elsif ($SubStituentsCount == 2) {
	if ($SubstituentAtEndOfChain) {
	  if ($Index == 0) {
	    # First substituent left down or up...
	    $AtomX = $ChainAtomX - $XOffset;
	    $AtomY = ($ChainAtomY > ($CmpdAbbrevTemplateDataMapRef->{SnYMidPoints}[$ChainIndex] + $AlkyneOffset)) ? ($CmpdAbbrevTemplateDataMapRef->{SnY1}[$ChainIndex] + $AlkyneOffset) : ($CmpdAbbrevTemplateDataMapRef->{SnY2}[$ChainIndex] + $AlkyneOffset);
	  }
	  else {
	    # Second substituent straight up or down...
	    $YShift = $YOffset;
	    if ($CorrectForChainProximity) {
	      $YShift = $YOffset/2;
	    }
	    $AtomX = $ChainAtomX;
	    $AtomY = $ChainAtomYAboveMidPoint ? ($ChainAtomY + $YShift) : ($ChainAtomY - $YShift);
	  }
	}
	else {
	  if ($Index == 0) {
	    # First substituent on the left...
	    $XShift =  -$XOffset;
	    $YShift = $YOffset;
	  }
	  else {
	    # Second substituent on the right...
	    $XShift =  $XOffset;
	    $YShift = $YOffset;
	  }
	  if ($CorrectForChainProximity) {
	    $YShift = $YOffset/2;
	  }
	  $AtomX = $ChainAtomX + $XShift;
	  $AtomY = $ChainAtomYAboveMidPoint ? ($ChainAtomY + $YShift) : ($ChainAtomY - $YShift);
	  # Move X half way close to the middle...
	  $AtomX = ($ChainAtomX + $AtomX)/2;
	}
      }
      else {
	if ($SubstituentAtEndOfChain) {
	  # Left down or up...
	  $AtomX = $ChainAtomX - $XOffset;
	  $AtomY = ($ChainAtomY > ($CmpdAbbrevTemplateDataMapRef->{SnYMidPoints}[$ChainIndex] + $AlkyneOffset)) ? ($CmpdAbbrevTemplateDataMapRef->{SnY1}[$ChainIndex] + $AlkyneOffset) : ($CmpdAbbrevTemplateDataMapRef->{SnY2}[$ChainIndex] + $AlkyneOffset);
	}
	else {
	  # Straight up or down...
	  $YShift = $YOffset;
	  if ($CorrectForChainProximity) {
	    $YShift = $YOffset/2;
	  }
	  $AtomX = $ChainAtomX;
	  $AtomY = $ChainAtomYAboveMidPoint ? ($ChainAtomY + $YShift) : ($ChainAtomY - $YShift);
	}
      }
    }

    # Setup substituent coordinates...
    $SubstituentAtomX = $AtomX;
    $SubstituentAtomY = $AtomY;

    # Setup atom data line...
    my($AtomSpec) = $Abbrev;
    if ($Abbrev =~ /^Ke$/i) {
      $AtomSpec = 'O';
    }
    elsif ($Abbrev =~ /^My$/i) {
      $AtomSpec = 'C';
    }
    elsif ($Abbrev =~ /^Me$/i) {
      $AtomSpec = 'C';
    }
    elsif ($Abbrev =~ /^Ep$/i) {
      $AtomSpec = 'O';
    }
    elsif ($Abbrev =~ /^Cp$/i) {
      $AtomSpec = 'C';
    }
    elsif ($Abbrev =~ /^OH$/i) {
      $AtomSpec = 'O';
    }
    elsif ($Abbrev =~ /^NH2$/i) {
      $AtomSpec = 'N';
    }
    elsif ($Abbrev =~ /^SH$/i) {
      $AtomSpec = 'S';
    }

    my($AtomLine) = LMAPSStr::GenerateCmpdAtomLine($SubstituentAtomX, $SubstituentAtomY, "0.0000", $AtomSpec);
    push @SubstituentAtomBlockLines, $AtomLine;

    # Setup bond line...
    my($BondOrder) = ChainAbbrev::GetSubstituentBondOrder($Abbrev);
    my($BondLine) = LMAPSStr::GenerateCmpdBondLine($CurrentAtomNum, $SubstituentAtomNum, $BondOrder);
    push @SubstituentBondBlockLines, $BondLine;

    # Add another bond line for epoxy with next atom...
    if ($Abbrev =~ /^[CE]p$/i) {
      my($SecondBondLine) = LMAPSStr::GenerateCmpdBondLine(($CurrentAtomNum + 1), $SubstituentAtomNum, $BondOrder);
      push @SubstituentBondBlockLines, $SecondBondLine;
    }
  }
  return (\@SubstituentAtomBlockLines, \@SubstituentBondBlockLines);
}

# Setup positions of substituents attached to a ring...
sub _SetupRingSubstituentAtomPositions {
  my($ChainIndex, $ChainAtomNum, $SubstituentIndex, $ChainAtomX, $ChainAtomY, $CmpdAbbrevTemplateDataMapRef) = @_;
  my($SubstituentPosNum,  $SubstituentSpec, $SubstituentAbbrev, $SubStituentsCount,  $StartRingPos, $EndRingPos, $RingSize, $SubstituentAtomX, $SubstituentAtomY, $XOffSet, $YOffSet, $XShift, $YShift);

  $XOffSet = 0.72;
  $YOffSet = 0.825;

  $SubstituentPosNum = $ChainAtomNum + $CmpdAbbrevTemplateDataMapRef->{SnCarbonCount}[$ChainIndex];
  $SubstituentSpec = $CmpdAbbrevTemplateDataMapRef->{SnSubstituentsPosInfo}[$ChainIndex]{$SubstituentPosNum}{Spec}[$SubstituentIndex];
  $SubstituentAbbrev = $CmpdAbbrevTemplateDataMapRef->{SnSubstituentsPosInfo}[$ChainIndex]{$SubstituentPosNum}{Abbrev}[$SubstituentIndex];

  $SubStituentsCount = @{$CmpdAbbrevTemplateDataMapRef->{SnSubstituentsPosInfo}[$ChainIndex]{$SubstituentPosNum}{Spec}};

  $StartRingPos = $CmpdAbbrevTemplateDataMapRef->{SnRingPosInfo}[$ChainIndex]{StartRingPos};
  $EndRingPos = $CmpdAbbrevTemplateDataMapRef->{SnRingPosInfo}[$ChainIndex]{EndRingPos};
  $RingSize = $CmpdAbbrevTemplateDataMapRef->{SnRingPosInfo}[$ChainIndex]{RingSize};

  if ($SubstituentAbbrev =~ /^[CE]p$/i) {
    if ($RingSize == 5) {
      if ($SubstituentPosNum == $StartRingPos) {
	$SubstituentAtomX = $ChainAtomX - 0.2190;
	$SubstituentAtomY = $ChainAtomY + 0.8591;
      }
      elsif ($SubstituentPosNum == ($StartRingPos + 1)) {
	$SubstituentAtomX = $ChainAtomX - 0.9012;
	$SubstituentAtomY = $ChainAtomY - 0.0510;
      }
      elsif ($SubstituentPosNum == ($StartRingPos + 2)) {
	$SubstituentAtomX = $ChainAtomX - 0.4725;
	$SubstituentAtomY = $ChainAtomY - 0.7145;
      }
      elsif ($SubstituentPosNum == ($StartRingPos + 3)) {
	$SubstituentAtomX = $ChainAtomX + 0.5721;
	$SubstituentAtomY = $ChainAtomY - 0.5435;
      }
    }
    elsif ($RingSize == 6) {
      if ($SubstituentPosNum == $StartRingPos) {
	$SubstituentAtomX = $ChainAtomX;
	$SubstituentAtomY = $ChainAtomY + 0.8225;
      }
      elsif ($SubstituentPosNum == ($StartRingPos + 1)) {
	$SubstituentAtomX = $ChainAtomX - 0.7159;
	$SubstituentAtomY = $ChainAtomY + 0.4150;
      }
      elsif ($SubstituentPosNum == ($StartRingPos + 2)) {
	$SubstituentAtomX = $ChainAtomX - 0.7144;
	$SubstituentAtomY = $ChainAtomY - 0.4145;
      }
      elsif ($SubstituentPosNum == ($StartRingPos + 3)) {
	$SubstituentAtomX = $ChainAtomX;
	$SubstituentAtomY = $ChainAtomY - 0.8226;
      }
      elsif ($SubstituentPosNum == ($StartRingPos + 4)) {
	$SubstituentAtomX = $ChainAtomX + 0.7145;
	$SubstituentAtomY = $ChainAtomY - 0.4125;
      }
    }
  }
  else {
    if ($RingSize == 5) {
      if ($SubstituentPosNum == $StartRingPos) {
	$SubstituentAtomX = $ChainAtomX;
	$SubstituentAtomY = $ChainAtomY + $YOffSet;
      }
      elsif ($SubstituentPosNum == ($StartRingPos + 1)) {
	if ($SubStituentsCount == 1) {
	  $SubstituentAtomX = $ChainAtomX - 0.2135;
	  $SubstituentAtomY = $ChainAtomY + 0.7949;
	}
	else {
	  if ($SubstituentIndex == 0) {
	    # First substituent: Top right...
	    $SubstituentAtomX = $ChainAtomX + 0.2136;
	    $SubstituentAtomY = $ChainAtomY + 0.7949;
	  }
	  else {
	    # Second substituent: Top left...
	    $SubstituentAtomX = $ChainAtomX - 0.5833;
	    $SubstituentAtomY = $ChainAtomY + 0.5834;
	  }
	}
      }
      elsif ($SubstituentPosNum == ($StartRingPos + 2)) {
	if ($SubStituentsCount == 1) {
	  $SubstituentAtomX = $ChainAtomX - $XOffSet;
	  $SubstituentAtomY = $ChainAtomY;
	}
	else {
	  if ($SubstituentIndex == 0) {
	    # First substituent: Left up
	    $XShift = -$XOffSet;
	    $YShift = $YOffSet;
	  }
	  else {
	    # Second substituent: Left down
	    $XShift = -$XOffSet;
	    $YShift = -$YOffSet;
	  }
	  $SubstituentAtomX = $ChainAtomX + $XShift;
	  $SubstituentAtomY = $ChainAtomY + $YShift;
	  # Move Y half way to the middle...
	  $SubstituentAtomY = ($SubstituentAtomY + $ChainAtomY)/2;
	}
      }
      elsif ($SubstituentPosNum == ($StartRingPos + 3)) {
	if ($SubStituentsCount == 1) {
	  $SubstituentAtomX = $ChainAtomX - 0.2135;
	  $SubstituentAtomY = $ChainAtomY - 0.7968;
	}
	else {
	  if ($SubstituentIndex == 0) {
	    # First substituent: Bottom right...
	    $SubstituentAtomX = $ChainAtomX + 0.2136;
	    $SubstituentAtomY = $ChainAtomY - 0.7968;
	  }
	  else {
	    # Second substituent: Bottom left...
	    $SubstituentAtomX = $ChainAtomX - 0.7144;
	    $SubstituentAtomY = $ChainAtomY - 0.4145;
	  }
	}
      }
      elsif ($SubstituentPosNum == ($StartRingPos + 4)) {
	$SubstituentAtomX = $ChainAtomX;
	$SubstituentAtomY = $ChainAtomY - $YOffSet;
      }
    }
    elsif ($RingSize == 6) {
      if ($SubstituentPosNum == $StartRingPos) {
	$SubstituentAtomX = $ChainAtomX;
	$SubstituentAtomY = $ChainAtomY + $YOffSet;
      }
      elsif ($SubstituentPosNum == ($StartRingPos + 1)) {
	if ($SubStituentsCount == 1) {
	  $SubstituentAtomX = $ChainAtomX;
	  $SubstituentAtomY = $ChainAtomY + $YOffSet;
	}
	else {
	  if ($SubstituentIndex == 0) {
	    # First substituent: Top right...
	    $XShift = $XOffSet;
	    $YShift = $YOffSet;
	  }
	  else {
	    # Second substituent: Top left...
	    $XShift = -$XOffSet;
	    $YShift = $YOffSet;
	  }
	  $SubstituentAtomX = $ChainAtomX + $XShift;
	  $SubstituentAtomY = $ChainAtomY + $YShift;
	  # Move X half way to the middle...
	  $SubstituentAtomX = ($SubstituentAtomX + $ChainAtomX)/2;
	}
      }
      elsif ($SubstituentPosNum == ($StartRingPos + 2)) {
	if ($SubStituentsCount == 1) {
	  $SubstituentAtomX = $ChainAtomX - 0.7145;
	  $SubstituentAtomY = $ChainAtomY + 0.4125;
	}
	else {
	  if ($SubstituentIndex == 0) {
	    # First substituent: Left top...
	    $SubstituentAtomX = $ChainAtomX - 0.4125;
	    $SubstituentAtomY = $ChainAtomY + 0.7144;
	  }
	  else {
	    # Second substituent: Left down
	    $SubstituentAtomX = $ChainAtomX - 0.8250;
	    $SubstituentAtomY = $ChainAtomY;
	  }
	}
      }
      elsif ($SubstituentPosNum == ($StartRingPos + 3)) {
	if ($SubStituentsCount == 1) {
	  $SubstituentAtomX = $ChainAtomX - 0.7145;
	  $SubstituentAtomY = $ChainAtomY - 0.4125;
	}
	else {
	  if ($SubstituentIndex == 0) {
	    # First substituent: Left down...
	    $SubstituentAtomX = $ChainAtomX - 0.4125;
	    $SubstituentAtomY = $ChainAtomY - 0.7144;
	  }
	  else {
	    # Second substituent: Left top...
	    $SubstituentAtomX = $ChainAtomX - 0.8250;
	    $SubstituentAtomY = $ChainAtomY;
	  }
	}
      }
      elsif ($SubstituentPosNum == ($StartRingPos + 4)) {
	if ($SubStituentsCount == 1) {
	  $SubstituentAtomX = $ChainAtomX;
	  $SubstituentAtomY = $ChainAtomY - $YOffSet;
	}
	else {
	  if ($SubstituentIndex == 0) {
	    # First substituent: Bottom right...
	    $XShift = $XOffSet;
	    $YShift = -$YOffSet;
	  }
	  else {
	    # Second substituent: Bottom left...
	    $XShift = -$XOffSet;
	    $YShift = -$YOffSet;
	  }
	  $SubstituentAtomX = $ChainAtomX + $XShift;
	  $SubstituentAtomY = $ChainAtomY + $YShift;
	  # Move X half way to the middle...
	  $SubstituentAtomX = ($SubstituentAtomX + $ChainAtomX)/2;
	}
      }
      elsif ($SubstituentPosNum == ($StartRingPos + 5)) {
	$SubstituentAtomX = $ChainAtomX;
	$SubstituentAtomY = $ChainAtomY - $YOffSet;
      }
    }
  }

  return ($SubstituentAtomX, $SubstituentAtomY);
}

# Modify ring atom positions for five and six membered rings. For five membered ring,
# four atom positions are modified: start + 1 to end. And for six membvered ring, only
# one position need to be changed: start + 3
#
sub _ModifyRingChainAtomPosition {
  my($ChainIndex, $ChainAtomNum, $PreviousSnX, $PreviousSnY, $CmpdAbbrevTemplateDataMapRef) = @_;
  my($ModifiedSnX, $ModifiedSnY, $RingSize, $RingPos, $StartRingPos, $EndRingPos, $PosXOffSet, $PosYOffSet);

  ($ModifiedSnX, $ModifiedSnY) = (0, 0);

  $RingPos = $ChainAtomNum + $CmpdAbbrevTemplateDataMapRef->{SnCarbonCount}[$ChainIndex];

  $StartRingPos = $CmpdAbbrevTemplateDataMapRef->{SnRingPosInfo}[$ChainIndex]{StartRingPos};
  $EndRingPos = $CmpdAbbrevTemplateDataMapRef->{SnRingPosInfo}[$ChainIndex]{EndRingPos};
  $RingSize = $CmpdAbbrevTemplateDataMapRef->{SnRingPosInfo}[$ChainIndex]{RingSize};

  $PosXOffSet = $CmpdAbbrevTemplateDataMapRef->{SnRingPosInfo}[$ChainIndex]{PosXOffSet}{$RingPos};
  $PosYOffSet = $CmpdAbbrevTemplateDataMapRef->{SnRingPosInfo}[$ChainIndex]{PosYOffSet}{$RingPos};

  if (exists($CmpdAbbrevTemplateDataMapRef->{SnReverseChainGrowth}) && $CmpdAbbrevTemplateDataMapRef->{SnReverseChainGrowth}[$ChainIndex]) {
    $PosXOffSet = - $PosXOffSet;
  }

  if ($RingSize == 5) {
    if ($RingPos == $StartRingPos) {
      $ModifiedSnX = $PreviousSnX  - $PosXOffSet;
      $ModifiedSnY = $PreviousSnY - $PosYOffSet;
    }
    elsif ($RingPos == ($StartRingPos + 1)) {
      $ModifiedSnX = $PreviousSnX  - $PosXOffSet;
      $ModifiedSnY = $PreviousSnY + $PosYOffSet;
    }
    elsif ($RingPos == ($StartRingPos + 2)) {
      $ModifiedSnX = $PreviousSnX  - $PosXOffSet;
      $ModifiedSnY = $PreviousSnY - $PosYOffSet;
    }
    elsif ($RingPos == ($StartRingPos + 3)) {
      $ModifiedSnX = $PreviousSnX  + $PosXOffSet;
      $ModifiedSnY = $PreviousSnY - $PosYOffSet;
    }
    elsif ($RingPos == ($StartRingPos + 4)) {
      $ModifiedSnX = $PreviousSnX  + $PosXOffSet;
      $ModifiedSnY = $PreviousSnY + $PosYOffSet;
    }
    elsif ($RingPos == ($StartRingPos + 5)) {
      # Position after last ring position...
      $ModifiedSnX = $PreviousSnX  + $PosXOffSet;
      $ModifiedSnY = $PreviousSnY - $PosYOffSet;
    }
  }
  elsif ($RingSize == 6) {
    # Only one position to modify...
    if ($RingPos == $StartRingPos) {
      $ModifiedSnX = $PreviousSnX - $PosXOffSet;
      $ModifiedSnY = $PreviousSnY - $PosYOffSet;
    }
    elsif ($RingPos == ($StartRingPos + 1)) {
      $ModifiedSnX = $PreviousSnX - $PosXOffSet;
      $ModifiedSnY = $PreviousSnY + $PosYOffSet;
    }
    elsif ($RingPos == ($StartRingPos + 2)) {
      $ModifiedSnX = $PreviousSnX - $PosXOffSet;
      $ModifiedSnY = $PreviousSnY - $PosYOffSet;
    }
    elsif ($RingPos == ($StartRingPos + 3)) {
      $ModifiedSnX = $PreviousSnX;
      $ModifiedSnY = $PreviousSnY - $PosYOffSet;
    }
    elsif ($RingPos == ($StartRingPos + 4)) {
      $ModifiedSnX = $PreviousSnX + $PosXOffSet;
      $ModifiedSnY = $PreviousSnY - $PosYOffSet;
    }
    elsif ($RingPos == ($StartRingPos + 5)) {
      $ModifiedSnX = $PreviousSnX + $PosXOffSet;
      $ModifiedSnY = $PreviousSnY + $PosYOffSet;
    }
    elsif ($RingPos == ($StartRingPos + 6)) {
      $ModifiedSnX = $PreviousSnX + $PosXOffSet;
      $ModifiedSnY = $PreviousSnY - $PosYOffSet;
    }
  }
  return ($ModifiedSnX, $ModifiedSnY);
}

sub GenerateCmpdCountsLine {
    my($CmpdAbbrevTemplateDataMapRef, $CmpdAtomCount, $ChainIndex, $CmpdBondCount);

    ($CmpdAbbrevTemplateDataMapRef) = @_;

    # Chain atom counts...
    $CmpdAtomCount = $CmpdAbbrevTemplateDataMapRef->{AtomCount};
    $CmpdBondCount = $CmpdAbbrevTemplateDataMapRef->{BondCount};
    CHAIN: for $ChainIndex (0 .. $#{$CmpdAbbrevTemplateDataMapRef->{SnChainAdd}}) {
      if (!$CmpdAbbrevTemplateDataMapRef->{SnChainAdd}[$ChainIndex]) {
	next CHAIN;
      }
      $CmpdAtomCount += $CmpdAbbrevTemplateDataMapRef->{SnChainLength}[$ChainIndex];
      $CmpdBondCount += $CmpdAbbrevTemplateDataMapRef->{SnChainLength}[$ChainIndex];
    }
    $CmpdAtomCount -= $CmpdAbbrevTemplateDataMapRef->{TotalSnCarbons};
    $CmpdBondCount -= $CmpdAbbrevTemplateDataMapRef->{TotalSnCarbons};

    #$CmpdAtomCount = $CmpdAbbrevTemplateDataMapRef->{AtomCount} + $CmpdAbbrevTemplateDataMapRef->{SnChainLength}[1] + $CmpdAbbrevTemplateDataMapRef->{SnChainLength}[0] + $CmpdAbbrevTemplateDataMapRef->{SnChainLength}[2] - $CmpdAbbrevTemplateDataMapRef->{TotalSnCarbons};
    #$CmpdBondCount = $CmpdAbbrevTemplateDataMapRef->{BondCount} + $CmpdAbbrevTemplateDataMapRef->{SnChainLength}[1] + $CmpdAbbrevTemplateDataMapRef->{SnChainLength}[0] + $CmpdAbbrevTemplateDataMapRef->{SnChainLength}[2] - $CmpdAbbrevTemplateDataMapRef->{TotalSnCarbons};

    # Substituents...
    my($SubstituentAtomsCount, $SubstituentBondsCount);
    $SubstituentAtomsCount = 0; $SubstituentBondsCount = 0;
    CHAIN: for $ChainIndex (0 .. $#{$CmpdAbbrevTemplateDataMapRef->{SnChainAdd}}) {
      if (!$CmpdAbbrevTemplateDataMapRef->{SnChainAdd}[$ChainIndex]) {
	next CHAIN;
      }
      if ($CmpdAbbrevTemplateDataMapRef->{SnSubstituentsCount}[$ChainIndex]) {
	$SubstituentAtomsCount += $CmpdAbbrevTemplateDataMapRef->{SnSubstituentsCount}[$ChainIndex];
	$SubstituentBondsCount += $CmpdAbbrevTemplateDataMapRef->{SnSubstituentsBondCount}[$ChainIndex];
      }
    }
    $CmpdAtomCount +=  $SubstituentAtomsCount;
    $CmpdBondCount += $SubstituentBondsCount;

    # Rings...
    my($RingBondCount);
    $RingBondCount = 0;
    CHAIN: for $ChainIndex (0 .. $#{$CmpdAbbrevTemplateDataMapRef->{SnChainAdd}}) {
      if (!$CmpdAbbrevTemplateDataMapRef->{SnChainAdd}[$ChainIndex]) {
	next CHAIN;
      }
      if ($CmpdAbbrevTemplateDataMapRef->{SnRing}[$ChainIndex]) {
	$RingBondCount += $CmpdAbbrevTemplateDataMapRef->{SnRing}[$ChainIndex];
      }
    }
    $CmpdBondCount += $RingBondCount;

    return LMAPSStr::GenerateCmpdCountsLine($CmpdAtomCount, $CmpdBondCount);
}

# Is any substituent specified along with its stereochemistry...
sub IsAnySubstituentSpecifiedWithStereoChemistry {
  my($ChainIndex, $SubstituentPosNum, $SubstituentsCount, $SubstituentIndex, $SubstituentSpec, $SubstituentAbbrev, $SubstituentStereoChemistry, $CmpdAbbrevTemplateDataMapRef);

  ($CmpdAbbrevTemplateDataMapRef) = @_;

  # Any substituent specified...
  CHAIN: for $ChainIndex (0 .. $#{$CmpdAbbrevTemplateDataMapRef->{SnChainAdd}}) {
    if (!$CmpdAbbrevTemplateDataMapRef->{SnSubstituentsCount}[$ChainIndex]) {
      next CHAIN;
    }
    for $SubstituentPosNum (keys %{$CmpdAbbrevTemplateDataMapRef->{SnSubstituentsPosInfo}[$ChainIndex]}) {
      $SubstituentsCount = @{$CmpdAbbrevTemplateDataMapRef->{SnSubstituentsPosInfo}[$ChainIndex]{$SubstituentPosNum}{Spec}};
      for $SubstituentIndex (0 .. ($SubstituentsCount - 1)) {
	$SubstituentSpec = $CmpdAbbrevTemplateDataMapRef->{SnSubstituentsPosInfo}[$ChainIndex]{$SubstituentPosNum}{Spec}[$SubstituentIndex];
	$SubstituentAbbrev = $CmpdAbbrevTemplateDataMapRef->{SnSubstituentsPosInfo}[$ChainIndex]{$SubstituentPosNum}{Abbrev}[$SubstituentIndex];
	$SubstituentStereoChemistry = $CmpdAbbrevTemplateDataMapRef->{SnSubstituentsPosInfo}[$ChainIndex]{$SubstituentPosNum}{StereoChemistry}[$SubstituentIndex];
	if ($SubstituentStereoChemistry) {
	  return 1;
	}
      }
    }
  }
  return 0;
}

# Set up template data...
sub SetupTemplateDataMap {
  my($TemplateType, $AbbrevTemplateDataMapRef, $TemplateData) = @_;

  my($AbbrevID, $HeadGroup, $HeadGroupNameBeforeBase, $HeadGroupAbbrev,$Y1Sn1, $Y2Sn1, $Y1Sn2, $Y2Sn2, $Y1Sn3, $Y2Sn3, $Y1Sn4, $Y2Sn4, $Sn1AtomNum, $Sn2AtomNum, $Sn3AtomNum, $Sn4AtomNum, $Sn1CarbonCount, $Sn2CarbonCount, $Sn3CarbonCount, $Sn4CarbonCount, $Sn2CAtomNum, $Sn2OAtomNum, $Sn2HAtomNum, $LMCategory, $LMMainClass, $LMSubClass, $CmpdString, $Sn1ReverseChainGrowth, $Sn2ReverseChainGrowth, $Sn3ReverseChainGrowth, $Sn4ReverseChainGrowth);

  ($Sn1ReverseChainGrowth, $Sn2ReverseChainGrowth, $Sn3ReverseChainGrowth, $Sn4ReverseChainGrowth) = (undef) x 4;

  if ($TemplateType =~ /^FA$/) {
    ($AbbrevID, $Y1Sn1, $Y2Sn1, $Y1Sn2, $Y2Sn2, $Y1Sn3, $Y2Sn3, $Sn1AtomNum, $Sn2AtomNum, $Sn3AtomNum, $Sn1CarbonCount, $Sn2CarbonCount, $Sn3CarbonCount, $Sn2CAtomNum, $Sn2OAtomNum, $Sn2HAtomNum, $LMCategory, $LMMainClass, $LMSubClass, $CmpdString) = split /\|/, $TemplateData;
  }
  elsif ($TemplateType =~ /^GL$/) {
    ($AbbrevID, $Y1Sn1, $Y2Sn1, $Y1Sn2, $Y2Sn2, $Y1Sn3, $Y2Sn3, $Sn1AtomNum, $Sn2AtomNum, $Sn3AtomNum, $Sn1CarbonCount, $Sn2CarbonCount, $Sn3CarbonCount, $Sn2CAtomNum, $Sn2OAtomNum, $Sn2HAtomNum, $LMCategory, $LMMainClass, $LMSubClass, $CmpdString) = split /\|/, $TemplateData;
  }
  elsif ($TemplateType =~ /^GP$/) {
    ($AbbrevID, $HeadGroup, $Y1Sn1, $Y2Sn1, $Y1Sn2, $Y2Sn2, $Y1Sn3, $Y2Sn3, $Sn1AtomNum, $Sn2AtomNum, $Sn3AtomNum, $Sn1CarbonCount, $Sn2CarbonCount, $Sn3CarbonCount, $Sn2CAtomNum, $Sn2OAtomNum, $Sn2HAtomNum, $LMCategory, $LMMainClass, $LMSubClass, $CmpdString) = split /\|/, $TemplateData;
    $AbbrevTemplateDataMapRef->{HeadGroupName} = $HeadGroup;
  }
  elsif ($TemplateType =~ /^CL$/) {
    ($AbbrevID, $HeadGroup, $Y1Sn1, $Y2Sn1, $Y1Sn2, $Y2Sn2, $Y1Sn3, $Y2Sn3, $Y1Sn4, $Y2Sn4, $Sn1AtomNum, $Sn2AtomNum, $Sn3AtomNum, $Sn4AtomNum, $Sn1CarbonCount, $Sn2CarbonCount, $Sn3CarbonCount, $Sn4CarbonCount, $Sn2CAtomNum, $Sn2OAtomNum, $Sn2HAtomNum, $LMCategory, $LMMainClass, $LMSubClass, $CmpdString) = split /\|/, $TemplateData;
    $AbbrevTemplateDataMapRef->{HeadGroupName} = $HeadGroup;
  }
  elsif ($TemplateType =~ /^SP$/) {
    ($AbbrevID, $HeadGroup, $HeadGroupNameBeforeBase, $HeadGroupAbbrev, $Y1Sn1, $Y2Sn1, $Y1Sn2, $Y2Sn2, $Y1Sn3, $Y2Sn3, $Sn1AtomNum, $Sn2AtomNum, $Sn3AtomNum, $Sn1CarbonCount, $Sn2CarbonCount, $Sn3CarbonCount, $Sn2CAtomNum, $Sn2OAtomNum, $Sn2HAtomNum, $LMCategory, $LMMainClass, $LMSubClass, $CmpdString) = split /\|/, $TemplateData;
    $AbbrevTemplateDataMapRef->{HeadGroupName} = $HeadGroup;
    $AbbrevTemplateDataMapRef->{HeadGroupNameBeforeBase} = $HeadGroupNameBeforeBase;
    $AbbrevTemplateDataMapRef->{HeadGroupAbbrev} = (defined($HeadGroupAbbrev) && $HeadGroupAbbrev) ? $HeadGroupAbbrev: '';
  }
  elsif ($TemplateType =~ /^LM$/) {
    ($AbbrevID, $HeadGroup, $Y1Sn1, $Y2Sn1, $Y1Sn2, $Y2Sn2, $Y1Sn3, $Y2Sn3, $Y1Sn4, $Y2Sn4, $Sn1AtomNum, $Sn2AtomNum, $Sn3AtomNum, $Sn4AtomNum, $Sn1CarbonCount, $Sn2CarbonCount, $Sn3CarbonCount, $Sn4CarbonCount, $Sn2CAtomNum, $Sn2OAtomNum, $Sn2HAtomNum, $LMCategory, $LMMainClass, $LMSubClass, $CmpdString, $Sn1ReverseChainGrowth, $Sn2ReverseChainGrowth, $Sn3ReverseChainGrowth, $Sn4ReverseChainGrowth) = split /\|/, $TemplateData;
    $AbbrevTemplateDataMapRef->{HeadGroupName} = $HeadGroup;
  }

  $AbbrevTemplateDataMapRef->{AbbrevID} = $AbbrevID;

  @{$AbbrevTemplateDataMapRef->{SnCarbonCount}} = ();
  push @{$AbbrevTemplateDataMapRef->{SnCarbonCount}}, ($Sn1CarbonCount, $Sn2CarbonCount, $Sn3CarbonCount);
  if ($TemplateType =~ /^(CL|LM)$/) {
    push @{$AbbrevTemplateDataMapRef->{SnCarbonCount}}, ($Sn4CarbonCount);
  }

  $AbbrevTemplateDataMapRef->{TotalSnCarbons} = $Sn1CarbonCount + $Sn2CarbonCount + $Sn3CarbonCount;
  if ($TemplateType =~ /^(CL|LM)$/) {
    $AbbrevTemplateDataMapRef->{TotalSnCarbons} += $Sn4CarbonCount;
  }

  $AbbrevTemplateDataMapRef->{LMCategory} = $LMCategory;
  $AbbrevTemplateDataMapRef->{LMMainClass} = $LMMainClass;
  $AbbrevTemplateDataMapRef->{LMSubClass} = $LMSubClass;

  @{$AbbrevTemplateDataMapRef->{SnY1}} = ();
  push @{$AbbrevTemplateDataMapRef->{SnY1}}, ($Y1Sn1, $Y1Sn2, $Y1Sn3);
  if ($TemplateType =~ /^(CL|LM)$/) {
    push @{$AbbrevTemplateDataMapRef->{SnY1}}, ($Y1Sn4);
  }

  @{$AbbrevTemplateDataMapRef->{SnY2}} = ();
  push @{$AbbrevTemplateDataMapRef->{SnY2}}, ($Y2Sn1, $Y2Sn2, $Y2Sn3);
  if ($TemplateType =~ /^(CL|LM)$/) {
    push @{$AbbrevTemplateDataMapRef->{SnY2}}, ($Y2Sn4);
  }

  my($YSn1MidPoint, $YSn2MidPoint, $YSn3MidPoint);
  $YSn1MidPoint = ($Y1Sn1 + $Y2Sn1)/2;
  $YSn2MidPoint = ($Y1Sn2 + $Y2Sn2)/2;
  $YSn3MidPoint = ($Y1Sn3 + $Y2Sn3)/2;
  @{$AbbrevTemplateDataMapRef->{SnYMidPoints}} = ();
  push @{$AbbrevTemplateDataMapRef->{SnYMidPoints}}, ($YSn1MidPoint, $YSn2MidPoint, $YSn3MidPoint);
  if ($TemplateType =~ /^(CL|LM)$/) {
    my($YSn4MidPoint);
    $YSn4MidPoint = ($Y1Sn4 + $Y2Sn4)/2;
    push @{$AbbrevTemplateDataMapRef->{SnYMidPoints}}, ($YSn4MidPoint);
  }

  my(@SnAtomNums) = ($Sn1AtomNum, $Sn2AtomNum, $Sn3AtomNum);
  if ($TemplateType =~ /^(CL|LM)$/) {
    push @SnAtomNums, ($Sn4AtomNum);
  }
  @{$AbbrevTemplateDataMapRef->{SnAtomNums}} = ();
  push @{$AbbrevTemplateDataMapRef->{SnAtomNums}}, @SnAtomNums;

  my(@Sn2AtomNums) = ($Sn2CAtomNum, $Sn2OAtomNum, $Sn2HAtomNum);
  @{$AbbrevTemplateDataMapRef->{Sn2AtomNums}} = ();
  push @{$AbbrevTemplateDataMapRef->{Sn2AtomNums}}, @Sn2AtomNums;

  # Setup chain growth direcion for each Sn chain. The default is to grow chains from left to right...
  my($ReverseChainGrowth, @SnReverseChainGrowth);

  @{$AbbrevTemplateDataMapRef->{SnReverseChainGrowth}} = ();
  @SnReverseChainGrowth = ();
  for $ReverseChainGrowth ($Sn1ReverseChainGrowth, $Sn2ReverseChainGrowth, $Sn3ReverseChainGrowth, $Sn4ReverseChainGrowth) {
    push @SnReverseChainGrowth, ((defined($ReverseChainGrowth) && $ReverseChainGrowth) ? 1 : 0);
  }
  @{$AbbrevTemplateDataMapRef->{SnReverseChainGrowth}} = @SnReverseChainGrowth;

  my(@CmpdLines) = ();
  @CmpdLines = split /\n/, $CmpdString;

  @{$AbbrevTemplateDataMapRef->{CmpdLines}} = ();
  push @{$AbbrevTemplateDataMapRef->{CmpdLines}}, @CmpdLines;

  my($AtomCount, $BondCount) = LMAPSStr::ParseCmpdCountsLine($CmpdLines[3]);
  $AbbrevTemplateDataMapRef->{AtomCount} = $AtomCount;
  $AbbrevTemplateDataMapRef->{BondCount} = $BondCount;

  my($Index, @SnAtomLines, @SnChainAdd);
  @SnAtomLines = ();
  @SnChainAdd = ();
  push @SnChainAdd, @{$AbbrevTemplateDataMapRef->{SnChainAdd}};
  for $Index (0 .. $#SnChainAdd) {
    $SnAtomLines[$Index] = '';
    if ($SnChainAdd[$Index]) {
      # Line numbers start from 0 and atom number from 1...
      my($SnLineNum) = $SnAtomNums[$Index] + 4 - 1;
      $SnAtomLines[$Index] = $CmpdLines[$SnLineNum];
    }
  }
  @{$AbbrevTemplateDataMapRef->{SnAtomLines}} = ();
  push @{$AbbrevTemplateDataMapRef->{SnAtomLines}}, @SnAtomLines;

  my(@AtomBlockLines) = ();
  for $Index (4 .. ($AtomCount + 3)) {
    push @AtomBlockLines, $CmpdLines[$Index];
  }
  @{$AbbrevTemplateDataMapRef->{AtomBlockLines}} = ();
  push @{$AbbrevTemplateDataMapRef->{AtomBlockLines}}, @AtomBlockLines;

  my(@BondBlockLines) = ();
  my($FirstAtomNum, $SecondAtomNum, $BondType, $BondStereo, $BondLine, $AbbrevModifier);
  $AbbrevModifier = $AbbrevTemplateDataMapRef->{AbbrevModifier};
  for $Index ( ($AtomCount + 4) .. ($AtomCount + $BondCount + 3)) {
    # Parsing the bond block lines and regenerate it to make the bond lines consistent...
    ($FirstAtomNum, $SecondAtomNum, $BondType, $BondStereo) = LMAPSStr::ParseCmpdBondLine($CmpdLines[$Index]);
    if ($AbbrevModifier && $AbbrevModifier !~ /^R$/i) {
      # Modify bond stereochemistry designation between Sn2CAtom-Sn2OAtom and Sn2CAtom-Sn2HAtom...
      # Default: R => Sn2CAtom-Sn2OAtom : 6 (down); Sn2CAtom-Sn2HAtom: 1 (up)
      #
      # So, S => Sn2CAtom-Sn2OAtom : 1 (up); Sn2CAtom-Sn2HAtom: 6 (down)
      # rac/U => Sn2CAtom-Sn2OAtom : 0; Sn2CAtom-Sn2HAtom: 0
      #
      if ($FirstAtomNum =~ /^($Sn2CAtomNum|$Sn2OAtomNum)$/i && $SecondAtomNum =~ /^($Sn2CAtomNum|$Sn2OAtomNum)$/i) {
	$BondStereo = ($AbbrevModifier =~ /^S$/i) ? 1 : 0;
      }
      elsif ($FirstAtomNum =~ /^($Sn2CAtomNum|$Sn2HAtomNum)$/i && $SecondAtomNum =~ /^($Sn2CAtomNum|$Sn2HAtomNum)$/i) {
	$BondStereo = ($AbbrevModifier =~ /^S$/i) ? 6 : 0;
      }
    }
    $BondLine = LMAPSStr::GenerateCmpdBondLine($FirstAtomNum, $SecondAtomNum, $BondType, $BondStereo);
    push @BondBlockLines, $BondLine;
  }
  @{$AbbrevTemplateDataMapRef->{BondBlockLines}} = ();
  push @{$AbbrevTemplateDataMapRef->{BondBlockLines}}, @BondBlockLines;

  my(@DataBlockLines) = ();
  LINE: for $Index ( ($AtomCount + $BondCount + 4) .. $#CmpdLines) {
    push @DataBlockLines, $CmpdLines[$Index];
    if ($CmpdLines[$Index] =~ /M  END/) {
      last LINE;
    }
  }
  @{$AbbrevTemplateDataMapRef->{DataBlockLines}} = ();
  push @{$AbbrevTemplateDataMapRef->{DataBlockLines}}, @DataBlockLines;

  # Based on SnAbbrev, Setup the double bond position and geometry for each chain...
  my($ChainAbbrev, $ChainLength, $DoubleBondCount, $DoubleBondGeometry, $Position, $Geometry, $Substituents, $Rings);
  my(@SnChainLength, @SnDblBondCount, @SnDbleBondPosToGeometry, @SnSubstituentsCount, @SnSubstituentsBondCount, @SnSubstituentsPosInfo, @SnRing, @SnRingPosInfo, @SnAbbrev);

  @SnChainLength = ();
  @SnDblBondCount = ();
  @SnDbleBondPosToGeometry = ();
  @SnSubstituentsCount = ();
  @SnSubstituentsBondCount = ();
  @SnSubstituentsPosInfo = ();
  @SnRing = ();
  @SnRingPosInfo = ();

  push @SnAbbrev, @{$AbbrevTemplateDataMapRef->{SnAbbrev}};

  CHAIN: for $Index (0 .. $#SnChainAdd) {
    $SnChainLength[$Index] = 0;

    # Multiple bonds info...
    $SnDblBondCount[$Index] = 0;
    %{$SnDbleBondPosToGeometry[$Index]} = ();

    # Substituents info...
    $SnSubstituentsCount[$Index] = 0;
    $SnSubstituentsBondCount[$Index] = 0;
    %{$SnSubstituentsPosInfo[$Index]} = ();

    $SnRing[$Index] = 0;
    %{$SnRingPosInfo[$Index]} = ();

    if (!$SnChainAdd[$Index]) {
      next CHAIN;
    }
    $ChainAbbrev = $SnAbbrev[$Index];
    ($ChainLength, $DoubleBondCount, $DoubleBondGeometry, $Substituents, $Rings) = ChainAbbrev::ParseChainAbbrev($ChainAbbrev);

    if ($ChainLength =~ /^O-/) {
      # Alkyl and alkenyl chain formats: O-16:0 and O-16:1
      $ChainLength =~ s/^O-//;
    }
    elsif ($ChainLength =~ /^P-/) {
      # Alkenyl chain formats: P-16:0 corresponding to O-16:1
      $ChainLength =~ s/^P-//;
    }

    $SnChainLength[$Index] = $ChainLength;
    if (ChainAbbrev::IsAlkylChainAbbrev($ChainAbbrev)) {
      $SnDblBondCount[$Index] = 0;
    }
    elsif (ChainAbbrev::IsAlkenylChainAbbrevWithImplicitDoubleBond($ChainAbbrev)) {
      # Add implicit double bond at 1Z...
      $SnDblBondCount[$Index] = 1;
      $SnDbleBondPosToGeometry[$Index]{1} = "Z";

      # Add any other double bond specifications...
      my(@DoubleBondGeometryList) = split /,/, $DoubleBondGeometry;
      GEOMETRY: for $DoubleBondGeometry (@DoubleBondGeometryList) {
	($Position, $Geometry) = ($DoubleBondGeometry =~ /^([0-9]+)([a-zA-Z]+)$/);
	if ($Position == 1) {
	  next GEOMETRY;
	}
	$SnDblBondCount[$Index] += 1;
	$SnDbleBondPosToGeometry[$Index]{$Position} = $Geometry;
      }

    }
    else {
      $SnDblBondCount[$Index] = $DoubleBondCount;
      my(@DoubleBondGeometryList) = split /,/, $DoubleBondGeometry;
      for $DoubleBondGeometry (@DoubleBondGeometryList) {
	($Position, $Geometry) = ($DoubleBondGeometry =~ /^([0-9]+)([a-zA-Z]+)$/);
	$SnDbleBondPosToGeometry[$Index]{$Position} = $Geometry;
      }
    }
    # Setup substituents information...
    if (!IsEmpty($Substituents)) {
      my($Substituent, $SubstituentPos, $ValidStartSubstituentPos, $SubstituentAbbrev, $SubstituentStereoChemistry, @SubstituentsList);
      $ValidStartSubstituentPos = $AbbrevTemplateDataMapRef->{SnCarbonCount}[$Index];
      (@SubstituentsList) = split /\,/, $Substituents;
      SUBSTITUENT: for $Substituent (@SubstituentsList) {
	($SubstituentPos, $SubstituentAbbrev, $SubstituentStereoChemistry) = ChainAbbrev::ParseSubstituentAbbrev($Substituent);
	if ($SubstituentPos <=  $ValidStartSubstituentPos) {
	  # Sustituent position number for a chain must be > number of chain atoms in the template;
	  # otherwise # substituent cann't be added during structure generation. The specific template
	  # must be modified.
	  #
	  print "Warning: Ignoring specified substituent $Substituent: Substituent position number $SubstituentPos, is not valid for chain number " .  ($Index + 1) . ". During structure generation for template type $TemplateType and template ID $AbbrevID, the substituent position must be > $ValidStartSubstituentPos, the number of template chain atoms.\n";
	  next SUBSTITUENT;
	}

	$SnSubstituentsCount[$Index] += 1;
	$SnSubstituentsBondCount[$Index] += 1;

	if ($SubstituentAbbrev =~ /^[CE]p$/i) {
	  # One more bond for epoxy..
	  $SnSubstituentsBondCount[$Index] += 1;
	}
	if ($SubstituentStereoChemistry) {
	  $SubstituentStereoChemistry = LMAPSStr::StandardizeStereochemistrySpec($SubstituentStereoChemistry);
	}
	if (exists $SnSubstituentsPosInfo[$Index]{$SubstituentPos}) {
	  # Multiple substituents at one position...
	  $SnSubstituentsPosInfo[$Index]{$SubstituentPos}{Count} += 1;
	  push @{$SnSubstituentsPosInfo[$Index]{$SubstituentPos}{Spec}}, $Substituent;
	  push @{$SnSubstituentsPosInfo[$Index]{$SubstituentPos}{Abbrev}}, $SubstituentAbbrev;
	  push @{$SnSubstituentsPosInfo[$Index]{$SubstituentPos}{StereoChemistry}}, $SubstituentStereoChemistry;
	}
	else {
	  %{$SnSubstituentsPosInfo[$Index]{$SubstituentPos}} = ();

	  $SnSubstituentsPosInfo[$Index]{$SubstituentPos}{Count} = 1;
	  @{$SnSubstituentsPosInfo[$Index]{$SubstituentPos}{Spec}} = ();
	  push @{$SnSubstituentsPosInfo[$Index]{$SubstituentPos}{Spec}}, $Substituent;

	  @{$SnSubstituentsPosInfo[$Index]{$SubstituentPos}{Abbrev}} = ();
	  push @{$SnSubstituentsPosInfo[$Index]{$SubstituentPos}{Abbrev}}, $SubstituentAbbrev;

	  @{$SnSubstituentsPosInfo[$Index]{$SubstituentPos}{StereoChemistry}} = ();
	  push @{$SnSubstituentsPosInfo[$Index]{$SubstituentPos}{StereoChemistry}}, $SubstituentStereoChemistry;

	  # Set up the place holders for chiral center and substituent atom number. Chiral center atom
	  # number corresponds to what's...
	  @{$SnSubstituentsPosInfo[$Index]{$SubstituentPos}{ChiralCenterAtomNum}} = ();
	  @{$SnSubstituentsPosInfo[$Index]{$SubstituentPos}{SubstituentAtomNum}} = ();
	}
      }
    }
    # Setup rings information...
    if (!IsEmpty($Rings)) {
      my($RingAbbrev, $StartRingPos, $EndRingPos, $RingSize, $StartRingStereoChemistry, $EndRingStereoChemistry, $StartRingAbbrev, $EndRingAbbrev, @RingWords);
      $SnRing[$Index] = 1;
      @RingWords = quotewords(',', 0, $Rings);

      # Start ring position...
      $StartRingAbbrev = $RingWords[0];
      ($StartRingPos, $StartRingStereoChemistry) = ChainAbbrev::ParseRingAbbrev($StartRingAbbrev);
      $StartRingStereoChemistry = LMAPSStr::StandardizeRingStereochemistrySpec($StartRingStereoChemistry);
      $SnRingPosInfo[$Index]{StartRingPos} = $StartRingPos;
      $SnRingPosInfo[$Index]{StartRingStereoChemistry} = $StartRingStereoChemistry;

      # End ring position...
      $EndRingAbbrev = $RingWords[1];
      ($EndRingPos, $EndRingStereoChemistry) = ChainAbbrev::ParseRingAbbrev($EndRingAbbrev);
      $EndRingStereoChemistry = LMAPSStr::StandardizeRingStereochemistrySpec($EndRingStereoChemistry);
      $SnRingPosInfo[$Index]{EndRingPos} = $EndRingPos;
      $SnRingPosInfo[$Index]{EndRingStereoChemistry} = $EndRingStereoChemistry;

      $RingSize = $EndRingPos - $StartRingPos + 1;
      $SnRingPosInfo[$Index]{RingSize} = $RingSize;

      # Setup ring positions to modify, offset for each position relative to previoud position and for reversing direction of structure growth...
      my($RingPos);
      %{$SnRingPosInfo[$Index]{RingPosToModify}} = ();
      %{$SnRingPosInfo[$Index]{PosXOffSet}} = ();
      %{$SnRingPosInfo[$Index]{PosYOffSet}} = ();
      %{$SnRingPosInfo[$Index]{ReverseStructureGrowth}} = ();
      if ($RingSize == 5) {
	# Tasks:
	# . In addition to all marking all ring positions, the atom position after the last ring position as well;
	#   otherwise, five member ring doesn't fit.
	# . Set offset for each position is set relative to the previous position.
	# . Offset for last three positions is the same as first three positions except for x and y direction.
	#
	my($XOffSet, $YOffSet);
	$XOffSet = 0.7200;
	$YOffSet = $AbbrevTemplateDataMapRef->{SnY2}[$Index] - $AbbrevTemplateDataMapRef->{SnY1}[$Index];
	$YOffSet = abs $YOffSet;
	for $RingPos ($StartRingPos ... ($StartRingPos + 5)) {
	  $SnRingPosInfo[$Index]{RingPosToModify}{$RingPos} = $RingPos;
	  if ( ($RingPos == $StartRingPos) || ($RingPos == ($StartRingPos + 5)) ) {
	    $SnRingPosInfo[$Index]{PosXOffSet}{$RingPos} = $XOffSet  - 0.0470;
	    $SnRingPosInfo[$Index]{PosYOffSet}{$RingPos} = $YOffSet  + 0.0724;
	  }
	  elsif ( ($RingPos == ($StartRingPos + 1))  || ($RingPos == ($StartRingPos + 4)) ) {
	    $SnRingPosInfo[$Index]{PosXOffSet}{$RingPos} = $XOffSet  + 0.0702;
	    $SnRingPosInfo[$Index]{PosYOffSet}{$RingPos} = $YOffSet  - 0.1576;
	  }
	  elsif ( ($RingPos == ($StartRingPos + 2))  || ($RingPos == ($StartRingPos + 3)) ) {
	    $SnRingPosInfo[$Index]{PosXOffSet}{$RingPos} = $XOffSet  - 0.2285;
	    $SnRingPosInfo[$Index]{PosYOffSet}{$RingPos} = $YOffSet  + 0.2549;
	  }
	}
	$RingPos = $EndRingPos + 1;
	$SnRingPosInfo[$Index]{ReverseStructureGrowth}{$RingPos} = $RingPos;
      }
      elsif ($RingSize == 6) {
	my($XOffSet, $YOffSet);
	$XOffSet = 0.7200;
	$YOffSet = $AbbrevTemplateDataMapRef->{SnY2}[$Index] - $AbbrevTemplateDataMapRef->{SnY1}[$Index];
	$YOffSet = abs $YOffSet;

	for $RingPos ($StartRingPos ... ($StartRingPos + 6)) {
	  $SnRingPosInfo[$Index]{RingPosToModify}{$RingPos} = $RingPos;
	  if ( $RingPos == ($StartRingPos + 3) ) {
	    $SnRingPosInfo[$Index]{PosXOffSet}{$RingPos} = 0;
	    $SnRingPosInfo[$Index]{PosYOffSet}{$RingPos} = $YOffSet + 0.3376;
	  }
	  else {
	    $SnRingPosInfo[$Index]{PosXOffSet}{$RingPos} = $XOffSet;
	    $SnRingPosInfo[$Index]{PosYOffSet}{$RingPos} = $YOffSet;
	  }
	}
	$RingPos = $EndRingPos + 1;
	$SnRingPosInfo[$Index]{ReverseStructureGrowth}{$RingPos} = $RingPos;
      }
    }
  }
  @{$AbbrevTemplateDataMapRef->{SnChainLength}} = ();
  push @{$AbbrevTemplateDataMapRef->{SnChainLength}}, @SnChainLength;

  @{$AbbrevTemplateDataMapRef->{SnDblBondCount}} = ();
  push @{$AbbrevTemplateDataMapRef->{SnDblBondCount}}, @SnDblBondCount;

  @{$AbbrevTemplateDataMapRef->{SnDbleBondPosToGeometry}} = ();
  push @{$AbbrevTemplateDataMapRef->{SnDbleBondPosToGeometry}}, @SnDbleBondPosToGeometry;

  @{$AbbrevTemplateDataMapRef->{SnSubstituentsCount}} = ();
  push @{$AbbrevTemplateDataMapRef->{SnSubstituentsCount}}, @SnSubstituentsCount;

  @{$AbbrevTemplateDataMapRef->{SnSubstituentsBondCount}} = ();
  push @{$AbbrevTemplateDataMapRef->{SnSubstituentsBondCount}}, @SnSubstituentsBondCount;

  @{$AbbrevTemplateDataMapRef->{SnSubstituentsPosInfo}} = ();
  push @{$AbbrevTemplateDataMapRef->{SnSubstituentsPosInfo}}, @SnSubstituentsPosInfo;

  @{$AbbrevTemplateDataMapRef->{SnRing}} = ();
  push @{$AbbrevTemplateDataMapRef->{SnRing}}, @SnRing;

  @{$AbbrevTemplateDataMapRef->{SnRingPosInfo}} = ();
  push @{$AbbrevTemplateDataMapRef->{SnRingPosInfo}}, @SnRingPosInfo;
}

# Assign bond stereochemistry to substituent bond: up or down...
sub _AssignStereoChemistry {
  my($CmpdAbbrevTemplateDataMapRef, $CmdDataLinesRef, $ChiralCenterAtomNum, $SubstituentAtomNum, $SubstituentStereoChemistry) = @_;

  my($AtomsDataRef);
  # Retrieve atom data including bonds, atom symbol and coordinates...
  $AtomsDataRef = _GetAtomDataFromCompoundDataLines($CmdDataLinesRef);

  # To assign chirality appropriately, ChiralCenterAtom must have only four or three (fourth one being implicit hydrogen)
  # substituent without any multiple bonds...
  if (!exists($AtomsDataRef->{BondAtomNums}{$ChiralCenterAtomNum})) {
    print "Warning: Internal inconsistency: Skipping chirality assignment...\n";
    return;
  }

  my($ChiralCenterNeighborsCount, $ChiralCenterMultipleBondCount, $NeighborAtomNum);
  $ChiralCenterNeighborsCount = scalar @{$AtomsDataRef->{BondAtomNums}{$ChiralCenterAtomNum}};
  $ChiralCenterMultipleBondCount = 0;
  for $NeighborAtomNum (@{$AtomsDataRef->{BondAtomNums}{$ChiralCenterAtomNum}}) {
    if ($AtomsDataRef->{BondType}{$ChiralCenterAtomNum}{$NeighborAtomNum} != 1) {
      $ChiralCenterMultipleBondCount++;
    }
  }
  if (!(($ChiralCenterNeighborsCount == 3 && $ChiralCenterMultipleBondCount == 0) || $ChiralCenterNeighborsCount == 4)) {
    # No chirality assignment...
    print "Warning: Skipping chirality assignment: Invalid number of substituents or multiple bond types around chiral center...\n";
    return;
  }

  # Initialize chirality info before going into a recursive loop for specific number
  # of neighborhood levels to resolve chirality assignment...
  my($ChiralityDataRef);
  $ChiralityDataRef = _InitializeChiralityData($AtomsDataRef, $ChiralCenterAtomNum, $SubstituentAtomNum, $SubstituentStereoChemistry);

  # Assign ligand preceence at zero for immediate neighbors...
  _AssignLigandsPrecedence($AtomsDataRef, $ChiralityDataRef);
  $ChiralityDataRef->{CurrentTopologicalDistance} += 1;

  _DetermineChirality($AtomsDataRef, $ChiralityDataRef);

  if ($ChiralityDataRef->{ChiralityAssignmentOkay}) {
    # Update connect record for substituent atom/chiral center atom to reflect
    # substituent bond value...
     _UpdateSubstituentBondType($CmdDataLinesRef, $ChiralityDataRef);
  }
  else {
    print "Warning: Skipping chirality assignment: Couldn't determine chirality after exploring upto topolgical distance ",  $ChiralityDataRef->{CurrentTopologicalDistance}," from chiral center...\n";
  }

}

# Assign and check precedence to resolve any conflicts in ranking of ligand groups around chiral center
# for chirality assignment before moving on to get next set of neighbors and repeating the process
# recursively upto a maximum neighorhood topological distance.
#
# CIP sequence rule applied:
# . Ligands with higher atomic mass number preceeds lower
# . Ligands with higher multiple bond order and atomic mass number similar preceeds ones with lower multiple bond order
# .
sub _DetermineChirality {
  my($AtomsDataRef, $ChiralityDataRef) = @_;

  #print "\nTopological distance: $ChiralityDataRef->{CurrentTopologicalDistance} \n";
  if ($ChiralityDataRef->{CurrentTopologicalDistance} > $ChiralityDataRef->{MaxTopologicalDistance}) {
    return;
  }
  if ($ChiralityDataRef->{Done}) {
    return;
  }
  _AssignLigandsPrecedence($AtomsDataRef, $ChiralityDataRef);

  _DetermineChiralityUsingLigandsPrecedence($AtomsDataRef, $ChiralityDataRef);

  $ChiralityDataRef->{CurrentTopologicalDistance} += 1;

  # Get next set of neighbors...
  _GetNextTopologicalNeighbors($AtomsDataRef, $ChiralityDataRef);

  # And recursively go over the next set of neighbors...
  _DetermineChirality($AtomsDataRef, $ChiralityDataRef);
}

# Assign ligands precedence using atomic masses and bond orders...
sub _AssignLigandsPrecedence {
  my($AtomsDataRef, $ChiralityDataRef) = @_;

  my($NbrAtomIndex, $NbrAtomNum, $CurrentTopologicalDistance, $CurrentDistanceIndex, $BondedAtomNum);

  $CurrentTopologicalDistance = $ChiralityDataRef->{CurrentTopologicalDistance};
  $CurrentDistanceIndex = $CurrentTopologicalDistance;

  # Go over the ligand positions around chiral center...
  for $NbrAtomIndex (0 .. $ChiralityDataRef->{ChiralCenterNbrsMaxIndex}) {
    # Go over the neighbors at each position...
    for $NbrAtomNum (@{$ChiralityDataRef->{NbrsAtomNums}[$NbrAtomIndex]}) {
      if ($CurrentTopologicalDistance == 0) {
	# Use the atomic weights of the atoms next to central chiral center...
	$ChiralityDataRef->{ChiralCenterNbrsScores}[$NbrAtomIndex][$CurrentDistanceIndex] += LMAPSStr::RoundToNextInteger($AtomsDataRef->{AtomicWeight}{$NbrAtomNum});
      }
      else {
	# Go over the bonded atoms of $NbrAtomNum and scale atomic weights by bond order...
	for $BondedAtomNum (@{$AtomsDataRef->{BondAtomNums}{$NbrAtomNum}}) {
	  $ChiralityDataRef->{ChiralCenterNbrsScores}[$NbrAtomIndex][$CurrentDistanceIndex] += LMAPSStr::RoundToNextInteger($AtomsDataRef->{AtomicWeight}{$BondedAtomNum} * $AtomsDataRef->{BondType}{$BondedAtomNum}{$NbrAtomNum});
	}
      }
    }
  }

}

# Assign ligands precedence using atomic masses and bond orders...
sub _DetermineChiralityUsingLigandsPrecedence {
  my($AtomsDataRef, $ChiralityDataRef) = @_;

  #print "\n_DetermineChiralityUsingLigandsPrecedence: Toplogical distance: $ChiralityDataRef->{CurrentTopologicalDistance}...\n";
  my($NbrAtomIndex, $DistanceIndex, $TotalNbrsScore, $NbrsScore, $ChiralCenterNbrAtomNum);

  %{$ChiralityDataRef->{ChiralCenterNbrsTotalScoresToAtomNumMap}} = ();
  for $NbrAtomIndex (0 .. $ChiralityDataRef->{ChiralCenterNbrsMaxIndex}) {
    $ChiralCenterNbrAtomNum = $ChiralityDataRef->{ChiralCenterNbrsAtomNums}[$NbrAtomIndex];
    $TotalNbrsScore = '';
    for $DistanceIndex (0 .. ($ChiralityDataRef->{CurrentTopologicalDistance} - 1)) {
      # Convert neighbors score into six digit integer with 0 as left padding...
      $NbrsScore = sprintf "%06.6s", $ChiralityDataRef->{ChiralCenterNbrsScores}[$NbrAtomIndex][$DistanceIndex];
      $TotalNbrsScore .= $NbrsScore;
    }
    #print "$ChiralCenterNbrAtomNum - $TotalNbrsScore\n";
    if (exists $ChiralityDataRef->{ChiralCenterNbrsTotalScoresToAtomNumMap}{$TotalNbrsScore}) {
      # Same score for ligand groups. So chirality determination is not possible at the current topological distance.
      return;
    }
    $ChiralityDataRef->{ChiralCenterNbrsTotalScoresToAtomNumMap}{$TotalNbrsScore} = $ChiralCenterNbrAtomNum;
  }

  # Sort chiral center neighbor atom number from highest to lowest using score...
  @{$ChiralityDataRef->{ChiralCenterNbrsSortedAtomNums}} = ();
  for $TotalNbrsScore (sort {$b cmp $a} keys %{$ChiralityDataRef->{ChiralCenterNbrsTotalScoresToAtomNumMap}}) {
    $ChiralCenterNbrAtomNum = $ChiralityDataRef->{ChiralCenterNbrsTotalScoresToAtomNumMap}{$TotalNbrsScore};
    push @{$ChiralityDataRef->{ChiralCenterNbrsSortedAtomNums}}, $ChiralCenterNbrAtomNum;
    #print "$ChiralCenterNbrAtomNum\n";
  }

  # Clockwise: ClockPosition > 0; Anti Clockwise: ClockPosition < 0; On a line: ClockPosition: 0
   _DetermineChiralCenterNbrsClockPosition($AtomsDataRef, $ChiralityDataRef);
  if ($ChiralityDataRef->{ChiralCenterNbrsClockPosition} == 0) {
    # No chirality determination at this step...
    return;
  }

  _AssignSubstituentBondPosition($AtomsDataRef, $ChiralityDataRef);

}

# Based on specified stereo chemistry and clock position of chiral center neighbors,
# assign substituents bond position for SD files: 1 - up; 6 - down.
#
# Notes:
#  Clockwise:  R - up; S - down;
#  Anticlockwise: R - down; S - up
#
sub _AssignSubstituentBondPosition {
  my($AtomsDataRef, $ChiralityDataRef) = @_;

  if ($ChiralityDataRef->{ChiralCenterNbrsClockPosition} > 0) {
    # Clockwise...
    #print "_AssignSubstituentBondPosition: Clockwise...\n";
    $ChiralityDataRef->{SubstituentBondTypeAssigned} = ($ChiralityDataRef->{SubstituentStereoChmistry} =~ /^R$/i) ? 1 : 6;
    $ChiralityDataRef->{ChiralityAssignmentOkay} = 1;
    $ChiralityDataRef->{Done} = 1;
  }
  elsif ($ChiralityDataRef->{ChiralCenterNbrsClockPosition} < 0) {
    # Anti clockwise...
    #print "_AssignSubstituentBondPosition: Anticlockwise...\n";
    $ChiralityDataRef->{SubstituentBondTypeAssigned} = ($ChiralityDataRef->{SubstituentStereoChmistry} =~ /^R$/i) ? 6 : 1;
    $ChiralityDataRef->{ChiralityAssignmentOkay} = 1;
    $ChiralityDataRef->{Done} = 1;
  }
  else {
    return;
  }
}


# Figure out the clock position of using first three atom in sorted chiral center neighbor list...
# Clockwise: ClockPosition > 0; Anti Clockwise: ClockPosition < 0; On a line: ClockPosition: 0
sub _DetermineChiralCenterNbrsClockPosition {
  my($AtomsDataRef, $ChiralityDataRef) = @_;
  my($FirstAtomNum, $SecondAtomNum, $ThirdAtomNum, $X1, $Y1, $X2, $Y2, $X3, $Y3);

  $ChiralityDataRef->{ChiralCenterNbrsClockPosition} = 0;
  if (@{$ChiralityDataRef->{ChiralCenterNbrsSortedAtomNums}} < 3) {
    return;
  }
  $FirstAtomNum = $ChiralityDataRef->{ChiralCenterNbrsSortedAtomNums}[0];
  $SecondAtomNum = $ChiralityDataRef->{ChiralCenterNbrsSortedAtomNums}[1];
  $ThirdAtomNum = $ChiralityDataRef->{ChiralCenterNbrsSortedAtomNums}[2];

  $X1 = $AtomsDataRef->{XYZ}{$FirstAtomNum}[0]; $Y1 = $AtomsDataRef->{XYZ}{$FirstAtomNum}[1];
  $X2 = $AtomsDataRef->{XYZ}{$SecondAtomNum}[0]; $Y2 = $AtomsDataRef->{XYZ}{$SecondAtomNum}[1];
  $X3 = $AtomsDataRef->{XYZ}{$ThirdAtomNum}[0]; $Y3 = $AtomsDataRef->{XYZ}{$ThirdAtomNum}[1];

  $ChiralityDataRef->{ChiralCenterNbrsClockPosition} = ($X3 - $X1) * ($Y2 - $Y1) - ($X2 - $X1) * ($Y3 - $Y1);
}

# Get next set of toplogical neigbors for each position around chiral center...
sub _GetNextTopologicalNeighbors {
  my($AtomsDataRef, $ChiralityDataRef) = @_;
  my($NbrAtomNum, $NextNbrAtomNum, $NbrAtomIndex);

  # Mark current neighbors already processed...
  for $NbrAtomIndex (0 .. $ChiralityDataRef->{ChiralCenterNbrsMaxIndex}) {
    for $NbrAtomNum (@{$ChiralityDataRef->{NbrsAtomNums}[$NbrAtomIndex]}) {
      $ChiralityDataRef->{AtomNumsNbrsAlreadyProcessed}{$NbrAtomNum} = $NbrAtomNum;
    }
  }

  # Get next set of neigbors without the earlier atoms...
  for $NbrAtomIndex (0 .. $ChiralityDataRef->{ChiralCenterNbrsMaxIndex}) {
    @{$ChiralityDataRef->{NextNbrsAtomNums}[$NbrAtomIndex]} = ();
    for $NbrAtomNum (@{$ChiralityDataRef->{NbrsAtomNums}[$NbrAtomIndex]}) {
      for $NextNbrAtomNum (@{$AtomsDataRef->{BondAtomNums}{$NbrAtomNum}}) {
	if (! exists $ChiralityDataRef->{AtomNumsNbrsAlreadyProcessed}{$NextNbrAtomNum}) {
	  push @{$ChiralityDataRef->{NextNbrsAtomNums}[$NbrAtomIndex]}, $NextNbrAtomNum;
	}
      }
    }
  }

  # Copy next set of neighbors back to neighbors...
  for $NbrAtomIndex (0 .. $ChiralityDataRef->{ChiralCenterNbrsMaxIndex}) {
    @{$ChiralityDataRef->{NbrsAtomNums}[$NbrAtomIndex]} = ();
    push @{$ChiralityDataRef->{NbrsAtomNums}[$NbrAtomIndex]}, @{$ChiralityDataRef->{NextNbrsAtomNums}[$NbrAtomIndex]};
  }
}

# Update substituent connection data to reflect new bond type assignment...
sub _UpdateSubstituentBondType {
  my($CmdDataLinesRef, $ChiralityDataRef) = @_;
  my($ChiralCenterAtomNum, $SubstituentAtomNum, $SubstituentBondStereoType, $LineIndex, $AtomCount, $BondCount, $FirstAtomNum, $SecondAtomNum, $BondType, $BondStereo);

  $ChiralCenterAtomNum = $ChiralityDataRef->{ChiralCenterAtomNum};
  $SubstituentAtomNum = $ChiralityDataRef->{SubstituentAtomNum};
  $SubstituentBondStereoType = $ChiralityDataRef->{SubstituentBondTypeAssigned};

  ($AtomCount, $BondCount) = LMAPSStr::ParseCmpdCountsLine($CmdDataLinesRef->[3]);

  # Go over the bonds data to find $SubstituentAtomNum...
  LINE: for ($LineIndex = 4 + $AtomCount; $LineIndex < (4 + $AtomCount + $BondCount); $LineIndex++) {
    ($FirstAtomNum, $SecondAtomNum, $BondType, $BondStereo) = LMAPSStr::ParseCmpdBondLine($CmdDataLinesRef->[$LineIndex]);
    if (($FirstAtomNum == $SubstituentAtomNum && $SecondAtomNum == $ChiralCenterAtomNum) || ($FirstAtomNum == $ChiralCenterAtomNum && $SecondAtomNum == $SubstituentAtomNum)) {
      $BondStereo = $SubstituentBondStereoType;
      $CmdDataLinesRef->[$LineIndex] = LMAPSStr::GenerateCmpdBondLine($FirstAtomNum, $SecondAtomNum, $BondType, $BondStereo);
      last LINE;
    }
  }
}

# Initialize chirality information...
sub _InitializeChiralityData {
  my($AtomsDataRef, $ChiralCenterAtomNum, $SubstituentAtomNum, $SubstituentStereoChemistry) = @_;
  my(%ChiralityDataMap);

  %ChiralityDataMap = ();

  $ChiralityDataMap{ChiralCenterAtomNum} = $ChiralCenterAtomNum;
  $ChiralityDataMap{SubstituentAtomNum} = $SubstituentAtomNum;
  $ChiralityDataMap{SubstituentStereoChmistry} = $SubstituentStereoChemistry;

  # Maximum neighborhood level to explore for resolving conflicts in chirality assignment...
  $ChiralityDataMap{MaxTopologicalDistance} = 20;
  $ChiralityDataMap{CurrentTopologicalDistance} = 0;

  # Flag to indicate either chirality was assigned or had to give up..
  $ChiralityDataMap{Done} = 0;

  # Flag to indicate whether chirality was successfully ascertained...
  $ChiralityDataMap{ChiralityAssignmentOkay} = 0;

  # Sustituent bond type assigned: 0 - no assignment; 1 - up; 6 - down.
  $ChiralityDataMap{SubstituentBondTypeAssigned} = 0;

  # Setup neighbors of chiral center...
  @{$ChiralityDataMap{ChiralCenterNbrsAtomNums}} = ();
  push @{$ChiralityDataMap{ChiralCenterNbrsAtomNums}}, @{$AtomsDataRef->{BondAtomNums}{$ChiralCenterAtomNum}};
  $ChiralityDataMap{ChiralCenterNbrsCount} = scalar @{$ChiralityDataMap{ChiralCenterNbrsAtomNums}};
  $ChiralityDataMap{ChiralCenterNbrsMaxIndex} =  $ChiralityDataMap{ChiralCenterNbrsCount} - 1;

  $ChiralityDataMap{ChiralCenterNbrsClockPosition} = 0;

  # To hold scores for each chiral center neighbore at different neighborhood levels...
  my($NbrAtomIndex);
  @{$ChiralityDataMap{ChiralCenterNbrsScores}} = ();
  for $NbrAtomIndex (0 .. $ChiralityDataMap{ChiralCenterNbrsMaxIndex}) {
    # Each neighbor atom can hold score for all possible levels...
    @{$ChiralityDataMap{ChiralCenterNbrsScores}[$NbrAtomIndex]} = ();
    push @{$ChiralityDataMap{ChiralCenterNbrsScores}[$NbrAtomIndex]}, (0) x $ChiralityDataMap{MaxTopologicalDistance};
  }
  %{$ChiralityDataMap{ChiralCenterNbrsTotalScoresToAtomNumMap}} = ();
  @{$ChiralityDataMap{ChiralCenterNbrsSortedAtomNums}} = ();

  # Setup arrays for holding neighbors at any level for downstream neighbors of chiral center neighbors...
  @{$ChiralityDataMap{NbrsAtomNums}} = ();
  for $NbrAtomIndex (0 .. $ChiralityDataMap{ChiralCenterNbrsMaxIndex}) {
    @{$ChiralityDataMap{NbrsAtomNums}[$NbrAtomIndex]} = ();
  }
  @{$ChiralityDataMap{NextNbrsAtomNums}} = ();
  for $NbrAtomIndex (0 .. $ChiralityDataMap{ChiralCenterNbrsMaxIndex}) {
    @{$ChiralityDataMap{NextNbrsAtomNums}[$NbrAtomIndex]} = ();
  }

  # Set up a map to hold atom numbers whose neighbors has already been processed: It'll also
  # help not to go back to the previous atom numbers during neighbors collection...
  %{$ChiralityDataMap{AtomNumsNbrsAlreadyProcessed}} = ();

  # Assign the neighbors of chiral center to start things off...
  for $NbrAtomIndex (0 .. $ChiralityDataMap{ChiralCenterNbrsMaxIndex}) {
    push @{$ChiralityDataMap{NbrsAtomNums}[$NbrAtomIndex]}, $ChiralityDataMap{ChiralCenterNbrsAtomNums}[$NbrAtomIndex];
  }
  $ChiralityDataMap{AtomNumsNbrsAlreadyProcessed}{$ChiralCenterAtomNum} = $ChiralCenterAtomNum;

  return \%ChiralityDataMap;

}

# Setup a hash to facilitate access to atom data...
sub _GetAtomDataFromCompoundDataLines {
  my($CmdDataLinesRef) = @_;
  my($LineIndex, $AtomCount, $BondCount, $AtomX, $AtomY, $AtomZ, $AtomSymbol, $AtomNum, $FirstAtomNum, $SecondAtomNum, $BondType, $BondStereo, %AtomsDataMap);

  %AtomsDataMap = ();
  @{$AtomsDataMap{AtomNums}} = ();
  %{$AtomsDataMap{XYZ}} = ();
  %{$AtomsDataMap{AtomSymbol}} = ();
  %{$AtomsDataMap{AtomicWeight}} = ();
  %{$AtomsDataMap{BondAtomNums}} = ();
  %{$AtomsDataMap{BondType}} = ();
  %{$AtomsDataMap{BondStereo}} = ();

  ($AtomCount, $BondCount) = LMAPSStr::ParseCmpdCountsLine($CmdDataLinesRef->[3]);

  # Process atom data...
  $AtomNum = 0;
  for ($LineIndex = 4; $LineIndex < (4 + $AtomCount); $LineIndex++) {
    $AtomNum++;

    # Read atom data...
    ($AtomX, $AtomY, $AtomZ, $AtomSymbol) = LMAPSStr::ParseCmpdAtomLine($CmdDataLinesRef->[$LineIndex]);

    # Save atom data...
    push @{$AtomsDataMap{AtomNums}}, $AtomNum;

    @{$AtomsDataMap{XYZ}{$AtomNum}} = ();
    push @{$AtomsDataMap{XYZ}{$AtomNum}}, ($AtomX, $AtomY, $AtomZ);

    $AtomsDataMap{AtomSymbol}{$AtomNum} = $AtomSymbol;
    $AtomsDataMap{AtomicWeight}{$AtomNum} = _GetAtomicWeightFromAtomSymbol($AtomSymbol);
  }

  # Process bonds data...
  for ($LineIndex = 4 + $AtomCount; $LineIndex < (4 + $AtomCount + $BondCount); $LineIndex++) {
    # Read bond data lines...
    ($FirstAtomNum, $SecondAtomNum, $BondType, $BondStereo) = LMAPSStr::ParseCmpdBondLine($CmdDataLinesRef->[$LineIndex]);

    # Store bond data for $FirstAtomNum...
    if (exists $AtomsDataMap{BondType}{$FirstAtomNum}) {
      push @{$AtomsDataMap{BondAtomNums}{$FirstAtomNum}}, $SecondAtomNum;

      $AtomsDataMap{BondType}{$FirstAtomNum}{$SecondAtomNum} = $BondType;
      $AtomsDataMap{BondStereo}{$FirstAtomNum}{$SecondAtomNum} = $BondStereo;
    }
    else {
      # New bond data for $FirstAtomNum...
      @{$AtomsDataMap{BondAtomNums}{$FirstAtomNum}} = ();
      push @{$AtomsDataMap{BondAtomNums}{$FirstAtomNum}}, $SecondAtomNum;

      %{$AtomsDataMap{BondType}{$FirstAtomNum}} = ();
      $AtomsDataMap{BondType}{$FirstAtomNum}{$SecondAtomNum} = $BondType;

      %{$AtomsDataMap{BondStereo}{$FirstAtomNum}} = ();
      $AtomsDataMap{BondStereo}{$FirstAtomNum}{$SecondAtomNum} = $BondStereo;
    }

    # Store bond data for $SecondAtomNum...
    if (exists $AtomsDataMap{BondType}{$SecondAtomNum}) {
      push @{$AtomsDataMap{BondAtomNums}{$SecondAtomNum}}, $FirstAtomNum;

      $AtomsDataMap{BondType}{$SecondAtomNum}{$FirstAtomNum} = $BondType;
      $AtomsDataMap{BondStereo}{$SecondAtomNum}{$FirstAtomNum} = $BondStereo;
    }
    else {
      # New bond data for $FirstAtomNum...
      @{$AtomsDataMap{BondAtomNums}{$SecondAtomNum}} = ();
      push @{$AtomsDataMap{BondAtomNums}{$SecondAtomNum}}, $FirstAtomNum;

      %{$AtomsDataMap{BondType}{$SecondAtomNum}} = ();
      $AtomsDataMap{BondType}{$SecondAtomNum}{$FirstAtomNum} = $BondType;

      %{$AtomsDataMap{BondStereo}{$SecondAtomNum}} = ();
      $AtomsDataMap{BondStereo}{$SecondAtomNum}{$FirstAtomNum} = $BondStereo;
    }
  }
  return \%AtomsDataMap;
}

# Get atomic weight from atom symbols used in SD file.
#
# Note:
#  . Atom symbols can correspond to substituents other than element symbols...
#
sub _GetAtomicWeightFromAtomSymbol {
  my($AtomSymbol) = @_;
  my($AtomicWeight) = undef;

  WEIGHT: {
    if ($AtomSymbol =~ /^H$/i) { $AtomicWeight = 1.00794; last WEIGHT;}
    if ($AtomSymbol =~ /^C$/i) { $AtomicWeight = 12.0107; last WEIGHT;}
    if ($AtomSymbol =~ /^N$/i) { $AtomicWeight = 14.0067; last WEIGHT;}
    if ($AtomSymbol =~ /^O$/i) { $AtomicWeight = 15.9994; last WEIGHT;}
    if ($AtomSymbol =~ /^F$/i) { $AtomicWeight = 18.998403; last WEIGHT;}
    if ($AtomSymbol =~ /^P$/i) { $AtomicWeight = 30.973761; last WEIGHT;}
    if ($AtomSymbol =~ /^S$/i) { $AtomicWeight = 32.065; last WEIGHT;}
    if ($AtomSymbol =~ /^Cl$/i) { $AtomicWeight = 35.453; last WEIGHT;}
    if ($AtomSymbol =~ /^Br$/i) { $AtomicWeight = 79.904; last WEIGHT;}
    if ($AtomSymbol =~ /^I$/i) { $AtomicWeight = 126.9045; last WEIGHT;}

    if ($AtomSymbol =~ /^OH$/i) { $AtomicWeight = 15.9994 + 1.00794; last WEIGHT;}
    if ($AtomSymbol =~ /^NH2$/i) { $AtomicWeight =  14.0067 + 2*1.00794; last WEIGHT;}
    if ($AtomSymbol =~ /^SH$/i) { $AtomicWeight = 32.065 + 1.00794; last WEIGHT;}
    if ($AtomSymbol =~ /^Me$/i) { $AtomicWeight = 12.0107 + 3*1.00794; last WEIGHT;}
    if ($AtomSymbol =~ /^Et$/i) { $AtomicWeight = 2*12.0107 + 5*1.00794; last WEIGHT;}
    if ($AtomSymbol =~ /^Pr$/i) { $AtomicWeight = 3*12.0107 + 7*1.00794; last WEIGHT;}
    if ($AtomSymbol =~ /^OMe$/i) { $AtomicWeight = 15.9994 + 12.0107 + 3*1.00794; last WEIGHT;}
    if ($AtomSymbol =~ /^OAc$/i) { $AtomicWeight = 2*15.9994 + 2*12.0107 + 3*1.00794; last WEIGHT;}
    if ($AtomSymbol =~ /^OOH$/i) { $AtomicWeight = 2*15.9994 + 1.00794; last WEIGHT;}
    if ($AtomSymbol =~ /^CN$/i) { $AtomicWeight = 12.0107 + 14.0067; last WEIGHT;}
    if ($AtomSymbol =~ /^NO2$/i) { $AtomicWeight = 14.0067 + 2*15.9994; last WEIGHT;}

    $AtomicWeight = undef;
  }
  return $AtomicWeight;
}

1;

__END__

=head1 NAME

ChainStr - Fatty Acyls (FA), Glycerolipids (GL) and Glycerophospholipids (GP) structure data generation methods

=head1 SYNOPSIS

use ChainStr;

use ChainStr qw(:all);

=head1 DESCRIPTION

ChainStr module provides these methods:

    AssignSubstituentStereoChemistry - Assign stereochemistry to substituents
    GenerateAtomBlockLines - Generate SD file atom data lines
    GenerateBondBlockLines - Generate SD file bond data lines
    GenerateChainStrData - Generate structure data for chains
    GenerateCmpdCountsLine - Generate SD file count data line
    IsAnySubstituentSpecifiedWithStereoChemistry - Check stereochemistry of
                                                   substituents
    SetupTemplateDataMap - Set up template data for a compound abbreviation

=head1 METHODS

=over 4

=item B<AssignSubstituentStereoChemistry>

    AssignSubstituentStereoChemistry($CmpdAbbrevTemplateDataMapRef,
        $CmdDataLinesRef);

Assign stereochemistry to substituents using structure and stereochemistry data available via
$CmpdAbbrevTemplateDataMapRef and $CmdDataLinesRef. And add new lines to existing
structure data using $CmdDataLinesRef.

=item B<GenerateAtomBlockLines>

    $AtomDataLines = GenerateAtomBlockLines($CmpdAbbrevTemplateDataMapRef,
        $Sn1AtomLinesRef, $Sn2AtomLinesRef, $Sn3AtomLinesRef);

Return atom data lines suitable for writing to SD file. Atom data for all approrpriate chains is
merged into a single string using new line character as delimiter.

=item B<GenerateBondBlockLines>

    $BondDataLines = GenerateBondBlockLines($CmpdAbbrevTemplateDataMapRef,
        $Sn1BondLinesRef, $Sn2AtomLinesRef, $Sn3AtomLinesRef);

Return bond data lines suitable for writing to SD file. Bond data for all approrpriate chains is
merged into a single string using new line character as delimiter.

=item B<GenerateChainStrData>

    ($AtomLinesArrayRef, $BondLinesArrayRef) = GenerateChainStrData($ChainType,
        $CmpdAbbrevTemplateDataMapRef);

Return array references containing atom and bond data lines for SD file. Appropriate atom
and bond data lines are generated using chain type and abbreviation template data.

=item B<GenerateCmpdCountsLine>

    $RetValue = GenerateCmpdCountsLine($CmpdAbbrevTemplateDataMapRef);

Return a formatted count data line for SD file.

=item B<IsAnySubstituentSpecifiedWithStereoChemistry>

    $Status = IsAnySubstituentSpecifiedWithStereoChemistry(
        $CmpdAbbrevTemplateDataMapRef);

Return 1 or 0 based on whether stereochemistry is specifed for any substituent.

=item B<SetupTemplateDataMap>

    SetupTemplateDataMap($TemplateType, $AbbrevTemplateDataMapRef,
        $TemplateData);

Setup compound abbreviation template data using a supported template.

=back

=head1 AUTHOR

Manish Sud

=head1 CONTRIBUTOR

Eoin Fahy

=head1 SEE ALSO

FAStr.pm, GLStr.pm, GPStr.pm, LMAPSStr.pm, SPStr.pm

=head1 COPYRIGHT

Copyright (C) 2006-2012. The Regents of the University of California. All Rights Reserved.

=head1 LICENSE

Modified BSD License

=cut
