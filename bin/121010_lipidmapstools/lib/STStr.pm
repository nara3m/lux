package STStr;
#
# File: STStr.pm
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
@EXPORT_OK = qw(ExpandSTCmpdAbbrevs GenerateCmpdOntologyData GenerateCmpdOntologySDDataLines GenerateSTStrData GenerateSDFile IsSTAbbrevSupported IsSTSubstituentsNameSupported IsSTDoubleBondsAbbrevOkay IsSTSubstituentsAbbrevOkay IsWildCardInSTAbbrev ParseSTAbbrev ParseSTSubstituentAbbrev ParseSTDoubleBondAbbrev SetupSTCmpdAbbrevTemplateDataMap ValidateSTAbbrev);
%EXPORT_TAGS = (
		all  => [@EXPORT, @EXPORT_OK],
		);

# Initialize ST data...
my(%STTemplatesDataMap, %STSupportedTemplateTypeMap);
_InitializeData();

#
# Parse template abbreviations with warnings for wild cards in templates...
#
sub ExpandSTCmpdAbbrevs {
  my($CmpdAbbrevsRef, $Abbrev, $ExpandedAbbrevRef, @ExpandedCmpdAbbrevs);

  ($CmpdAbbrevsRef) = @_;

  @ExpandedCmpdAbbrevs = ();

 ABBREV: for $Abbrev (@{$CmpdAbbrevsRef}) {
    if (!ValidateSTAbbrev($Abbrev)) {
      next ABBREV;
    }

    if (IsWildCardInSTAbbrev($Abbrev)) {
      warn "Warning: Ignoring compound abbreviation $Abbrev: Wild card is not supported.\n";
      next ABBREV;
    }

    if (!IsSTAbbrevSupported($Abbrev)) {
      warn "Warning: Ignoring compound abbreviation $Abbrev\n";
      next ABBREV;
    }

    push @ExpandedCmpdAbbrevs, $Abbrev;
  }

  return \@ExpandedCmpdAbbrevs;
}

# Generate ontology data lines for SD file...
sub GenerateCmpdOntologySDDataLines {
  my($CmpdAbbrevTemplateDataMapRef) = @_;
  my(@OntologyDataLines) = ();
  my($OntologyDataMapRef) = GenerateCmpdOntologyData($CmpdAbbrevTemplateDataMapRef);

  # Setup SD data lines...
  my($DataLabel, $DataValue);
  for $DataLabel (sort keys %{$OntologyDataMapRef}) {
    $DataValue = $OntologyDataMapRef->{$DataLabel};
    push @OntologyDataLines, ">  <$DataLabel>";
    push @OntologyDataLines, "$DataValue";
    push @OntologyDataLines, '';
  }
  return \@OntologyDataLines;
}

# Generate compound ontology data and return a reference to hash with property/value pairs...
sub GenerateCmpdOntologyData {
  my($CmpdAbbrevTemplateDataMapRef) = @_;

  my(%OntologyDataMap) = ();

  SetupAbbrevAndSystematicName($CmpdAbbrevTemplateDataMapRef, \%OntologyDataMap);
  SetupLMClassificationInfo($CmpdAbbrevTemplateDataMapRef, \%OntologyDataMap);
  SetupMultipleBondsCountInfo($CmpdAbbrevTemplateDataMapRef, \%OntologyDataMap);
  SetupChemicalGroupsCountInfo($CmpdAbbrevTemplateDataMapRef, \%OntologyDataMap);

  return \%OntologyDataMap;
}

# Multiple bond counts...
sub SetupMultipleBondsCountInfo {
  my($CmpdAbbrevTemplateDataMapRef, $OntologyDataMapRef) = @_;
  my($DoubleBondCount);

  $DoubleBondCount = keys %{$CmpdAbbrevTemplateDataMapRef->{DoubleBondsPosInfo}};
  $OntologyDataMapRef->{'Double Bonds'} = $DoubleBondCount;

}

# Miscellaneous chemical group count...
sub SetupChemicalGroupsCountInfo {
  my($CmpdAbbrevTemplateDataMapRef, $OntologyDataMapRef) = @_;

  # Count various substituents...
  my($SubstituentPosNum, $SubstituentsCount, $SubstituentAbbrev, $SubstituentName, $SubstituentIndex, $SubstituentAbbrevToNameRef);

  $SubstituentAbbrevToNameRef = ChainAbbrev::GetSubstituentsAbbrevToNameMap();

  for $SubstituentPosNum (sort {$a <=> $b} keys %{$CmpdAbbrevTemplateDataMapRef->{SubstituentsPosInfo}}) {
    $SubstituentsCount = @{$CmpdAbbrevTemplateDataMapRef->{SubstituentsPosInfo}{$SubstituentPosNum}{Abbrev}};
    for $SubstituentIndex (0 .. ($SubstituentsCount - 1)) {
      $SubstituentAbbrev = $CmpdAbbrevTemplateDataMapRef->{SubstituentsPosInfo}{$SubstituentPosNum}{Abbrev}[$SubstituentIndex];
      $SubstituentName = $SubstituentAbbrevToNameRef->{$SubstituentAbbrev};
      if ($SubstituentAbbrev !~ /^Hydro$/i) {
	if (exists $OntologyDataMapRef->{ucfirst($SubstituentName)}) {
	  $OntologyDataMapRef->{ucfirst($SubstituentName)} += 1;
	}
	else {
	  $OntologyDataMapRef->{ucfirst($SubstituentName)} = 1;
	}
      }
    }
  }
}

# Abbreviation and systematic name...
#
sub SetupAbbrevAndSystematicName {
  my($CmpdAbbrevTemplateDataMapRef, $OntologyDataMapRef) = @_;
  my($Abbrev, $SystematicName, $TemplateID);

  ($Abbrev, $SystematicName) = ('') x 2;
  $Abbrev = $CmpdAbbrevTemplateDataMapRef->{Abbrev};
  $TemplateID = $CmpdAbbrevTemplateDataMapRef->{TemplateID};

  if ($TemplateID =~ /^(CHOLESTANE|ERGOSTANE|CAMPESTANE|STIGMASTANE)$/i) {
    $SystematicName = _SetupSystematicName($TemplateID, $CmpdAbbrevTemplateDataMapRef);
  }

  $OntologyDataMapRef->{Abbrev} = $Abbrev;
  $OntologyDataMapRef->{'Systematic Name'} = $SystematicName;
}

#
# To do: Setup systematic name...
#
# Examples:
#   cholest-5-en-3b-ol
#   cholest-5-en-3b,4b-diol
#   5,6b-epoxy-5b-cholestan-3b-ol
#   cholest-4-en-3-one
#   cholest-3,5-dien-7-one
#   4a-methyl-cholest-7-en-3b-ol
#
sub _SetupSystematicName {
  my($TemplateID, $CmpdAbbrevTemplateDataMapRef) = @_;
  my($SystematicName, $SubstituentsAbbrev, $DoubleBondsAbbrev, $DoubleBondsName);

  $SubstituentsAbbrev = $CmpdAbbrevTemplateDataMapRef->{SubstituentsAbbrev};
  $DoubleBondsAbbrev = $CmpdAbbrevTemplateDataMapRef->{DoubleBondsAbbrev};

  $SystematicName = '';

  # Setup substituent names and name suffix...
  my($SubstituentsName, $NameSuffix) = ('') x 2;
  if (!IsEmpty($SubstituentsAbbrev)) {
    my($FirstPosNum, $SubstituentNameWithCount, $SubstituentAbbrev, $HydroxyFound, $SubstituentPosNumAndStereoPrefix, $SubstituentsNameDataMapRef);
    $SubstituentsNameDataMapRef = _SetupSubstituentsNameDataMap($CmpdAbbrevTemplateDataMapRef);

    # Setup the substituents name...
    $HydroxyFound = 0;
    for $FirstPosNum (sort {$a <=> $b} keys %{$SubstituentsNameDataMapRef->{FirstPosNum}}) {
      for $SubstituentNameWithCount (sort keys %{$SubstituentsNameDataMapRef->{NameWithCount}{$FirstPosNum}}) {
	$SubstituentAbbrev =  $SubstituentsNameDataMapRef->{Abbrev}{$FirstPosNum};
	$SubstituentPosNumAndStereoPrefix = $SubstituentsNameDataMapRef->{PosNumAndStereoPrefix}{$FirstPosNum}{$SubstituentNameWithCount};
	if ($SubstituentAbbrev =~ /^OH$/i) {
	  $NameSuffix = "${SubstituentPosNumAndStereoPrefix}-${SubstituentNameWithCount}";
	  $NameSuffix =~ s/hydroxy$/ol/g;
	  $HydroxyFound = 1;
	}
	elsif (!$HydroxyFound && $SubstituentAbbrev =~ /^Ke$/i) {
	  $NameSuffix = "${SubstituentPosNumAndStereoPrefix}-${SubstituentNameWithCount}";
	  $NameSuffix =~ s/oxo$/one/g;
	}
	else {
	  if ($SubstituentsName) {
	    $SubstituentsName .= "-${SubstituentPosNumAndStereoPrefix}";
	  }
	  else {
	    $SubstituentsName = "${SubstituentPosNumAndStereoPrefix}";
	  }
	  if (!IsEmpty($SubstituentNameWithCount)) {
	    $SubstituentsName .= "-${SubstituentNameWithCount}";
	  }
	}
      }
    }
  }

  # Double bond position and count suffix...
  $DoubleBondsName = '';
  if ($DoubleBondsAbbrev) {
    my($DoubleBondPos, $DoubleBondCount, $DoubleBondsPosInfo, $CountToNamePrefixMap, $CountPrefix);

    $DoubleBondsPosInfo = $CmpdAbbrevTemplateDataMapRef->{DoubleBondsPosInfo};
    $CountToNamePrefixMap = ChainAbbrev::GetCountToNamePrefixMap();

    $DoubleBondCount = 0;
    for $DoubleBondPos (sort {$a <=> $b} keys %{$DoubleBondsPosInfo}) {
      $DoubleBondCount++;
      if ($DoubleBondsName) {
	$DoubleBondsName .= ",${DoubleBondPos}";
      }
      else {
	$DoubleBondsName .= "${DoubleBondPos}";
      }
    }
    $CountPrefix = (exists $CountToNamePrefixMap->{$DoubleBondCount}) ? $CountToNamePrefixMap->{$DoubleBondCount} : '';
    $DoubleBondsName .= "-${CountPrefix}en";
  }

  # Setup name with any double bond info...
  my($NamePrefix);
  $NamePrefix = _SetupSystematicNamePrefix($TemplateID, $NameSuffix, $DoubleBondsName);

  # And the systematic name...
  if ($SubstituentsName) {
    $SystematicName = "${SubstituentsName}-${NamePrefix}";
  }
  else {
    $SystematicName = "${NamePrefix}";
  }
  if ($NameSuffix) {
    $SystematicName .= "-${NameSuffix}";
  }

  return $SystematicName;
}

# Setup systematic name prefix with any double bond info...
#
sub _SetupSystematicNamePrefix {
  my($TemplateID, $NameSuffix, $DoubleBondsName) = @_;
  my($NamePrefix, $Prefix);

  $NamePrefix = '';

  if ($TemplateID =~ /^CHOLESTANE$/i) {
    $Prefix = "cholest";
  }
  elsif ($TemplateID =~ /^ERGOSTANE$/i) {
    $Prefix = "ergost";
  }
  elsif ($TemplateID =~ /^CAMPESTANE$/i) {
    $Prefix = "campest";
  }
  elsif ($TemplateID =~ /^STIGMASTANE$/i) {
    $Prefix = "stigmast";
  }
  else {
    warn "Warning: Unknown sterol template: Systematic name prefix cann't be assigned\n";
  }

  if ($DoubleBondsName) {
    $NamePrefix = ($DoubleBondsName =~ /\,/) ? "${Prefix}a-${DoubleBondsName}" : "${Prefix}-${DoubleBondsName}";
  }
  else {
    $NamePrefix = $NameSuffix ? "${Prefix}an" : "${Prefix}ane";
  }

  return $NamePrefix;
}

# Setup substituents name...
sub _SetupSubstituentsNameDataMap {
  my($CmpdAbbrevTemplateDataMapRef) = @_;
  my($SubstituentPosNum, $SubstituentsCount, $SubstituentIndex, $SubstituentAbbrev, $SubstituentStereoChemistry, $SubstituentName, $SubtituentCount, $CountPrefix, $SubstituentPosPrefix, %SubstituentsDataMap, %FirstPosNumDataMap);

  # Setup a data map based on first position of the substituents and names along with count
  # prefix at that position for further sorting by names...
  %FirstPosNumDataMap = ();
  %{$FirstPosNumDataMap{FirstPosNum}} = ();
  %{$FirstPosNumDataMap{Abbrev}} = ();
  %{$FirstPosNumDataMap{NameWithCount}} = ();
  %{$FirstPosNumDataMap{PosNumAndStereoPrefix}} = ();

  if (!$CmpdAbbrevTemplateDataMapRef->{SubstituentsCount}) {
    return \%FirstPosNumDataMap;
  }

  # To allow sorting and counting at the current position, setup a hash using substituent abbreviation...
  %SubstituentsDataMap = ();
  @{$SubstituentsDataMap{Name}} = ();
  %{$SubstituentsDataMap{Abbrev}} = ();
  %{$SubstituentsDataMap{Count}} = ();
  %{$SubstituentsDataMap{NameWithCount}} = ();
  %{$SubstituentsDataMap{PosNumAndStereoPrefix}} = ();
  %{$SubstituentsDataMap{FirstPosNum}} = ();

  my($SubstituentAbbrevToNameMapRef) = ChainAbbrev::GetSubstituentsAbbrevToNameMap();

  for $SubstituentPosNum (sort {$a <=> $b} keys %{$CmpdAbbrevTemplateDataMapRef->{SubstituentsPosInfo}}) {
    $SubstituentsCount = @{$CmpdAbbrevTemplateDataMapRef->{SubstituentsPosInfo}{$SubstituentPosNum}{Spec}};

    for $SubstituentIndex (0 .. ($SubstituentsCount - 1)) {
      $SubstituentAbbrev = $CmpdAbbrevTemplateDataMapRef->{SubstituentsPosInfo}{$SubstituentPosNum}{Abbrev}[$SubstituentIndex];
      $SubstituentStereoChemistry = $CmpdAbbrevTemplateDataMapRef->{SubstituentsPosInfo}{$SubstituentPosNum}{StereoChemistry}[$SubstituentIndex];
      if (exists $SubstituentAbbrevToNameMapRef->{$SubstituentAbbrev}) {
	$SubstituentName = $SubstituentAbbrevToNameMapRef->{$SubstituentAbbrev};
      }
      else {
	if ($SubstituentAbbrev =~ /^H$/i) {
	  $SubstituentName = "Hydro";
	}
      }
      if (($SubstituentPosNum >= 17) && (!IsEmpty($SubstituentStereoChemistry))) {
	if ($SubstituentStereoChemistry =~ /^alpha$/i) {
	  # Assume S...
	  $SubstituentStereoChemistry = 'S';
	}
	elsif ($SubstituentStereoChemistry =~ /^beta$/i) {
	  # Assume R...
	  $SubstituentStereoChemistry = 'R';
	}
      }

      if (exists $SubstituentsDataMap{Abbrev}{$SubstituentName}) {
	$SubstituentsDataMap{Count}{$SubstituentName} += 1;
	$SubstituentsDataMap{PosNumAndStereoPrefix}{$SubstituentName} .= ",${SubstituentPosNum}${SubstituentStereoChemistry}";
      }
      else {
	push @{$SubstituentsDataMap{Name}}, $SubstituentName;
	$SubstituentsDataMap{Abbrev}{$SubstituentName} = $SubstituentAbbrev;
	$SubstituentsDataMap{Count}{$SubstituentName} = 1;
	if ($SubstituentAbbrev =~ /^Ep$/i) {
	  my($NextSubstituentPosNum);
	  $NextSubstituentPosNum = $SubstituentPosNum + 1;
	  $SubstituentsDataMap{PosNumAndStereoPrefix}{$SubstituentName} = "${SubstituentPosNum}${SubstituentStereoChemistry},${NextSubstituentPosNum}${SubstituentStereoChemistry}";
	}
	else {
	  $SubstituentsDataMap{PosNumAndStereoPrefix}{$SubstituentName} = "${SubstituentPosNum}${SubstituentStereoChemistry}";
	}
	$SubstituentsDataMap{FirstPosNum}{$SubstituentName} = $SubstituentPosNum;
      }
    }
  }
  # Setup name with count prefix...
  my($SubstituentNameWithCount, $CountToNamePrefixMapRef);
  $CountToNamePrefixMapRef = ChainAbbrev::GetCountToNamePrefixMap();

  for $SubstituentName (@{$SubstituentsDataMap{Name}}) {
    $SubtituentCount = $SubstituentsDataMap{Count}{$SubstituentName};
    $CountPrefix = $CountToNamePrefixMapRef->{$SubtituentCount};
    $SubstituentNameWithCount = "";
    if ($SubstituentName !~ /^Hydro$/i) {
      $SubstituentNameWithCount = "${CountPrefix}${SubstituentName}";
    }
    $SubstituentsDataMap{NameWithCount}{$SubstituentName} = $SubstituentNameWithCount;
  }

  # Setup a data map based on first position of the substituents and names along with count
  # prefix at that position for further sorting by names...
  my($FirstPosNum, $SubstituentPosNumAndStereoPrefix);
  %FirstPosNumDataMap = ();
  %{$FirstPosNumDataMap{FirstPosNum}} = ();
  %{$FirstPosNumDataMap{NameWithCount}} = ();
  %{$FirstPosNumDataMap{PosNumAndStereoPrefix}} = ();

  for $SubstituentName (@{$SubstituentsDataMap{Name}}) {
    $FirstPosNum = $SubstituentsDataMap{FirstPosNum}{$SubstituentName};
    $SubstituentAbbrev = $SubstituentsDataMap{Abbrev}{$SubstituentName};
    $SubstituentNameWithCount = $SubstituentsDataMap{NameWithCount}{$SubstituentName};
    $SubstituentPosNumAndStereoPrefix = $SubstituentsDataMap{PosNumAndStereoPrefix}{$SubstituentName};

    if (exists $FirstPosNumDataMap{FirstPosNum}{$FirstPosNum}) {
      $FirstPosNumDataMap{NameWithCount}{$FirstPosNum}{$SubstituentNameWithCount} = $SubstituentNameWithCount;
      $FirstPosNumDataMap{PosNumAndStereoPrefix}{$FirstPosNum}{$SubstituentNameWithCount} = $SubstituentPosNumAndStereoPrefix;
    }
    else {
      $FirstPosNumDataMap{FirstPosNum}{$FirstPosNum} = $FirstPosNum;
      $FirstPosNumDataMap{Abbrev}{$FirstPosNum} = $SubstituentAbbrev;
      %{$FirstPosNumDataMap{NameWithCount}{$FirstPosNum}} = ();
      $FirstPosNumDataMap{NameWithCount}{$FirstPosNum}{$SubstituentNameWithCount} = $SubstituentNameWithCount;
      %{$FirstPosNumDataMap{PosNumAndStereoPrefix}{$FirstPosNum}} = ();
      $FirstPosNumDataMap{PosNumAndStereoPrefix}{$FirstPosNum}{$SubstituentNameWithCount} = $SubstituentPosNumAndStereoPrefix;
    }
  }

  return \%FirstPosNumDataMap;
}

# LM classification info...
sub SetupLMClassificationInfo {
  my($CmpdAbbrevTemplateDataMapRef, $OntologyDataMapRef) = @_;
  my($LMCategory, $LMMainClass, $LMSubClass, $STAbbrev, $TemplateID,);

  $LMCategory = $CmpdAbbrevTemplateDataMapRef->{LMCategory};
  $LMMainClass = $CmpdAbbrevTemplateDataMapRef->{LMMainClass};
  $LMSubClass = $CmpdAbbrevTemplateDataMapRef->{LMSubClass};

  $OntologyDataMapRef->{'LM Category'} = $LMCategory;
  $OntologyDataMapRef->{'LM Main Class'} = $LMCategory . $LMMainClass;
  $OntologyDataMapRef->{'LM Sub Class'} = $LMCategory . $LMMainClass . $LMSubClass;

}


# Generate a SD file...
sub GenerateSDFile {
  my($SDFileName, $CmdAbbrevsRef, $Abbrev, $CmpdDataString, $CmpdAbbrevTemplateDataMapRef);

  ($SDFileName, $CmdAbbrevsRef) = @_;

  open SDFILE, ">$SDFileName" or die "Error: Couldn't open $SDFileName: $! \n";

  ABBREV: for $Abbrev (@$CmdAbbrevsRef) {

    # Set up a compound data string...
    $CmpdDataString = '';
    $CmpdDataString = _GenerateCmdDataString($Abbrev);

    # Write it out to the file...
    print SDFILE "$CmpdDataString\n";
  }
  close SDFILE;
}

# Generate atom and bond block lines for template and substituents...
sub GenerateSTStrData {
  my($CmpdAbbrevTemplateDataMapRef) = @_;
  my($Index, $SubstituentAtomNum, $SubstituentsPosInfo, $SubstituentPos, $SubstituentCount, $SustituentSpec, $SustituentAbbrev, $SustituentStereoChemistry, $CurrentAtomNum, $CurrentAtomX, $CurrentAtomY, $SubstituentAtomXOffset, $SubstituentAtomYOffset, $SubstituentAtomX, $SubstituentAtomY, @AtomBlockLines, @BondBlockLines );

  @AtomBlockLines = ();
  @BondBlockLines = ();

  # Generate structure data for sustituents...
  $SubstituentAtomNum = $CmpdAbbrevTemplateDataMapRef->{AtomCount};

  $SubstituentsPosInfo = $CmpdAbbrevTemplateDataMapRef->{SubstituentsPosInfo};
  for $SubstituentPos (sort {$a <=> $b} keys %{$SubstituentsPosInfo}) {
    $CurrentAtomNum = $SubstituentPos;
    ($CurrentAtomX, $CurrentAtomY) = _GetAtomPositions($CurrentAtomNum, $CmpdAbbrevTemplateDataMapRef);
    $SubstituentCount = $SubstituentsPosInfo->{$SubstituentPos}{Count};
    for $Index (0 ... ($SubstituentCount - 1)) {
      $SubstituentAtomNum++;
      $SustituentSpec = $SubstituentsPosInfo->{$SubstituentPos}{Spec}[$Index];
      $SustituentAbbrev = $SubstituentsPosInfo->{$SubstituentPos}{Abbrev}[$Index];
      $SustituentStereoChemistry = $SubstituentsPosInfo->{$SubstituentPos}{StereoChemistry}[$Index];

      ($SubstituentAtomXOffset, $SubstituentAtomYOffset) = _GetSubstituentAtomOffset($CurrentAtomNum, $SustituentAbbrev, $CmpdAbbrevTemplateDataMapRef);

      # Set up final value...
      $SubstituentAtomX = $CurrentAtomX + $SubstituentAtomXOffset;
      $SubstituentAtomY = $CurrentAtomY + $SubstituentAtomYOffset;

      # Setup atom data line...
      my($AtomSpec) = $SustituentAbbrev;
      if ($SustituentAbbrev =~ /^Ke$/i) {
	$AtomSpec = 'O';
      }
      elsif ($SustituentAbbrev =~ /^My$/i) {
	$AtomSpec = 'C';
      }
      elsif ($SustituentAbbrev =~ /^Me$/i) {
	$AtomSpec = 'C';
      }
      elsif ($SustituentAbbrev =~ /^Ep$/i) {
	$AtomSpec = 'O';
      }
      elsif ($SustituentAbbrev =~ /^OH$/i) {
	$AtomSpec = 'O';
      }
      elsif ($SustituentAbbrev =~ /^NH2$/i) {
	$AtomSpec = 'N';
      }
      elsif ($SustituentAbbrev =~ /^SH$/i) {
	$AtomSpec = 'S';
      }

      my($AtomLine) = LMAPSStr::GenerateCmpdAtomLine($SubstituentAtomX, $SubstituentAtomY, "0.0000", $AtomSpec);
      push @AtomBlockLines, $AtomLine;

      # Setup bond line along with any stereochemistry spec...
      my($BondStereo) = ($SustituentStereoChemistry =~ /^(alpha)$/i) ? 6 : (($SustituentStereoChemistry =~ /^(beta)$/i) ? 1 : '');
      my($BondOrder) = ChainAbbrev::GetSubstituentBondOrder($SustituentAbbrev);
      my($BondLine);
      if (IsEmpty($BondStereo)) {
	$BondLine = LMAPSStr::GenerateCmpdBondLine($CurrentAtomNum, $SubstituentAtomNum, $BondOrder);
      }
      else {
	$BondLine = LMAPSStr::GenerateCmpdBondLine($CurrentAtomNum, $SubstituentAtomNum, $BondOrder, $BondStereo);
      }
      push @BondBlockLines, $BondLine;

      # Add another bond line for epoxy with next atom...
      if ($SustituentAbbrev =~ /^Ep$/i) {
	my($SecondBondLine);
	if (IsEmpty($BondStereo)) {
	  ($SecondBondLine) = LMAPSStr::GenerateCmpdBondLine(($CurrentAtomNum + 1), $SubstituentAtomNum, $BondOrder);
	}
	else {
	  ($SecondBondLine) = LMAPSStr::GenerateCmpdBondLine(($CurrentAtomNum + 1), $SubstituentAtomNum, $BondOrder, $BondStereo);
	}
	push @BondBlockLines, $SecondBondLine;
      }
    }
  }

  return (\@AtomBlockLines, \@BondBlockLines);
}

# Parse abbreviation...
sub ParseSTAbbrev {
  my($Abbrev) = @_;
  my($STType,  $SubstituentsAbbrev, $DoubleBondsAbbrev);

  $DoubleBondsAbbrev = '';
  if ($Abbrev =~ /\:/) {
    ($STType,  $SubstituentsAbbrev, $DoubleBondsAbbrev) = $Abbrev =~ /^(.*?)\((.*?)\:(.*?)\)$/;
  }
  else {
    ($STType,  $SubstituentsAbbrev) = $Abbrev =~ /^(.*?)\((.*?)\)$/;
    $DoubleBondsAbbrev = '';
    if (IsEmpty($SubstituentsAbbrev)) {
      if ($STType =~ /^(CHOLESTANE|ERGOSTANE|CAMPESTANE|STIGMASTANE)$/i) {
	$SubstituentsAbbrev = "5,a,H";
      }
    }
  }
  return ($STType,  $SubstituentsAbbrev, $DoubleBondsAbbrev);
}

# Does template exist to handle this abbreviation?
#
sub IsSTAbbrevSupported {
  my($STAbbrev) = @_;
  my($STType, $SubstituentsAbbrev, $DoubleBondsAbbrev);

  ($STType, $SubstituentsAbbrev, $DoubleBondsAbbrev) = ParseSTAbbrev($STAbbrev);

  if ($SubstituentsAbbrev) {
    if (!IsSTSubstituentsAbbrevOkay($STAbbrev, $STType, $SubstituentsAbbrev, $DoubleBondsAbbrev)) {
      return 0;
    }
  }

  if ($DoubleBondsAbbrev) {
    if (!IsSTDoubleBondsAbbrevOkay($STAbbrev, $STType, $SubstituentsAbbrev, $DoubleBondsAbbrev)) {
      return 0;
    }
  }
  return 1;
}

# Does abbreviation contain any wild cards?
sub IsWildCardInSTAbbrev {
  my($Abbrev) = @_;

  my($Status, $STType, $SubstituentsAbbrev, $DoubleBondsAbbrev);

  ($STType, $SubstituentsAbbrev, $DoubleBondsAbbrev) = ParseSTAbbrev($Abbrev);

  $Status = (($STType =~ /\*/) || ($SubstituentsAbbrev =~ /\*/) || ($DoubleBondsAbbrev =~ /\*/)) ? 1 : 0;
  return $Status;
}

# Make sure it's a valid substituent abbreviation...
#
#
sub IsSTSubstituentsAbbrevOkay {
  my($STAbbrev, $STType, $SubstituentsAbbrev, $DoubleBondsAbbrev) = @_;

  if ($SubstituentsAbbrev =~ /\*/) {
    warn "Warning: Ignoring abbreviation $STAbbrev : Wild card found in substituents abbreviation\n";
    return 0;
  }
  my($Substituent, $SubstituentPos, $SubstituentAbbrev, $SubstituentStereoChemistry, @SubstituentsList, %SubstituentsDataMap);
  %SubstituentsDataMap = ();
  %{$SubstituentsDataMap{SubstituentsPos}} = ();
  %{$SubstituentsDataMap{SubstituentsAtPosCount}} = ();
  %{$SubstituentsDataMap{SubstituentsAtPos}} = ();
  %{$SubstituentsDataMap{SubstituentsAtPosStereoSpecCount}} = ();
  %{$SubstituentsDataMap{SubstituentsAtPosStereoSpec}} = ();
  %{$SubstituentsDataMap{SubstituentsAtPosBondOrderCount}} = ();

  (@SubstituentsList) = split /\//, $SubstituentsAbbrev;
  SUBSTITUENT: for $Substituent (@SubstituentsList) {
    ($SubstituentPos, $SubstituentStereoChemistry, $SubstituentAbbrev) = ParseSTSubstituentAbbrev($Substituent);

    # Validate position value...
    if (IsEmpty($SubstituentPos)) {
      warn "Warning: Ignoring abbreviation $STAbbrev : Missing substituent position value\n";
      return 0;
    }
    if ($SubstituentPos =~ /[^0-9]+/) {
      warn "Warning: Ignoring abbreviation $STAbbrev : Substituent position value must be an integer\n";
      return 0;
    }
    # Check substituent name...
    if (!IsSTSubstituentsNameSupported($SubstituentAbbrev)) {
      warn "Warning: Ignoring abbreviation $STAbbrev : Unknown substituent $SubstituentAbbrev\n";
      return 0;
    }
    if ($STType =~ /^CHOLESTANE$/i) {
      if ($SubstituentPos !~ /^(1|2|3|4|5|6|7|11|12|15|16|22|23|24|25|26|27)$/i) {
	warn "Warning: Ignoring abbreviation $STAbbrev : Substituent position, $SubstituentPos, is not valid. Supported values: 1-7, 11, 12, 15, 16, 22-27\n";
      return 0;
      }
    }
    elsif ($STType =~ /^(ERGOSTANE|CAMPESTANE|STIGMASTANE)$/i) {
      if ($SubstituentPos !~ /^(1|2|3|4|5|6|7|11|12|15|16|22|23|25|26|27)$/i) {
	warn "Warning: Ignoring abbreviation $STAbbrev : Substituent position, $SubstituentPos, is not valid. Supported values: 1-7, 11, 12, 15, 16, 22, 23, 25-27\n";
      return 0;
      }
    }
    # Allow only a/b or alpha/beta for SubstituentStereoChemistry
    if (!IsEmpty($SubstituentStereoChemistry)) {
      if ($SubstituentStereoChemistry !~ /^(a|alpha|b|beta)$/i) {
	warn "Warning: Ignoring abbreviation $STAbbrev: Unknown stereochemistry specification $SubstituentStereoChemistry. Supported values: a, alpha, b, or beta\n";
	return 0;
      }
    }
    if (exists $SubstituentsDataMap{SubstituentsPos}{$SubstituentPos}) {
      $SubstituentsDataMap{SubstituentsAtPos}{$SubstituentPos} .= ",${SubstituentAbbrev}";
      $SubstituentsDataMap{SubstituentsAtPosCount}{$SubstituentPos} += 1;
      $SubstituentsDataMap{SubstituentsAtPosStereoSpec}{$SubstituentPos} .= ",${SubstituentStereoChemistry}";
      if (!IsEmpty($SubstituentStereoChemistry)) {
	$SubstituentsDataMap{SubstituentsAtPosStereoSpecCount}{$SubstituentPos} += 1;
      }
      $SubstituentsDataMap{SubstituentsAtPosBondOrderCount}{$SubstituentPos} += GetSTSubstituentBondOrder($SubstituentAbbrev);
    }
    else {
      $SubstituentsDataMap{SubstituentsPos}{$SubstituentPos} = $SubstituentPos;
      $SubstituentsDataMap{SubstituentsAtPos}{$SubstituentPos} = $SubstituentAbbrev;
      $SubstituentsDataMap{SubstituentsAtPosCount}{$SubstituentPos} = 1;
      $SubstituentsDataMap{SubstituentsAtPosStereoSpec}{$SubstituentPos} = $SubstituentStereoChemistry;
      $SubstituentsDataMap{SubstituentsAtPosStereoSpecCount}{$SubstituentPos} = IsEmpty($SubstituentStereoChemistry) ? 0 : 1;
      $SubstituentsDataMap{SubstituentsAtPosBondOrderCount}{$SubstituentPos} = GetSTSubstituentBondOrder($SubstituentAbbrev);
    }
    #  Maximum number of substituents allowed at specific positions...
    if ($SubstituentsDataMap{SubstituentsAtPosCount}{$SubstituentPos} > 1) {
      warn "Warning: Ignoring abbreviation $STAbbrev : More than one substituent found at position $SubstituentPos\n";
      return 0;
    }
  }
  return 1;
}

# Make sure it's a valid double bond abbreviation...
#
sub IsSTDoubleBondsAbbrevOkay {
  my($STAbbrev, $STType, $SubstituentsAbbrev, $DoubleBondsAbbrev) = @_;
  my($DoubleBondPos1, $DoubleBondPos2, $DoubleBond, @DoubleBondsList, %DoubleBondsPosInfo);

  if ($DoubleBondsAbbrev =~ /\*/) {
    warn "Warning: Ignoring abbreviation $STAbbrev : Wild card found in double bond abbreviation\n";
    return 0;
  }
  (@DoubleBondsList) = split /\//, $DoubleBondsAbbrev;
  for $DoubleBond (@DoubleBondsList) {
    $DoubleBondPos1 = ''; $DoubleBondPos2 = '';
    ($DoubleBondPos1, $DoubleBondPos2) = ParseSTDoubleBondAbbrev($DoubleBond);
    if ($DoubleBondPos1 =~ /[^0-9]/) {
      warn "Warning: Ignoring abbreviation $STAbbrev : Invalid charactar in double bond, $DoubleBond, specification\n";
      return 0;
    }
    if (!IsEmpty($DoubleBondPos2)) {
      if ($DoubleBondPos2 =~ /[^0-9]/) {
	warn "Warning: Ignoring abbreviation $STAbbrev : Invalid charactar in double bond, $DoubleBond, specification\n";
	return 0;
      }
    }
    if (IsEmpty($DoubleBondPos2)) {
      $DoubleBondPos2 = $DoubleBondPos1 + 1;
    }
    else {
      if (($DoubleBondPos1 + 1) != $DoubleBondPos2) {
	warn "Warning: Ignoring abbreviation $STAbbrev : Invalid double bond, $DoubleBond, specification: Second atom position, $DoubleBondPos2, is not valid.\n";
	return 0;
      }
    }
    if (exists $DoubleBondsPosInfo{$DoubleBondPos1}) {
      warn "Warning: Ignoring abbreviation $STAbbrev : Invalid double bond, $DoubleBond, specification: Atom position, $DoubleBondPos1, is already involved in a double bond.\n";
      return 0;
    }
    if ($STType =~ /^CHOLESTANE$/i) {
      if ($DoubleBondPos1 !~ /^(1|2|3|4|5|6|11|15|22|23|24|25)$/i) {
	warn "Warning: Ignoring abbreviation $STAbbrev : Invalid double bond, $DoubleBond, specification: First atom position, $DoubleBondPos1, is not valid. Supported values: 1-6, 11, 15, 16, 22-25\n";
      return 0;
      }
    }
    elsif ($STType =~ /^(ERGOSTANE|CAMPESTANE|STIGMASTANE)$/i) {
      if ($DoubleBondPos1 !~ /^(1|2|3|4|5|6|11|15|22|25)$/i) {
	warn "Warning: Ignoring abbreviation $STAbbrev : Invalid double bond, $DoubleBond, specification: First atom position, $DoubleBondPos1, is not valid. Supported values: 1-6, 11, 15, 16, 22, 25\n";
      return 0;
      }
    }
    %{$DoubleBondsPosInfo{$DoubleBondPos1}} = ();
    $DoubleBondsPosInfo{$DoubleBondPos1}{$DoubleBondPos2} = "${DoubleBondPos1}(${DoubleBondPos1})";
  }

  return 1;
}

# Is it a supported ST substituent name?
sub IsSTSubstituentsNameSupported {
  my($SubstituentAbbrev) = @_;

  my($SubstituentsAbbrevToNameMapRef) = ChainAbbrev::GetSubstituentsAbbrevToNameMap();

  if (!exists $SubstituentsAbbrevToNameMapRef->{$SubstituentAbbrev}) {
    if ($SubstituentAbbrev !~ /^[H]$/i) {
      return 0;
    }
  }
  if ($SubstituentAbbrev =~ /^(COOH|CHO)$/i) {
    return 0;
  }
  return 1;
}

# Get the bond order for substituent...
sub GetSTSubstituentBondOrder {
  my($SubstituentAbbrev) = @_;
  my($BondOrder, $SubstituentsAbbrevToNameMapRef);

  $BondOrder = 1;
  $SubstituentsAbbrevToNameMapRef = ChainAbbrev::GetSubstituentsAbbrevToNameMap();
  if (exists $SubstituentsAbbrevToNameMapRef->{$SubstituentAbbrev}) {
    $BondOrder = ChainAbbrev::GetSubstituentBondOrder($SubstituentAbbrev)
  }

  return $BondOrder;
}

# Parse substituent abbreviation...
#
# Format:
#   . SubstituentPos,[SubstituentSterochemistry],Substituent/SubstituentPos,[SubstituentSterochemistry],Substituent/...
#
# Caveats:
#  . SubstituentSterochemistry specification is optional. Default value: none. Possible values: a/b or alpha/beta.
#    alpha: Going down the plane (6); beta: Going up the plan (1).
#
# Examples:
#   . 3,a,OH/5,b,H
#   . 3,,O/5,b,H is equivalent to 3,,k/5,b,H
#   . 3,,K/5,a,H
#
#   . Cholesterol: CHOLESTANE(3,b,OH/8,b,H/9,a,H/14,a,H:5)
#
sub ParseSTSubstituentAbbrev {
  my($Abbrev) = @_;
  my($SubstituentPos, $SubstituentAbbrev, $SubstituentStereoChemistry, @AbbrevList);

  ($SubstituentPos, $SubstituentAbbrev, $SubstituentStereoChemistry) = ('') x 3;
  @AbbrevList = split /\,/, $Abbrev;
  if (@AbbrevList == 3) {
    ($SubstituentPos, $SubstituentAbbrev, $SubstituentStereoChemistry) = @AbbrevList;
  }

  return ($SubstituentPos, $SubstituentAbbrev, $SubstituentStereoChemistry);
}

# Format:
#   . AtomPos[(NextAtomPos)]
#
# Caveats:
#  . NextAtomPos involved in double bond is optional. Default value: AtomPos + 1
#
# Examples:
#   . 26
#   . 9(11)
sub ParseSTDoubleBondAbbrev {
  my($Abbrev) = @_;
  my($Pos1, $Pos2) = ('') x 2;

  if ($Abbrev =~ /\(/ && $Abbrev =~ /\)/) {
    ($Pos1, $Pos2) = $Abbrev =~ /^(.*?)\((.*?)\)$/;
    if (!defined $Pos1) {
      $Pos1 = '';
    }
    if (!defined $Pos2) {
      $Pos2 = '';
    }
  }
  else {
    $Pos1 = $Abbrev;
  }

  return ($Pos1, $Pos2);
}

#
# Check out the validity of ST abbreviation...
sub ValidateSTAbbrev {
  my($Abbrev) = @_;

  if (!($Abbrev =~ /\(/ && $Abbrev =~ /\)/)) {
    warn "Warning: Ignored ST compound abbreviation $Abbrev due to incorrect format. Valid format: SterolType(SubstituentsAbbrev:DoubleBondsAbbrev)\n";
    return 0;
  }

  # Make sure sterol type is these...
  if ($Abbrev =~ /^\(/) {
    warn "Warning: Ignored ST compound abbreviation $Abbrev due to incorrect format. Valid format: SterolType(SubstituentsAbbrev:DoubleBondsAbbrev)\n";
    return 0;
  }
  # Number of '(' and ')' parenthesis must match...
  my(@LeftParenthesis) = $Abbrev =~ /\(/g;
  my(@RightParenthesis) = $Abbrev =~ /\)/g;
  if ($#LeftParenthesis != $#RightParenthesis) {
    warn "Warning: Ignored ST compound abbreviation $Abbrev due to incorrect format: Number of ( and ) doesn't match. Valid format: SterolType(SubstituentsAbbrev:DoubleBondsAbbrev).\n";
    return 0;
  }

  my($STType, $SubstituentsAbbrev, $DoubleBondsAbbrev);
  ($STType, $SubstituentsAbbrev, $DoubleBondsAbbrev) = ParseSTAbbrev($Abbrev);

  if (!defined $STType) {
    warn "Warning: Ignored ST compound abbreviation $Abbrev due to incorrect format: SterolType is not specifed. Valid format: SterolType(SubstituentsAbbrev:DoubleBondsAbbrev).\n";
  }
  if (IsEmpty($SubstituentsAbbrev)) {
    warn "Warning: Ignored ST compound abbreviation $Abbrev due to incorrect format: SubstituentsAbbrev is not specifed. Valid format: SterolType(SubstituentsAbbrev:DoubleBondsAbbrev).\n";
  }

  if ($STType !~ /\*/ ) {
    if (!(exists $STSupportedTemplateTypeMap{$STType})) {
      warn "Warning: Ignored St compound abbreviation $Abbrev) : Unknown sterol template type $STType.\n";
      return 0;
    }
  }
  return 1;
}

# Generate ST compound count line...
sub GenerateSTCmpdCountsLine {
  my($CmpdAbbrevTemplateDataMapRef) = @_;
  my($CmpdAtomCount, $CmpdBondCount);

  # Template atom counts...
  $CmpdAtomCount = $CmpdAbbrevTemplateDataMapRef->{AtomCount};
  $CmpdBondCount = $CmpdAbbrevTemplateDataMapRef->{BondCount};

  $CmpdAtomCount += $CmpdAbbrevTemplateDataMapRef->{SubstituentsCount};
  $CmpdBondCount += $CmpdAbbrevTemplateDataMapRef->{SubstituentsBondCount};

  return LMAPSStr::GenerateCmpdCountsLine($CmpdAtomCount, $CmpdBondCount);
}

# Generate ST atom block lines...
sub GenerateSTAtomBlockLines {
  my($CmpdAbbrevTemplateDataMapRef, $AtomBlockLinesRef) = @_;
  my($StrDataLines);

  $StrDataLines = '';
  # Template atom block...
  $StrDataLines .= JoinWords(\@{$CmpdAbbrevTemplateDataMapRef->{AtomBlockLines}}, "\n", 0);

  $StrDataLines .= "\n" . JoinWords($AtomBlockLinesRef, "\n", 0);

  return $StrDataLines;
}

# Generate bondatom block lines...
sub GenerateSTBondBlockLines {
  my($CmpdAbbrevTemplateDataMapRef, $BondBlockLinesRef) = @_;
  my($StrDataLines);

  $StrDataLines = '';

  # Template bond block..
  $StrDataLines .= JoinWords(\@{$CmpdAbbrevTemplateDataMapRef->{BondBlockLines}}, "\n", 0);

  $StrDataLines .= "\n" . JoinWords($BondBlockLinesRef, "\n", 0);

  return $StrDataLines;
}

# Set up template data for a specific cmpd abbreviation and return a
# reference to the hash...
sub SetupSTCmpdAbbrevTemplateDataMap {
  my($Abbrev, $STType, $SubstituentsAbbrev, $DoubleBondsAbbrev, $TemplateID, %AbbrevTemplateDataMap);

  ($Abbrev) = @_;

  %AbbrevTemplateDataMap = ();

  ($STType, $SubstituentsAbbrev, $DoubleBondsAbbrev) = ParseSTAbbrev($Abbrev);
  $TemplateID = _GetTemplateID($Abbrev, $STType);

  $AbbrevTemplateDataMap{TemplateID} = $STType;
  $AbbrevTemplateDataMap{Abbrev} = $Abbrev;
  $AbbrevTemplateDataMap{SubstituentsAbbrev} = $SubstituentsAbbrev;
  $AbbrevTemplateDataMap{DoubleBondsAbbrev} = $DoubleBondsAbbrev;

  _SetupTemplateDataMap(\%AbbrevTemplateDataMap, $STTemplatesDataMap{$TemplateID});
  return \%AbbrevTemplateDataMap;
}

# Get atom positions...
sub _GetAtomPositions {
  my($AtomNum, $CmpdAbbrevTemplateDataMapRef) = @_;
  my($AtomX, $AtomY, $AtomLine);

  ($AtomX, $AtomY) = ('') x 2;
  if ($AtomNum <= $CmpdAbbrevTemplateDataMapRef->{AtomCount}) {
    $AtomLine = $CmpdAbbrevTemplateDataMapRef->{AtomBlockLines}[$AtomNum - 1];
    ($AtomX, $AtomY) = LMAPSStr::ParseCmpdAtomLine($AtomLine);
  }

  return ($AtomX, $AtomY);
}

# Setup atom positions for substituent...
sub _GetSubstituentAtomOffset {
  my($AtomNum, $SubstituentAbbrev, $CmpdAbbrevTemplateDataMapRef) = @_;
  my($AtomXOffset, $AtomYOffset, $XOffset, $YOffset, $NSYOffset);

  $XOffset = 0.6266; $YOffset = 0.3618; $NSYOffset = 0.6266;

  # Set up final offset based on which direction it needs to go: N, NE, E, SE, S, SW, W, NW
  if ($SubstituentAbbrev =~ /^Ep$/i) {
    ATOMPOS: {
      if ($AtomNum =~ /^(8|10|11|16)$/) { $AtomXOffset = 0;  $AtomYOffset = $NSYOffset; last ATOMPOS; } # N
      if ($AtomNum =~ /^(7|12|15|22|24)$/) { $AtomXOffset = $XOffset;  $AtomYOffset = $YOffset; last ATOMPOS; } # NE
      if ($AtomNum =~ /^(4|6|13|14|23|25)$/) { $AtomXOffset = $XOffset;  $AtomYOffset = -$YOffset; last ATOMPOS; } # SE
      if ($AtomNum =~ /^(3|5)$/) { $AtomXOffset = 0;  $AtomYOffset = -$NSYOffset; last ATOMPOS; } # S
      if ($AtomNum =~ /^(2)$/) { $AtomXOffset = -$XOffset;  $AtomYOffset = -$YOffset; last ATOMPOS; } # SW
      if ($AtomNum =~ /^(1|9|17)$/) { $AtomXOffset = -$XOffset;  $AtomYOffset = $YOffset; last ATOMPOS; } # NW
      $AtomXOffset = -$XOffset;  $AtomYOffset = $YOffset;
    }
  }
  else {
    ATOMPOS: {
      if ($AtomNum =~ /^(1|10|12|13|22|24)$/) { $AtomXOffset = 0;  $AtomYOffset = $NSYOffset; last ATOMPOS; } # N
      if ($AtomNum =~ /^(16|17)$/) { $AtomXOffset = $XOffset;  $AtomYOffset = $YOffset; last ATOMPOS; } # NE
      if ($AtomNum =~ /^(7|8|15|25|26|27)$/) { $AtomXOffset = $XOffset;  $AtomYOffset = -$YOffset; last ATOMPOS; } # SE
      if ($AtomNum =~ /^(4|5|6|9|14|23)$/) { $AtomXOffset = 0;  $AtomYOffset = -$NSYOffset; last ATOMPOS; } # S
      if ($AtomNum =~ /^(3)$/) { $AtomXOffset = -$XOffset;  $AtomYOffset = -$YOffset; last ATOMPOS; } # SW
      if ($AtomNum =~ /^(2|11)$/) { $AtomXOffset = -$XOffset;  $AtomYOffset = $YOffset; last ATOMPOS; } # NW
      $AtomXOffset = -$XOffset;  $AtomYOffset = $YOffset;
    }
  }
  return ($AtomXOffset, $AtomYOffset);
}

# Set up template data for ST...
sub _SetupTemplateDataMap {
  my($AbbrevTemplateDataMapRef, $TemplateData) = @_;
  my($Index, $AbbrevID, $NumOfCarbonAtoms, $RingJunction1AtomNum1, $RingJunction1AtomNum2, $RingJunction2AtomNum1, $RingJunction2AtomNum2, $RingJunction3AtomNum1, $RingJunction3AtomNum2, $LMCategory, $LMMainClass, $LMSubClass, $CmpdString);

  ($AbbrevID, $NumOfCarbonAtoms, $RingJunction1AtomNum1, $RingJunction1AtomNum2, $RingJunction2AtomNum1, $RingJunction2AtomNum2, $RingJunction3AtomNum1, $RingJunction3AtomNum2, $LMCategory, $LMMainClass, $LMSubClass, $CmpdString) = split /\|/, $TemplateData;

  $AbbrevTemplateDataMapRef->{AbbrevID} = $AbbrevID;
  $AbbrevTemplateDataMapRef->{NumOfCarbonAtoms} = $NumOfCarbonAtoms;

  $AbbrevTemplateDataMapRef->{LMCategory} = $LMCategory;
  $AbbrevTemplateDataMapRef->{LMMainClass} = $LMMainClass;
  $AbbrevTemplateDataMapRef->{LMSubClass} = $LMSubClass;

  # Setup ring junction atom numbers...
  my($RingJunctionAtom, @RingJunctionAtomsList, %RingJunctionAtoms);
  @RingJunctionAtomsList = ();
  %RingJunctionAtoms = ();
  push @RingJunctionAtomsList, ($RingJunction1AtomNum1, $RingJunction1AtomNum2, $RingJunction2AtomNum1, $RingJunction2AtomNum2, $RingJunction3AtomNum1, $RingJunction3AtomNum2);
  RINGJUNCTIONATOM: for $RingJunctionAtom (@RingJunctionAtomsList) {
    if (IsEmpty($RingJunctionAtom)) {
      next RINGJUNCTIONATOM;
    }
    $RingJunctionAtoms{$RingJunctionAtom} = $RingJunctionAtom;
  }
  $AbbrevTemplateDataMapRef->{RingJunctionAtoms} = \%RingJunctionAtoms;

  my($Abbrev, $Substituents, $DoubleBonds);
  $Abbrev = $AbbrevTemplateDataMapRef->{Abbrev};
  $Substituents = $AbbrevTemplateDataMapRef->{SubstituentsAbbrev};
  $DoubleBonds = $AbbrevTemplateDataMapRef->{DoubleBondsAbbrev};

  # Setup double bonds information...
  my($DoubleBond, $DoubleBondPos1, $DoubleBondPos2, @DoubleBondsList, %DoubleBondsPosInfo);
  %DoubleBondsPosInfo = (); @DoubleBondsList = ();
  (@DoubleBondsList) = split /\//, $DoubleBonds;
  for $DoubleBond (@DoubleBondsList) {
    $DoubleBondPos1 = ''; $DoubleBondPos2 = '';
    ($DoubleBondPos1, $DoubleBondPos2) = ParseSTDoubleBondAbbrev($DoubleBond);
    if (IsEmpty($DoubleBondPos2)) {
      $DoubleBondPos2 = $DoubleBondPos1 + 1;
    }
    if (!exists $DoubleBondsPosInfo{$DoubleBondPos1}) {
      %{$DoubleBondsPosInfo{$DoubleBondPos1}} = ();
    }
    $DoubleBondsPosInfo{$DoubleBondPos1}{$DoubleBondPos2} = "${DoubleBondPos1}(${DoubleBondPos2})";
  }
  $AbbrevTemplateDataMapRef->{DoubleBondsPosInfo} = \%DoubleBondsPosInfo;


  # Setup template structure data lines...
  my(@CmpdLines) = ();
  @CmpdLines = split /\n/, $CmpdString;

  @{$AbbrevTemplateDataMapRef->{CmpdLines}} = ();
  push @{$AbbrevTemplateDataMapRef->{CmpdLines}}, @CmpdLines;

  my($AtomCount, $BondCount) = LMAPSStr::ParseCmpdCountsLine($CmpdLines[3]);
  $AbbrevTemplateDataMapRef->{AtomCount} = $AtomCount;
  $AbbrevTemplateDataMapRef->{BondCount} = $BondCount;

  my(@AtomBlockLines) = ();
  for $Index (4 .. ($AtomCount + 3)) {
    push @AtomBlockLines, $CmpdLines[$Index];
  }
  @{$AbbrevTemplateDataMapRef->{AtomBlockLines}} = ();
  push @{$AbbrevTemplateDataMapRef->{AtomBlockLines}}, @AtomBlockLines;

  my(@BondBlockLines) = ();
  my($FirstAtomNum, $SecondAtomNum, $BondType, $BondStereo, $BondLine, $AbbrevModifier, $CmpdLine);
  for $Index ( ($AtomCount + 4) .. ($AtomCount + $BondCount + 3)) {
    # Parsing the bond block lines and regenerate it to make the bond lines consistent...
    $CmpdLine = $CmpdLines[$Index];
    ($FirstAtomNum, $SecondAtomNum, $BondType, $BondStereo) = LMAPSStr::ParseCmpdBondLine($CmpdLine);

    if (exists $DoubleBondsPosInfo{$FirstAtomNum}) {
      if (exists $DoubleBondsPosInfo{$FirstAtomNum}{$SecondAtomNum}) {
	$BondType = 2;
      }
    }
    elsif ($DoubleBondsPosInfo{$SecondAtomNum}) {
      if ($DoubleBondsPosInfo{$SecondAtomNum}{$FirstAtomNum}) {
	$BondType = 2;
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

  # Setup substituents information...
  my($SubstituentsCount, $SubstituentsBondCount, %SubstituentsPosInfo);
  $SubstituentsCount = 0; $SubstituentsBondCount = 0;
  %SubstituentsPosInfo = ();
  if (!IsEmpty($Substituents)) {
    my($Substituent, $SubstituentPos, $SubstituentAbbrev, $SubstituentStereoChemistry, @SubstituentsList);
    (@SubstituentsList) = split /\//, $Substituents;
    for $Substituent (@SubstituentsList) {
      ($SubstituentPos, $SubstituentStereoChemistry, $SubstituentAbbrev) = ParseSTSubstituentAbbrev($Substituent);
      $SubstituentsCount++;
      $SubstituentsBondCount++;
      if ($SubstituentAbbrev =~ /^Ep$/i) {
	$SubstituentsBondCount++;
      }
      if ($SubstituentStereoChemistry) {
	$SubstituentStereoChemistry = LMAPSStr::StandardizeRingStereochemistrySpec($SubstituentStereoChemistry);
      }
      if (exists $SubstituentsPosInfo{$SubstituentPos}) {
	  # Multiple substituents at one position...
	  $SubstituentsPosInfo{$SubstituentPos}{Count} += 1;
	  push @{$SubstituentsPosInfo{$SubstituentPos}{Spec}}, $Substituent;
	  push @{$SubstituentsPosInfo{$SubstituentPos}{Abbrev}}, $SubstituentAbbrev;
	  push @{$SubstituentsPosInfo{$SubstituentPos}{StereoChemistry}}, $SubstituentStereoChemistry;
      }
      else {
	  %{$SubstituentsPosInfo{$SubstituentPos}} = ();

	  $SubstituentsPosInfo{$SubstituentPos}{Count} = 1;
	  @{$SubstituentsPosInfo{$SubstituentPos}{Spec}} = ();
	  push @{$SubstituentsPosInfo{$SubstituentPos}{Spec}}, $Substituent;

	  @{$SubstituentsPosInfo{$SubstituentPos}{Abbrev}} = ();
	  push @{$SubstituentsPosInfo{$SubstituentPos}{Abbrev}}, $SubstituentAbbrev;

	  @{$SubstituentsPosInfo{$SubstituentPos}{StereoChemistry}} = ();
	  push @{$SubstituentsPosInfo{$SubstituentPos}{StereoChemistry}}, $SubstituentStereoChemistry;
      }
    }
  }

  $AbbrevTemplateDataMapRef->{SubstituentsCount} = $SubstituentsCount;
  $AbbrevTemplateDataMapRef->{SubstituentsBondCount} = $SubstituentsBondCount;
  $AbbrevTemplateDataMapRef->{SubstituentsPosInfo} = \%SubstituentsPosInfo;
}

# Get template ID...
sub _GetTemplateID {
  my($STAbbrev, $STType) = @_;
  my($TemplateID);

  $TemplateID = $STType;

  return $TemplateID;
}

# Generate appropriate compound data string containing structure data and ontology...
sub _GenerateCmdDataString {
  my($Abbrev) = @_;
  my($CmpdDataString, $AtomBlockLinesRef, $BondBlockLinesRef, $CmpdAbbrevTemplateDataMapRef, $OntologyDataLinesRef, @CmpdDataLines);

  $CmpdDataString = '';
  @CmpdDataLines = ();

  # Setup template data for a specific compound abbreviation...
  $CmpdAbbrevTemplateDataMapRef = SetupSTCmpdAbbrevTemplateDataMap($Abbrev);

  # Generate structure data ...
  ($AtomBlockLinesRef, $BondBlockLinesRef) = GenerateSTStrData($CmpdAbbrevTemplateDataMapRef);

  # Generate data block lines including various desriptors...
  ($OntologyDataLinesRef) = GenerateCmpdOntologySDDataLines($CmpdAbbrevTemplateDataMapRef);

  # Write out first four SD file lines: Name, MiscInfo, Comments, Count
  $CmpdDataString .= "$Abbrev\n";
  $CmpdDataString .= LMAPSStr::GenerateCmpdMiscInfoLine(). "\n";
  $CmpdDataString .= LMAPSStr::GenerateCmpdCommentsLine() . "\n";
  $CmpdDataString .= GenerateSTCmpdCountsLine($CmpdAbbrevTemplateDataMapRef) . "\n";

  # Atom lines for template...
  $CmpdDataString .= GenerateSTAtomBlockLines($CmpdAbbrevTemplateDataMapRef, $AtomBlockLinesRef) . "\n";

  # Bond lines for template...
  $CmpdDataString .= GenerateSTBondBlockLines($CmpdAbbrevTemplateDataMapRef, $BondBlockLinesRef) . "\n";

  # Write out any template data block lines: it contains "M  END" as well...
  $CmpdDataString .= JoinWords(\@{$CmpdAbbrevTemplateDataMapRef->{DataBlockLines}}, "\n", 0) . "\n";

  # Write out data block lines including various desriptors...
  $CmpdDataString .= JoinWords($OntologyDataLinesRef, "\n", 0) . "\n";

  $CmpdDataString .= "\$\$\$\$";

  return $CmpdDataString;
}


# Initialize ST data...
sub _InitializeData {
  _InitializeStrTemplates();
  _InitializeSupportedTemplateTypeDataMap();
}

# Initialize structure template data for these supported templates:
#
sub _InitializeStrTemplates {
  %STTemplatesDataMap = ();

  # Cholestane
  my($CholestaneTemplateString)=<<ENDTEMPLATE;
CHOLESTANE structure template
  LipdMAPS02280709152D

 32 35  0  0  0  0  0  0  0  0999 V2000
   -2.8581   -0.4756    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.5725   -0.8882    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.5725   -1.7132    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.8581   -2.1258    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1435   -1.7132    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4291   -2.1258    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7144   -1.7132    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7144   -0.8882    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4291   -0.4756    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1435   -0.8882    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4291    0.3494    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7144    0.7620    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0001    0.3494    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0001   -0.4756    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4291   -0.4756    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4291    0.3494    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7145    0.7620    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0001    0.9889    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1435   -0.3105    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7145    1.4220    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0426    1.6020    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2861    1.7520    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.8577    1.4220    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.4293    1.7520    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0009    1.4220    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.5725    1.7520    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0009    0.8648    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7144   -0.1868    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0001   -1.1158    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.1772    1.0291    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.7145    2.1258    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4291   -1.1022    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
  3  4  1  0  0  0  0
  4  5  1  0  0  0  0
  5  6  1  0  0  0  0
  6  7  1  0  0  0  0
  7  8  1  0  0  0  0
  8  9  1  0  0  0  0
  9 11  1  0  0  0  0
 11 12  1  0  0  0  0
 12 13  1  0  0  0  0
 13 14  1  0  0  0  0
 14 15  1  0  0  0  0
 15 16  1  0  0  0  0
 16 17  1  0  0  0  0
 17 20  1  0  0  0  0
 20 22  1  0  0  0  0
 22 23  1  0  0  0  0
 23 24  1  0  0  0  0
 24 25  1  0  0  0  0
 25 26  1  0  0  0  0
 25 27  1  0  0  0  0
  1 10  1  0  0  0  0
  5 10  1  0  0  0  0
  9 10  1  0  0  0  0
  8 14  1  0  0  0  0
 13 17  1  0  0  0  0
 13 18  1  1  0  0  0
 10 19  1  1  0  0  0
 20 21  1  6  0  0  0
  8 28  1  1  0  0  0
 14 29  1  6  0  0  0
 17 30  1  6  0  0  0
 20 31  1  1  0  0  0
  9 32  1  6  0  0  0
M  END

ENDTEMPLATE

  # Ergostane
  my($ErgostaneTemplateString)=<<ENDTEMPLATE;
ERGOSTANE structure template
  LipdMAPS02280709152D

 33 36  0  0  0  0  0  0  0  0999 V2000
   -2.8581   -0.6187    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.5729   -1.0313    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.5729   -1.8563    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.8581   -2.2689    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1435   -1.8563    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4291   -2.2689    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7144   -1.8563    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7144   -1.0313    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4291   -0.6187    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1435   -1.0313    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4291    0.2063    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7144    0.6189    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0001    0.2063    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0001   -0.6187    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4291   -0.6187    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4291    0.2063    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7145    0.6189    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0001    0.8459    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1435   -0.4536    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7145    1.2790    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0426    1.4590    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2861    1.6090    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.8577    1.2790    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.4293    1.6090    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0009    1.2790    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.5729    1.6090    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0009    0.7217    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7144   -0.3299    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0001   -1.2589    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.1772    0.8861    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.7145    1.9828    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4291   -1.2453    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.4293    2.2689    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  1  0      
  5  6  1  0      
  6  7  1  0      
  7  8  1  0      
  8  9  1  0      
  9 11  1  0      
 11 12  1  0      
 12 13  1  0      
 13 14  1  0      
 14 15  1  0      
 15 16  1  0      
 16 17  1  0      
 17 20  1  0      
 20 22  1  0      
 22 23  1  0      
 23 24  1  0      
 24 25  1  0      
 25 26  1  0      
 25 27  1  0      
  1 10  1  0      
  5 10  1  0      
  9 10  1  0      
  8 14  1  0      
 13 17  1  0      
 13 18  1  1      
 10 19  1  1      
 20 21  1  6      
  8 28  1  1      
 14 29  1  6      
 17 30  1  6      
 20 31  1  1      
  9 32  1  6      
 24 33  1  6      
M  END

ENDTEMPLATE

  # Campestane
  my($CampestaneTemplateString)=<<ENDTEMPLATE;
CAMPESTANE structure template
  LipdMAPS02280709152D

 33 36  0  0  0  0  0  0  0  0999 V2000
   -2.8581   -0.6187    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.5728   -1.0313    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.5728   -1.8563    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.8581   -2.2689    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1435   -1.8563    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4291   -2.2689    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7144   -1.8563    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7144   -1.0313    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4291   -0.6187    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1435   -1.0313    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4291    0.2063    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7144    0.6189    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0001    0.2063    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0001   -0.6187    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4291   -0.6187    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4291    0.2063    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7145    0.6189    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0001    0.8459    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1435   -0.4536    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7145    1.2790    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0426    1.4590    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2861    1.6090    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.8577    1.2790    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.4293    1.6090    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0009    1.2790    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.5728    1.6090    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0009    0.7217    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7144   -0.3299    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0001   -1.2589    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.1772    0.8861    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.7145    1.9828    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4291   -1.2453    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.4293    2.2689    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  1  0      
  5  6  1  0      
  6  7  1  0      
  7  8  1  0      
  8  9  1  0      
  9 11  1  0      
 11 12  1  0      
 12 13  1  0      
 13 14  1  0      
 14 15  1  0      
 15 16  1  0      
 16 17  1  0      
 17 20  1  0      
 20 22  1  0      
 22 23  1  0      
 23 24  1  0      
 24 25  1  0      
 25 26  1  0      
 25 27  1  0      
  1 10  1  0      
  5 10  1  0      
  9 10  1  0      
  8 14  1  0      
 13 17  1  0      
 13 18  1  1      
 10 19  1  1      
 20 21  1  6      
  8 28  1  1      
 14 29  1  6      
 17 30  1  6      
 20 31  1  1      
  9 32  1  6      
 24 33  1  1      
M  END

ENDTEMPLATE

  # Stigmastane
  my($StigmastaneTemplateString)=<<ENDTEMPLATE;
STIGMASTANE structure template
  LipdMAPS02280709152D

 34 37  0  0  0  0  0  0  0  0999 V2000
   -2.8581   -0.7940    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.5729   -1.2066    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.5729   -2.0316    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.8581   -2.4442    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1435   -2.0316    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4291   -2.4442    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7144   -2.0316    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7144   -1.2066    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4291   -0.7940    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1435   -1.2066    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4291    0.0310    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7144    0.4436    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0001    0.0310    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0001   -0.7940    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4291   -0.7940    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4291    0.0310    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7145    0.4436    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0001    0.6706    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1435   -0.6289    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7145    1.1037    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0426    1.2837    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2861    1.4337    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.8577    1.1037    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.4293    1.4337    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0009    1.1037    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.5729    1.4337    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0009    0.5464    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7144   -0.5052    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0001   -1.4342    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.1772    0.7108    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.7145    1.8075    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4291   -1.4206    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.4293    2.0936    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0200    2.4442    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  1  0      
  5  6  1  0      
  6  7  1  0      
  7  8  1  0      
  8  9  1  0      
  9 11  1  0      
 11 12  1  0      
 12 13  1  0      
 13 14  1  0      
 14 15  1  0      
 15 16  1  0      
 16 17  1  0      
 17 20  1  0      
 20 22  1  0      
 22 23  1  0      
 23 24  1  0      
 24 25  1  0      
 25 26  1  0      
 25 27  1  0      
  1 10  1  0      
  5 10  1  0      
  9 10  1  0      
  8 14  1  0      
 13 17  1  0      
 13 18  1  1      
 10 19  1  1      
 20 21  1  6      
  8 28  1  1      
 14 29  1  6      
 17 30  1  6      
 20 31  1  1      
  9 32  1  6      
 24 33  1  1      
 33 34  1  0      
M  END

ENDTEMPLATE

  # Format: ID => AbbrevID|NumOfCarbonAtoms|RingJunction1AtomNum1|RingJunction1AtomNum2|RingJunction2AtomNum1|RingJunction2AtomNum2|RingJunction3AtomNum1|RingJunction3AtomNum2|LMCategory|LMMainClass|LMSubClass|TemplateString
  % STTemplatesDataMap = (
		       "CHOLESTANE" => "CHOLESTANE|27|5|10|8|9|13|14|ST|01|01|$CholestaneTemplateString",

		       "ERGOSTANE" => "ERGOSTANE|28|5|10|8|9|13|14|ST|01|03|$ErgostaneTemplateString",

		       "CAMPESTANE" => "CAMPESTANE|28|5|10|8|9|13|14|ST|01|03|$CampestaneTemplateString",

		       "STIGMASTANE" => "STIGMASTANE|29|5|10|8|9|13|14|ST|01|04|$StigmastaneTemplateString",
		      );

}

# Initialize supported template groups...
sub _InitializeSupportedTemplateTypeDataMap {
  my($STTypeID, $STType);
  %STSupportedTemplateTypeMap = ();

  for $STTypeID (keys %STTemplatesDataMap) {
    ($STType) = split /\|/, $STTemplatesDataMap{$STTypeID};
    if (!(exists $STSupportedTemplateTypeMap{$STType})) {
      $STSupportedTemplateTypeMap{$STType} = $STType;
    }
  }
}

1;

__END__

=head1 NAME

STStr - Sterol (ST) structure generation methods

=head1 SYNOPSIS

use STStr;

use STStr qw(:all);

=head1 DESCRIPTION

STStr module provides these methods:

    ExpandSTCmpdAbbrevs - Expand ST abbreviation
    GenerateCmpdOntologyData - Generate ontology data
    GenerateCmpdOntologySDDataLines - Generate ontology data lines for
                                      SD file
    GenerateSTStrData - Generate structure data
    GenerateSDFile - Generate SD file
    IsSTAbbrevSupported - Is it a supported ST abbreviation
    IsSTSubstituentsNameSupported - Is it a supported ST substituent name
    IsSTDoubleBondsAbbrevOkay - Is it a valid ST double bond abbreviation
    IsSTSubstituentsAbbrevOkay - Is it a valid ST substituent abbreviation
    IsWildCardInSTAbbrev - Does ST abbreviatio contains a wild card
    ParseSTAbrev - Parse ST abbreviation
    ParseSTDoubleBondAbbrev - Parse ST double bond abbreviation
    ParseSTSubstituentAbbrev - Parse ST substituent abbreviation
    SetupSTCmpdAbbrevTemplateDataMap - Setup template structure data map
    ValidateSTAbbrev - Validate ST abbreviation

=head1 METHODS

=over 4

=item B<ExpandSTCmpdAbbrevs>

    $ExpandedAbbrevArrayRef = ExpandSTCmpdAbbrevs($CmpdAbbrev);

Return a reference to an array containing complete ST abbreviations. Wild card
characters in ST abbreviation name are expanded to generate fully qualified
ST abbreviations.

=item B<GenerateCmpdOntologyData>

    $DataHashRef = GenerateCmpdOntologyData($CmpdDataRef);

Return a reference to a hash containing ontology data with hash keys and values
corresponding to property names and values.

=item B<GenerateCmpdOntologySDDataLines>

    $DataLinesArrayRef =
        GenerateCmpdOntologySDDataLines($CmpDataRef);

Return a reference to an array containing ontology data lines suitable for
generate SD file data block.

=item B<GenerateSTStrData>

    ($AtomLinesArrayRef, $BondLinesArrayRef) =
       GenerateSTStrData($CmpdDataRef);

Return array references containing atom and bond data lines for SD file. Appropriate atom
and bond data lines are generated using abbreviation template data.

=item B<GenerateSDFile>

    GenerateSDFile($SDFileName, $CmdAbbrevsRef);

Generate a SD file for compound abbreviations. Structure data for specified abbreviation
is generated sequentially and written to SD file.

=item B<IsSTAbbrevSupported>

    $Status = IsSTAbbrevSupported($Abbrev);

Return 1 or 0 based on whether ST abbreviation is supported.

=item B<IsSTSubstituentsNameSupported>

    $Status = IsSTSubstituentsNameSupported($SubstituentAbbrev);

Return 1 or 0 based on whether ST substituent abbreviation is supported.

=item B<IsSTDoubleBondsAbbrevOkay>

    $Status = IsSTDoubleBondsAbbrevOkay($STAbbrev, $STType,
       $SubstituentsAbbrev, $DoubleBondsAbbrev);

Return 1 or 0 based on whether ST double bond abbreviation is valid.

=item B<IsSTSubstituentsAbbrevOkay>

    $Status = IsSTSubstituentsAbbrevOkay($STAbbrev, $STType,
       $SubstituentsAbbrev, $DoubleBondsAbbrev);

Return 1 or 0 based on whether ST substituent abbreviation is valid.

=item B<IsWildCardInSTAbbrev>

    $Status = IsSTAbbrevSupported($Abbrev);

Return 1 or 0 based on whether ST abbreviation contains wild card.

=item B<ParseSTAbbrev>

    ($STType,  $SubstituentsAbbrev, $DoubleBondsAbbrev) =
        ParseSTAbrev($Abbrev);

Parse ST abbreviation and return these values: STType,  SubstituentsAbbrev, and
DoubleBondsAbbrev.

=item B<ParseSTDoubleBondAbbrev>

    ($BondPos1, $BondPos1) = ParseSTDoubleBondAbbrev($Abbrev);

Parse ST double bond abbreviation and return these values: BondPos1 and BondPos2.

=item B<ParseSTSubstituentAbbrev>

    ($SubstituentPos, $SubstituentAbbrev, $StereoChemistry) =
       ParseSTSubstituentAbbrev($Abbrev);

Parse ST substituents abbreviation and return these values: SubstituentPos, SubstituentAbbrev,
and SubstituentStereoChemistry.

=item B<SetupSTCmpdAbbrevTemplateDataMap>

    $AbbrevTemplateDataMapRef =
       SetupSTCmpdAbbrevTemplateDataMap($Abbrev);

Return a reference to a hash containing template data for compound abbreviation. The
template data is used to generate SD file for compound abbreviation.

=item B<ValidateSTAbbrev>

    $Status = ValidateSTAbbrev($Abbrev);

Return 1 or 0 based on whether a ST abbreviation is valid.

=back

=head1 AUTHOR

Manish Sud

=head1 CONTRIBUTOR

Eoin Fahy

=head1 SEE ALSO

ChainStr.pm, LMAPSStr.pm

=head1 COPYRIGHT

Copyright (C) 2006-2012. The Regents of the University of California. All Rights Reserved.

=head1 LICENSE

Modified BSD License

=cut
