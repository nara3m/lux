package LMStr;
#
# File: LMStr.pm
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
use TextUtil;
use ChainAbbrev;
use ChainStr;

our ($VERSION, @ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);

$VERSION = '2.00';
@ISA = qw(Exporter);
@EXPORT = qw();
@EXPORT_OK = qw(ExpandLMCmpdAbbrevs GenerateCmpdOntologyData GenerateCmpdOntologySDDataLines GenerateLMChainStrData GenerateSDFile GetLMTemplatesData GetLMSupportedHeadGroupMap GetLMTemplateID IsLMChainsAbbrevSupported ParseLMAbbrev SetupLMCmpdAbbrevTemplateDataMap ValidateLMAbbrev);
%EXPORT_TAGS = (
		all  => [@EXPORT, @EXPORT_OK],
		);

# Initialize LM data...
my(%LMTemplatesDataMap, %LMSupportedHeadGroupMap);
_InitializeData();

#
# Process compound abbreviation containing upto four acyl chains with
# no support for wild card and return a valid list of abbreviations.
#
sub ExpandLMCmpdAbbrevs {
  my($CmpdAbbrevsRef, $Abbrev, $HeadGroupAbbrev, $ChainsAbbrev, $AbbrevModifier, $ChainAbbrev, $AllowSubstituents, $AllowRings, $AllowArbitraryChainLenSpec, $NewAbbrev, @ChainsAbbrevList, @ExpandedCmpdAbbrevs);

  ($CmpdAbbrevsRef) = @_;

  @ExpandedCmpdAbbrevs = ();

 ABBREV: for $Abbrev (@{$CmpdAbbrevsRef}) {

    if (!($Abbrev =~ /\(/ && $Abbrev =~ /\)/)) {
      warn "Warning: Ignored LM compound abbreviation $Abbrev due to incorrect format. Valid format: HeadGroupAbbrev(Chain1Abbrev), HeadGroupAbbrev(Chain1Abbrev/Chain2Abbrev), HeadGroupAbbrev(Chain1Abbrev/Chain2Abbrev/Chain3Abbrev) or HeadGroupAbbrev(Chain1Abbrev/Chain2Abbrev/Chain3Abbrev/Chain4Abbrev)\n";
      next ABBREV;
    }
    if (!ValidateLMAbbrev($Abbrev)) {
      next ABBREV;
    }

    ($HeadGroupAbbrev, $ChainsAbbrev, $AbbrevModifier) = ParseLMAbbrev($Abbrev);

    @ChainsAbbrevList = ();
    @ChainsAbbrevList = split /\//, $ChainsAbbrev;
    if (@ChainsAbbrevList < 1) {
      warn "Warning: Ignored LM compound Abbreviation $Abbrev due to incorrect format: Must contain chain abbreviation for at least one chain\n";
      next ABBREV;
    }

    if (@ChainsAbbrevList > 4)  {
      warn "Warning: Ignored LM compound Abbreviation $Abbrev due to incorrect format: Number of chain abbreviations must be <= 4\n";
      next ABBREV;
    }

    ($AllowSubstituents, $AllowRings, $AllowArbitraryChainLenSpec) = (1, 0, 1);
    for $ChainAbbrev (@ChainsAbbrevList) {
      if (!(ChainAbbrev::IsChainAbbrevOkay($ChainAbbrev, $AllowSubstituents, $AllowRings, $AllowArbitraryChainLenSpec))) {
	warn "Warning: Ignoring LM compound abbreviation $HeadGroupAbbrev($ChainsAbbrev)\n";
	next ABBREV;
      }
    }

    if (!IsLMChainsAbbrevSupported($ChainsAbbrev, 1)) {
      warn "Warning: Ignoring LM compound abbreviation $HeadGroupAbbrev($ChainsAbbrev)\n";
      next ABBREV;
    }

    my($WildCardInHeadGroup) = ($HeadGroupAbbrev =~ /\*/ ) ? 1 : 0;
    my($WildCardInChains) = 0;

    for $ChainAbbrev (@ChainsAbbrevList) {
      if (ChainAbbrev::IsWildCardInChainAbbrev($ChainAbbrev)) {
	$WildCardInChains = 1;
      }
    }
    if (!($WildCardInHeadGroup || $WildCardInChains)) {
      my($TemplateHeadGroup) = GetLMTemplateID($HeadGroupAbbrev, $ChainsAbbrev);
      if (exists $LMTemplatesDataMap{$TemplateHeadGroup}) {
	$NewAbbrev = "$HeadGroupAbbrev($ChainsAbbrev)$AbbrevModifier";
	push @ExpandedCmpdAbbrevs, $NewAbbrev;
      }
      else {
	warn "Warning: Ignored LM compound abbreviation $Abbrev : Abbreviation doesn't match any template\n";
      }
      next ABBREV;
    }

    #
    # Wild cards in abbreaviations are currently not supported.
    warn "Warning: Ignored LM compound abbreviation $Abbrev : Abbreviation contains wild cards\n";
    next ABBREV;

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

  return \%OntologyDataMap;
}

# Abbreviation and systematic name...
sub SetupAbbrevAndSystematicName {
  my($CmpdAbbrevTemplateDataMapRef, $OntologyDataMapRef) = @_;
  my($SnChainIndex, $SnAbbrev, $HeadGroupAbbrev, $AbbrevModifier, $Abbreviation, @NonEmptyAbbrevs, @SnChainPositions, @SnAbbrevs);

  # Skip setting of systematic names for LM generic structure module...

  @SnChainPositions = qw(1 2 3 4);
  @SnAbbrevs = ();

  for $SnChainIndex (0 .. $#SnChainPositions) {
    push @SnAbbrevs, $CmpdAbbrevTemplateDataMapRef->{SnChainAdd}[$SnChainIndex] ? $CmpdAbbrevTemplateDataMapRef->{SnAbbrev}[$SnChainIndex] : '0:0';
  }

  $HeadGroupAbbrev = $CmpdAbbrevTemplateDataMapRef->{AbbrevID};
  $AbbrevModifier = $CmpdAbbrevTemplateDataMapRef->{AbbrevModifier};

  @NonEmptyAbbrevs = ();
  for $SnChainIndex (0 .. $#SnChainPositions) {
    $SnAbbrev = $SnAbbrevs[$SnChainIndex];
    if ($CmpdAbbrevTemplateDataMapRef->{SnChainAdd}[$SnChainIndex] && $SnAbbrev !~ /^0:0$/i) {
      push @NonEmptyAbbrevs, $SnAbbrev;
    }
  }
  $Abbreviation = "$HeadGroupAbbrev(" . join("/", @NonEmptyAbbrevs) . ")$AbbrevModifier";

  $OntologyDataMapRef->{Name} = $Abbreviation;
  $OntologyDataMapRef->{Abbrev} = $Abbreviation;
}

# LM classification info...
sub SetupLMClassificationInfo {
  my($CmpdAbbrevTemplateDataMapRef, $OntologyDataMapRef) = @_;

  $OntologyDataMapRef->{'LM Category'} = $CmpdAbbrevTemplateDataMapRef->{LMCategory};
  if ($CmpdAbbrevTemplateDataMapRef->{LMMainClass}) {
    $OntologyDataMapRef->{'LM Main Class'} = $CmpdAbbrevTemplateDataMapRef->{LMCategory} . $CmpdAbbrevTemplateDataMapRef->{LMMainClass};
  }
  if ($CmpdAbbrevTemplateDataMapRef->{LMMainClass} && $CmpdAbbrevTemplateDataMapRef->{LMSubClass} ) {
    $OntologyDataMapRef->{'LM Sub Class'} = $CmpdAbbrevTemplateDataMapRef->{LMCategory} . $CmpdAbbrevTemplateDataMapRef->{LMMainClass} . $CmpdAbbrevTemplateDataMapRef->{LMSubClass};
  }
}

# Generate a SD file...
sub GenerateSDFile {
  my($SDFileName, $CmdAbbrevsRef) = @_;
  my($Abbrev, $CmpdDataString);

  open SDFILE, ">$SDFileName" or die "Error: Couldn't open $SDFileName: $! \n";

  for $Abbrev (@$CmdAbbrevsRef) {

    # Set up a compound data string...
    $CmpdDataString = '';
    $CmpdDataString = _GenerateCmdDataString($Abbrev);

    # Write it out to the file...
    print SDFILE "$CmpdDataString\n";

  }
  close SDFILE;
}

# Generate atom and bond block lines for a chain
sub GenerateLMChainStrData {
  my($ChainType, $CmpdAbbrevTemplateDataMapRef) = @_;

  return ChainStr::GenerateChainStrData($ChainType, $CmpdAbbrevTemplateDataMapRef);
}

# Return reference to LMTemplatesDataMap...
sub GetLMTemplatesData {
  return \%LMTemplatesDataMap;
}

# Return reference to LMSupportedHeadGroupMap...
sub GetLMSupportedHeadGroupMap {
  return \%LMSupportedHeadGroupMap;
}

# Get GP template ID...
sub GetLMTemplateID {
  my($HeadGroupAbbrev, $ChainsAbbrev) = @_;
  my($HeadGroupID, $ChainAbbrev, $ChainNum);

  $HeadGroupID = "";
  if (!$HeadGroupAbbrev) {
    return $HeadGroupID;
  }

  $ChainNum = 0;
  ABBREV: for $ChainAbbrev (split /\//, $ChainsAbbrev) {
    $ChainNum++;
    if ($ChainAbbrev eq "0:0") {
      next ABBREV;
    }
    $HeadGroupID .= "Sn${ChainNum}";
    if (ChainAbbrev::IsAlkylChainAbbrev($ChainAbbrev)) {
      $HeadGroupID .= "Alkyl";
    }
    elsif (ChainAbbrev::IsAlkenylChainAbbrev($ChainAbbrev)) {
      $HeadGroupID .= "Alkenyl";
    }
  }
  $HeadGroupID = "${HeadGroupAbbrev}${HeadGroupID}";

  return $HeadGroupID;
}

# Does template exist to handle this abbreviation?
#
sub IsLMChainsAbbrevSupported {
  my($Abbrev, $PrintWarning) = @_;
  my(@AbbrevList);

  @AbbrevList = split /\//, $Abbrev;
  if (@AbbrevList < 1) {
    if ($PrintWarning) {
      warn "Warning: Ignoring LM compound abbreviation $Abbrev : Must contain chain abbreviation for at least one chain\n";
    }
    return 0;
  }

  if (@AbbrevList > 4) {
    if ($PrintWarning) {
      warn "Warning: Ignoring LM compound abbreviation $Abbrev : Number of chain abbreviations must be <= 4\n";
    }
    return 0;
  }

  return 1;
}

# Parse abbreviation...
sub ParseLMAbbrev {
  my($Abbrev) = @_;
  my($HeadGroup, $ChainsAbbrev, $AbbrevModifier);

  if ($Abbrev =~ /\]$/) {
    ($HeadGroup, $ChainsAbbrev, $AbbrevModifier) = $Abbrev =~ /^(.+?)\((.+?)\)\[(.+?)\]$/;
    $AbbrevModifier = '[' . $AbbrevModifier . ']';
  }
  else {
    ($HeadGroup, $ChainsAbbrev) = $Abbrev =~ /^(.+?)\((.+?)\)$/;
    $AbbrevModifier = '';
  }

  return ($HeadGroup, $ChainsAbbrev, $AbbrevModifier);
}

#
# Check out the validity of LM abbreviation...
sub ValidateLMAbbrev {
  my($Abbrev) = @_;

  if (!($Abbrev =~ /\(/ && $Abbrev =~ /\)/)) {
    warn "Warning: Ignored LM compound abbreviation $Abbrev due to incorrect format. Valid format: HeadGroupAbbrev(Chain1Abbrev), HeadGroupAbbrev(Chain1Abbrev/Chain2Abbrev), HeadGroupAbbrev(Chain1Abbrev/Chain2Abbrev/Chain3Abbrev) or HeadGroupAbbrev(Chain1Abbrev/Chain2Abbrev/Chain3Abbrev/Chain4Abbrev)\n";
    return 0;
  }

  # Make sure head group is these...
  if ($Abbrev =~ /^\(/) {
    warn "Warning: Ignored LM compound abbreviation $Abbrev due to incorrect format. Valid format: HeadGroupAbbrev(Chain1Abbrev), HeadGroupAbbrev(Chain1Abbrev/Chain2Abbrev), HeadGroupAbbrev(Chain1Abbrev/Chain2Abbrev/Chain3Abbrev) or HeadGroupAbbrev(Chain1Abbrev/Chain2Abbrev/Chain3Abbrev/Chain4Abbrev)\n";
    return 0;
  }

  # Number of '(' and ')' parenthesis must match...
  my(@LeftParenthesis) = $Abbrev =~ /\(/g;
  my(@RightParenthesis) = $Abbrev =~ /\)/g;
  if ($#LeftParenthesis != $#RightParenthesis) {
    warn "Warning: Ignored LM compound abbreviation $Abbrev due to incorrect format: Number of ( and ) doesn't match. Valid format: HeadGroupAbbrev(Chain1Abbrev), HeadGroupAbbrev(Chain1Abbrev/Chain2Abbrev), HeadGroupAbbrev(Chain1Abbrev/Chain2Abbrev/Chain3Abbrev) or HeadGroupAbbrev(Chain1Abbrev/Chain2Abbrev/Chain3Abbrev/Chain4Abbrev)\n";
    return 0;
  }

  @LeftParenthesis = $Abbrev =~ /\[/g;
  @RightParenthesis = $Abbrev =~ /\]/g;
  if ($#LeftParenthesis != $#RightParenthesis) {
    warn "Warning: Ignored LM compound abbreviation $Abbrev due to incorrect format: Number of [ and ] doesn't match. Valid format: HeadGroupAbbrev(Chain1Abbrev), HeadGroupAbbrev(Chain1Abbrev/Chain2Abbrev), HeadGroupAbbrev(Chain1Abbrev/Chain2Abbrev/Chain3Abbrev) or HeadGroupAbbrev(Chain1Abbrev/Chain2Abbrev/Chain3Abbrev/Chain4Abbrev)\n";
    return 0;
  }

  my($HeadGroup, $ChainsAbbrev, $AbbrevModifier) = ParseLMAbbrev($Abbrev);

  if ($HeadGroup !~ /\*/ ) {
    if (!(exists $LMSupportedHeadGroupMap{$HeadGroup})) {
      warn "Warning: Ignored LM compound abbreviation $Abbrev) : Unknown head group $HeadGroup.\n";
      return 0;
    }
  }

  my(@ChainsAbbrevList) = ();
  @ChainsAbbrevList = split /\//, $ChainsAbbrev;

  if (@ChainsAbbrevList < 1) {
    warn "Warning: Ignoring LM compound abbreviation $Abbrev : Must contain chain abbreviation for at least one chain\n";
    return 0;
  }

  if (@ChainsAbbrevList > 4) {
    warn "Warning: Ignoring LM compound abbreviation $Abbrev : Number of chain abbreviations must be <= 4\n";
    return 0;
  }

  if ($AbbrevModifier) {
    if (!($AbbrevModifier =~ /\[/ && $AbbrevModifier =~ /\]/)) {
      warn "Warning: Ignored LM compound abbreviation $Abbrev due to incorrect stereochemistry abbreviation: It must be enclosed by []\n";
      return 0;
    }
    $AbbrevModifier =~ s/(\[|\])//g;
    if ($AbbrevModifier !~ /^(rac|R|S|U)$/) {
      warn "Warning: Ignored LM compound abbreviation $Abbrev due to incorrect stereochemistry abbreviation: Valid values: rac, R, S, U\n";
      return 0;
    }
  }
  return 1;
}


# Generate appropriate compound data string containing structure data...
sub _GenerateCmdDataString {
  my($Abbrev) = @_;
  my($CmpdDataString, $CmpdAbbrevTemplateDataMapRef);

  $CmpdDataString = '';

  # Setup template data for a specific compound abbreviation...
  $CmpdAbbrevTemplateDataMapRef = SetupLMCmpdAbbrevTemplateDataMap($Abbrev);

  # Generate structure data for chains...
  my($Sn1AtomBlockLinesRef, $Sn1BondBlockLinesRef) = LMStr::GenerateLMChainStrData('Sn1', $CmpdAbbrevTemplateDataMapRef);
  my($Sn2AtomBlockLinesRef, $Sn2BondBlockLinesRef) = LMStr::GenerateLMChainStrData('Sn2', $CmpdAbbrevTemplateDataMapRef);
  my($Sn3AtomBlockLinesRef, $Sn3BondBlockLinesRef) = LMStr::GenerateLMChainStrData('Sn3', $CmpdAbbrevTemplateDataMapRef);
  my($Sn4AtomBlockLinesRef, $Sn4BondBlockLinesRef) = LMStr::GenerateLMChainStrData('Sn4', $CmpdAbbrevTemplateDataMapRef);

  # Generate data block lines including various desriptors...
  my($OntologyDataLinesRef) = GenerateCmpdOntologySDDataLines($CmpdAbbrevTemplateDataMapRef);

  # Write out first four SD file lines: Name, MiscInfo, Comments, Count
  $CmpdDataString .= "$Abbrev\n";
  $CmpdDataString .= LMAPSStr::GenerateCmpdMiscInfoLine(). "\n";
  $CmpdDataString .= LMAPSStr::GenerateCmpdCommentsLine() . "\n";
  $CmpdDataString .= ChainStr::GenerateCmpdCountsLine($CmpdAbbrevTemplateDataMapRef) . "\n";

  # Atom lines for template and chains...
  $CmpdDataString .= ChainStr::GenerateAtomBlockLines($CmpdAbbrevTemplateDataMapRef, $Sn1AtomBlockLinesRef, $Sn2AtomBlockLinesRef, $Sn3AtomBlockLinesRef, $Sn4AtomBlockLinesRef) . "\n";

  # Bond lines for template and chains...
  $CmpdDataString .= ChainStr::GenerateBondBlockLines($CmpdAbbrevTemplateDataMapRef, $Sn1BondBlockLinesRef, $Sn2BondBlockLinesRef, $Sn3BondBlockLinesRef, $Sn4BondBlockLinesRef) . "\n";

  # Write out any template data block lines: it contains "M  END" as well...
  $CmpdDataString .= JoinWords(\@{$CmpdAbbrevTemplateDataMapRef->{DataBlockLines}}, "\n", 0) . "\n";

  # Write out data block lines including various desriptors...
  $CmpdDataString .= JoinWords($OntologyDataLinesRef, "\n", 0) . "\n";

  $CmpdDataString .= "\$\$\$\$";

  # Assign any specified stereochemistry to substituents...
  if (ChainStr::IsAnySubstituentSpecifiedWithStereoChemistry($CmpdAbbrevTemplateDataMapRef)) {
    my(@CmpdDataLines) = split /\n/, $CmpdDataString;

    ChainStr::AssignSubstituentStereoChemistry($CmpdAbbrevTemplateDataMapRef, \@CmpdDataLines);
    $CmpdDataString = '';
    $CmpdDataString = JoinWords(\@CmpdDataLines, "\n", 0);
  }

  return $CmpdDataString;
}

# Set up template data for a specific cmpd abbreviation and return a
# reference to the hash...
sub SetupLMCmpdAbbrevTemplateDataMap {
  my($AbbrevHeadGroup, $Abbrev, $ChainsAbbrev, $AbbrevModifier, $TemplateID, @SnAbbrev, @SnChainAdd, %AbbrevTemplateDataMap);

  ($Abbrev) = @_;
  %AbbrevTemplateDataMap = ();

  ($AbbrevHeadGroup, $ChainsAbbrev, $AbbrevModifier) = ParseLMAbbrev($Abbrev);

  $TemplateID = GetLMTemplateID($AbbrevHeadGroup, $ChainsAbbrev);

  # Set up up to four sn abbreviations...
  @SnAbbrev = split /\//, $ChainsAbbrev;
  if (@SnAbbrev < 4) {
    my($SnAbbrevsToAdd);
    $SnAbbrevsToAdd = 4 - @SnAbbrev;
    push @SnAbbrev, ("") x $SnAbbrevsToAdd;
  }

  @SnChainAdd = map { ($_ && $_ ne "0:0") ? 1 : 0 } @SnAbbrev;

  @{$AbbrevTemplateDataMap{SnChainAdd}} = ();
  push @{$AbbrevTemplateDataMap{SnChainAdd}}, @SnChainAdd;

  @{$AbbrevTemplateDataMap{SnAbbrev}} = ();
  push @{$AbbrevTemplateDataMap{SnAbbrev}}, @SnAbbrev;

  $AbbrevModifier =~ s/(\[|\])//g;
  $AbbrevTemplateDataMap{AbbrevModifier} = $AbbrevModifier;

  ChainStr::SetupTemplateDataMap('LM', \%AbbrevTemplateDataMap, $LMTemplatesDataMap{$TemplateID});

  return \%AbbrevTemplateDataMap;
}

# Initialize GP data...
sub _InitializeData {
  _InitializeStrTemplates();
  _InitializeSupportedHeadGroupsData();
}

# Initialize structure template data for these supported templates:
#
sub _InitializeStrTemplates {
  %LMTemplatesDataMap = ();

  my($PIM3Sn1Sn2TemplateString)=<<ENDPIM3SN1SN2TEMPLATE;
PIM3 sn1 acyl and sn2 acyl template structure
  LipdMAPS01111109572D

 61 64  0  0  0  0  0  0  0  0999 V2000
   -4.2672   -0.2270    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.9795    0.1831    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.6921   -0.2270    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -6.4042    0.1831    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.4042    1.0066    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.8553   -0.9393    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -4.6789   -0.9393    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -7.1166   -0.2270    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.5545    0.1844    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.8421   -0.2270    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8283    0.0613    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1912   -0.5679    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8283    0.8116    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.4251   -1.3604    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.4251   -2.1841    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -6.1376   -0.9488    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.2307   -0.5051    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.8263   -1.5409    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0721   -0.5051    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0771   -0.8363    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.5678   -1.2525    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.6786   -1.2126    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5207   -1.5409    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.1105   -0.7408    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.1163   -4.2924    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.5207   -3.2565    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.4191   -4.2924    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.2700   -3.9612    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.1267   -3.8682    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.6284   -3.4381    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.6684   -4.4086    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.5207   -2.4327    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.9541   -4.4217    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.6684   -3.5848    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.8263   -3.2565    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.6319   -5.0887    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.1085   -1.9570    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.6727    1.7559    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.0771    0.7202    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.9755    1.7559    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.8263    1.4248    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.6831    1.3319    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.1848    0.9019    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.2248    1.8721    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.5177   -0.2626    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.5106    1.8853    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.2248    1.0485    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.3827    0.7202    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.1063    4.5687    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.5106    3.5329    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.4090    4.5687    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.2599    4.2376    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.1166    4.1446    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    6.6183    3.7145    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.6583    4.6849    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.5106    2.7091    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.9441    4.6981    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.6583    3.8612    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.8163    3.5329    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.6759    5.0887    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1139   -0.3512    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  2  0      
  4  8  1  0      
  1  7  1  6      
  1  6  1  1      
  9  1  1  0      
 10  9  1  0      
 11 12  1  0      
 11 13  2  0      
 14 15  2  0      
 14 16  1  0      
 14  7  1  0      
 11 10  1  0      
 17 18  1  0      
 18 22  1  1      
 19 23  1  1      
 19 20  1  0      
 17 20  1  0      
 18 21  1  0      
 23 22  1  1      
 20 45  1  0      
 17 24  1  0      
 22 37  1  0      
 23 32  1  0      
 25 26  1  0      
 26 34  1  1      
 27 35  1  1      
 27 28  1  0      
 25 28  1  0      
 27 29  1  0      
 30 35  1  0      
 31 34  1  0      
 26 32  1  0      
 28 33  1  0      
 34 35  1  1      
 33 36  1  0      
 38 39  1  0      
 39 47  1  1      
 40 48  1  1      
 40 41  1  0      
 38 41  1  0      
 40 42  1  0      
 43 48  1  0      
 44 47  1  0      
 39 45  1  0      
 41 46  1  0      
 47 48  1  1      
 46 56  1  0      
 49 50  1  0      
 50 58  1  1      
 51 59  1  1      
 51 52  1  0      
 49 52  1  0      
 51 53  1  0      
 54 59  1  0      
 55 58  1  0      
 50 56  1  0      
 52 57  1  0      
 58 59  1  1      
 57 60  1  0      
 11 61  1  0      
 61 19  1  0      
M  END

ENDPIM3SN1SN2TEMPLATE

  my($DATSn1Sn2TemplateString)=<<ENDDATSN1SN2TEMPLATE;
DAT sn1 acyl and sn2 acyl template structure
  LipdMAPS02060609152D

 28 29  0  0  0  0  0  0  0  0999 V2000
   -0.6154    2.2903    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0645    1.3323    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.7450    2.2903    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6823    1.9840    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.3994    1.8980    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.9386    1.5003    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5643    1.1214    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.3151    2.4100    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1259    1.6360    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1969    1.3323    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.3151    2.9535    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.4864   -1.1494    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0645   -0.1914    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.6160   -1.1494    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5533   -0.8430    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.2704   -0.7571    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.8096   -0.3593    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0645    0.5704    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.9970   -0.4951    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.0679   -0.1914    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2695   -1.0540    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.3994   -0.0188    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9867   -0.6463    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9867   -0.0366    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6368   -1.0216    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0524   -1.8434    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.4512   -2.3034    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4512   -2.9535    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  9  1  1      
  3 10  1  1      
  3  4  1  0      
  1  4  1  0      
  3  5  1  0      
  6 10  1  0      
  2 18  1  0      
  7  9  1  0      
  4  8  1  0      
  9 10  1  1      
  8 11  1  0      
 13 12  1  1      
 13 19  1  0      
 14 20  1  0      
 14 15  1  1      
 12 15  1  1      
 14 16  1  0      
 17 20  1  0      
 13 18  1  0      
 19 20  1  0      
 12 21  1  0      
 15 26  1  0      
 17 22  1  0      
 21 23  1  0      
 23 24  2  0      
 23 25  1  0      
 26 27  1  0      
 27 28  2  0      
M  END

ENDDATSN1SN2TEMPLATE


  my($COASn1TemplateString)=<<ENDCOASN1TEMPLATE;
COA sn1 acyl template structure
  LipdMAPS02060609152D

 52 54  0  0  0  0  0  0  0  0999 V2000
   -3.2004   -1.8639    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.8336   -2.2763    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.0088   -1.9360    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    0.0088   -1.2762    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6494   -1.2451    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4600   -2.2763    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8105   -2.2763    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6353   -1.9053    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    1.9899   -2.2763    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.9899   -3.2665    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.6393   -3.2665    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.6393   -0.1344    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    4.1239   -0.8022    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.6393   -1.4699    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    2.1411    0.8474    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    2.1411    0.0225    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.8550   -0.3896    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.8550   -1.2144    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.1411   -1.6270    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    1.4265   -1.2146    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4265   -0.3896    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    2.2694   -3.7510    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.3570   -3.7550    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.7471   -2.8288    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.3570   -4.4658    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.2694   -4.4890    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.3828   -3.9773    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0655   -2.6506    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4206   -2.6806    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.3828   -3.2482    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.9304   -4.7606    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.5560   -4.3194    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6962    1.2457    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4280    0.7807    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0013    2.5771    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7330    2.1121    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.9485    1.0920    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2471    2.1438    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5054    2.5771    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.5054    3.4470    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2471    3.8803    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0013    3.4470    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2471    1.5251    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.2442   -0.6677    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4192   -0.0902    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.0415   -0.1728    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.1239   -1.0252    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9518   -0.8601    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5119   -0.2003    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8183    3.9404    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8183    4.7606    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5337    3.5305    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  6  1  1  0      
  9  2  1  0      
  7  3  1  0      
  3  2  1  0      
  4  3  2  0      
  8  5  2  0      
  8  6  1  0      
  8  7  1  0      
 10  9  1  0      
 14 11  1  0      
 18 14  1  0      
 14 13  1  0      
 13 12  2  0      
 12 17  1  0      
 16 15  1  0      
 21 20  1  0      
 20 19  2  0      
 19 18  1  0      
 18 17  2  0      
 17 16  1  0      
 16 21  2  0      
 23 25  1  0      
 22 26  1  0      
 26 27  1  0      
  8 28  1  0      
  3 29  1  0      
 27 30  2  0      
 27 31  1  0      
 27 32  1  0      
 10 24  1  0      
 24 11  1  0      
 22 23  1  0      
 10 22  1  1      
 11 23  1  1      
 33 34  1  0      
 35 36  1  0      
 36 33  1  0      
 34 37  2  0      
 35 38  1  0      
 38 39  1  0      
 39 40  1  0      
 40 41  1  0      
 41 42  1  0      
 38 43  2  0      
  1 44  1  0      
 44 45  1  0      
 44 46  1  1      
 44 47  1  6      
 45 48  1  1      
 45 49  1  6      
 45 34  1  0      
 50 51  2  0      
 50 52  1  0      
 50 42  1  0      
M  END

ENDCOASN1TEMPLATE


  my($DIMA22TemplateString)=<<ENDDIMA22TEMPLATE;
DIMA_22_Ethyl sn1 acyl and sn2 acyl template structure
  LipdMAPS12221016212D

 43 42  0  0  0  0  0  0  0  0999 V2000
    7.7703    2.7456    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.3795    2.0686    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.7703    1.3919    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.3795    0.7151    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.7703    0.0383    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.3795   -0.6385    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.7703   -1.3154    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.3795   -1.9922    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.7703   -2.6690    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.3795   -3.3459    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.7703   -4.0226    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.3795   -4.6995    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.6414    2.0727    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    6.5546    0.7050    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    8.4867   -2.6817    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.7066   -3.3547    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    6.0770   -3.8540    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.7731    0.3142    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.7731   -0.4022    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.8486    2.4451    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.8430    3.1728    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.3795    3.3535    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.8137    4.0265    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.2493    4.6560    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.5343    4.2436    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.6508    4.6995    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.7672    4.2436    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.8836    4.6995    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.9999    4.2436    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.1165    4.6995    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2328    4.2436    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3491    4.6995    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5343    4.2436    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4180    4.6995    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.3015    4.2436    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1852    4.6995    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.0687    4.2436    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.9523    4.6995    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.8359    4.2436    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.7196    4.6995    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.6030    4.2436    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -8.4867    4.6995    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -9.3867    4.2436    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  1  0      
  5  6  1  0      
  6  7  1  0      
  7  8  1  0      
  8  9  1  0      
  9 10  1  0      
 10 11  1  0      
 11 12  1  0      
  2 13  1  1      
  4 14  1  6      
  9 15  1  1      
 10 16  1  0      
 16 17  1  0      
 14 18  1  0      
 18 19  2  0      
 13 20  1  0      
 20 21  2  0      
  1 22  1  0      
 22 23  1  0      
 23 24  1  0      
 24 25  1  0      
 25 26  1  0      
 26 27  1  0      
 27 28  1  0      
 28 29  1  0      
 29 30  1  0      
 30 31  1  0      
 31 32  1  0      
 32 33  1  0      
 33 34  1  0      
 34 35  1  0      
 35 36  1  0      
 36 37  1  0      
 37 38  1  0      
 38 39  1  0      
 39 40  1  0      
 40 41  1  0      
 41 42  1  0      
 42 43  1  0      
 M  END

ENDDIMA22TEMPLATE

  my($DIMA20TemplateString)=<<ENDDIMA20TEMPLATE;
DIMA_20_Ethyl sn1 acyl and sn2 acyl template structure
  LipdMAPS12221016212D

 41 40  0  0  0  0  0  0  0  0999 V2000
    6.8870    2.7457    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.4962    2.0687    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.8870    1.3919    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.4962    0.7151    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.8870    0.0383    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.4962   -0.6385    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.8870   -1.3154    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.4962   -1.9923    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.8870   -2.6691    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.4962   -3.3460    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.8870   -4.0227    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.4962   -4.6997    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.7581    2.0728    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.6712    0.7051    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.6035   -2.6818    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.8233   -3.3548    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.1937   -3.8542    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.8897    0.3142    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.8897   -0.4022    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.9652    2.4452    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.9596    3.1729    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    6.4962    3.3536    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.9304    4.0266    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.3660    4.6562    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.6510    4.2438    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.7674    4.6997    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.8838    4.2438    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0002    4.6997    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.1164    4.2438    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2330    4.6997    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3493    4.2438    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5345    4.6997    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4179    4.2438    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.3016    4.6997    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1852    4.2438    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.0689    4.6997    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.9524    4.2438    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.8361    4.6997    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.7197    4.2438    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.6035    4.6997    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -8.5035    4.2438    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0 
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  1  0      
  5  6  1  0      
  6  7  1  0      
  7  8  1  0      
  8  9  1  0      
  9 10  1  0      
 10 11  1  0      
 11 12  1  0      
  2 13  1  1      
  4 14  1  6      
  9 15  1  1      
 10 16  1  0      
 16 17  1  0      
 14 18  1  0      
 18 19  2  0      
 13 20  1  0      
 20 21  2  0      
  1 22  1  0      
 22 23  1  0      
 23 24  1  0      
 24 25  1  0      
 25 26  1  0      
 26 27  1  0      
 27 28  1  0      
 28 29  1  0      
 29 30  1  0      
 30 31  1  0      
 31 32  1  0      
 32 33  1  0      
 33 34  1  0      
 34 35  1  0      
 35 36  1  0      
 36 37  1  0      
 37 38  1  0      
 38 39  1  0      
 39 40  1  0  
 40 41  1  0 
M  END

ENDDIMA20TEMPLATE

  my($DIMA22MeTemplateString)=<<ENDDIMA22METEMPLATE;
DIMA22Me sn1 acyl and sn2 acyl template structure
  LipdMAPS01101111592D

 42 41  0  0  0  0  0  0  0  0999 V2000
    7.7665    2.4060    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.3758    1.7293    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.7665    1.0529    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.3758    0.3765    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.7665   -0.3000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.3758   -0.9765    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.7665   -1.6530    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.3758   -2.3295    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.7665   -3.0060    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.3758   -3.6825    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.7665   -4.3589    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.6381    1.7334    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    6.5514    0.3664    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    8.4825   -3.0187    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.7033   -3.6913    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    6.0740   -4.1904    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.7702   -0.0242    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.7702   -0.7403    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.8457    2.1056    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.8401    2.8329    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.3758    3.0136    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.8098    3.6862    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.2457    4.3154    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.5311    3.9032    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.6480    4.3589    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.7648    3.9032    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.8817    4.3589    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.9984    3.9032    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.1155    4.3589    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2322    3.9032    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3489    4.3589    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5340    3.9032    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4173    4.3589    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.3004    3.9032    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1836    4.3589    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.0667    3.9032    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.9498    4.3589    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.8330    3.9032    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.7163    4.3589    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.5992    3.9032    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -8.4825    4.3589    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -9.3725    3.9032    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  1  0      
  5  6  1  0      
  6  7  1  0      
  7  8  1  0      
  8  9  1  0      
  9 10  1  0      
 10 11  1  0      
  2 12  1  1      
  4 13  1  6      
  9 14  1  1      
 10 15  1  0      
 15 16  1  0      
 13 17  1  0      
 17 18  2  0      
 12 19  1  0      
 19 20  2  0      
  1 21  1  0      
 21 22  1  0      
 22 23  1  0      
 23 24  1  0      
 24 25  1  0      
 25 26  1  0      
 26 27  1  0      
 27 28  1  0      
 28 29  1  0      
 29 30  1  0      
 30 31  1  0      
 31 32  1  0      
 32 33  1  0      
 33 34  1  0      
 34 35  1  0      
 35 36  1  0      
 36 37  1  0      
 37 38  1  0      
 38 39  1  0      
 39 40  1  0      
 40 41  1  0      
 41 42  1  0      
 M  END

ENDDIMA22METEMPLATE

  my($DIMA20MeTemplateString)=<<ENDDIMA20METEMPLATE;
DIMA20Me sn1 acyl and sn2 acyl template structure
  LipdMAPS01101112002D

 40 39  0  0  0  0  0  0  0  0999 V2000
    6.8870    2.4072    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.4962    1.7302    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.8870    1.0535    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.4962    0.3766    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.8870   -0.3002    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.4962   -0.9770    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.8870   -1.6539    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.4962   -2.3307    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.8870   -3.0076    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.4962   -3.6845    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.8870   -4.3612    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.7581    1.7343    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.6712    0.3666    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.6035   -3.0203    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.8233   -3.6933    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.1937   -4.1926    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.8897   -0.0243    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.8897   -0.7407    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.9652    2.1067    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.9596    2.8344    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    6.4962    3.0151    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.9304    3.6882    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.3660    4.3177    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.6510    3.9053    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.7674    4.3612    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.8838    3.9053    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0002    4.3612    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.1164    3.9053    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2330    4.3612    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3493    3.9053    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5345    4.3612    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4179    3.9053    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.3016    4.3612    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1852    3.9053    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.0689    4.3612    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.9524    3.9053    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.8361    4.3612    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.7197    3.9053    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.6035    4.3612    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -8.5035    3.9053    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  1  0      
  5  6  1  0      
  6  7  1  0      
  7  8  1  0      
  8  9  1  0      
  9 10  1  0      
 10 11  1  0      
  2 12  1  1      
  4 13  1  6      
  9 14  1  1      
 10 15  1  0      
 15 16  1  0      
 13 17  1  0      
 17 18  2  0      
 12 19  1  0      
 19 20  2  0      
  1 21  1  0      
 21 22  1  0      
 22 23  1  0      
 23 24  1  0      
 24 25  1  0      
 25 26  1  0      
 26 27  1  0      
 27 28  1  0      
 28 29  1  0      
 29 30  1  0      
 30 31  1  0      
 31 32  1  0      
 32 33  1  0      
 33 34  1  0      
 34 35  1  0      
 35 36  1  0      
 36 37  1  0      
 37 38  1  0      
 38 39  1  0      
 39 40  1  0      
 M  END

ENDDIMA20METEMPLATE

  my($DIMB20MeTemplateString)=<<ENDDIMB20METEMPLATE;
DIMB20_Me sn1 acyl and sn2 acyl template structure
  LipdMAPS01101112322D

 39 38  0  0  0  0  0  0  0  0999 V2000
    6.8836    2.4060    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.4930    1.7293    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.8836    1.0530    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.4930    0.3764    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.8836   -0.3001    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.4930   -0.9765    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.8836   -1.6531    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.4930   -2.3295    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.8836   -3.0061    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.4930   -3.6827    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.8836   -4.3590    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.7553    1.7334    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.6684    0.3664    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.5997   -3.0188    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.8204   -3.6915    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.8873   -0.0243    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.8873   -0.7403    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.9627    2.1057    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.9571    2.8330    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    6.4930    3.0136    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.9270    3.6864    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.3629    4.3156    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.6482    3.9034    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.7650    4.3590    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.8819    3.9034    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.9987    4.3590    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.1154    3.9034    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2324    4.3590    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3491    3.9034    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5342    4.3590    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4172    3.9034    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.3005    4.3590    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1836    3.9034    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.0669    4.3590    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.9500    3.9034    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.8332    4.3590    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.7164    3.9034    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.5997    4.3590    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -8.4897    3.9034    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  1  0      
  5  6  1  0      
  6  7  1  0      
  7  8  1  0      
  8  9  1  0      
  9 10  1  0      
 10 11  1  0      
  2 12  1  1      
  4 13  1  6      
  9 14  1  1      
 10 15  2  0      
 13 16  1  0      
 16 17  2  0      
 12 18  1  0      
 18 19  2  0      
  1 20  1  0      
 20 21  1  0      
 21 22  1  0      
 22 23  1  0      
 23 24  1  0      
 24 25  1  0      
 25 26  1  0      
 26 27  1  0      
 27 28  1  0      
 28 29  1  0      
 29 30  1  0      
 30 31  1  0      
 31 32  1  0      
 32 33  1  0      
 33 34  1  0      
 34 35  1  0      
 35 36  1  0      
 36 37  1  0      
 37 38  1  0      
 38 39  1  0      
 M  END

ENDDIMB20METEMPLATE

  my($DIMB22MeTemplateString)=<<ENDDIMB22METEMPLATE;
DIMB22_Me sn1 acyl and sn2 acyl template structure
  LipdMAPS01101112422D

 41 40  0  0  0  0  0  0  0  0999 V2000
    7.5829    2.3491    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.2014    1.6884    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.5829    1.0280    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.2014    0.3676    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.5829   -0.2929    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.2014   -0.9534    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.5829   -1.6139    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.2014   -2.2744    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.5829   -2.9349    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.2014   -3.5954    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.5829   -4.2559    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.4811    1.6924    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    6.3965    0.3577    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    8.2820   -2.9473    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.5448   -3.6040    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.6338   -0.0236    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.6338   -0.7228    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.7075    2.0559    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.7021    2.7659    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.2014    2.9423    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.6251    3.5991    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.0744    4.2134    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.3767    3.8109    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.5145    4.2559    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.6522    3.8109    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.7899    4.2559    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.9275    3.8109    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0655    4.2559    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2031    3.8109    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3406    4.2559    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5213    3.8109    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3838    4.2559    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2460    3.8109    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1084    4.2559    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.9705    3.8109    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.8328    4.2559    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.6951    3.8109    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.5575    4.2559    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.4195    3.8109    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -8.2820    4.2559    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -9.1720    3.8109    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  1  0      
  5  6  1  0      
  6  7  1  0      
  7  8  1  0      
  8  9  1  0      
  9 10  1  0      
 10 11  1  0      
  2 12  1  1      
  4 13  1  6      
  9 14  1  1      
 10 15  2  0      
 13 16  1  0      
 16 17  2  0      
 12 18  1  0      
 18 19  2  0      
  1 20  1  0      
 20 21  1  0      
 21 22  1  0      
 22 23  1  0      
 23 24  1  0      
 24 25  1  0      
 25 26  1  0      
 26 27  1  0      
 27 28  1  0      
 28 29  1  0      
 29 30  1  0      
 30 31  1  0      
 31 32  1  0      
 32 33  1  0      
 33 34  1  0      
 34 35  1  0      
 35 36  1  0      
 36 37  1  0      
 37 38  1  0      
 38 39  1  0      
 39 40  1  0      
 40 41  1  0      
 M  END

ENDDIMB22METEMPLATE

  my($DIMB20TemplateString)=<<ENDDIMB20TEMPLATE;
DIMB20 template structure
  LipdMAPS01101112442D

 40 39  0  0  0  0  0  0  0  0999 V2000
    6.8865    2.7455    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.4957    2.0685    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.8865    1.3918    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.4957    0.7150    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.8865    0.0383    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.4957   -0.6385    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.8865   -1.3153    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.4957   -1.9922    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.8865   -2.6689    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.4957   -3.3458    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.8865   -4.0224    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.4957   -4.6994    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.7577    2.0726    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.6708    0.7050    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.6029   -2.6816    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.8229   -3.3546    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.8893    0.3142    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.8893   -0.4022    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.9648    2.4450    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.9592    3.1727    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    6.4957    3.3534    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.9299    4.0263    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.3655    4.6559    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.6506    4.2435    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.7670    4.6994    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.8835    4.2435    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0000    4.6994    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.1162    4.2435    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2329    4.6994    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3493    4.2435    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5345    4.6994    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4178    4.2435    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.3014    4.6994    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1850    4.2435    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.0686    4.6994    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.9520    4.2435    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.8357    4.6994    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.7192    4.2435    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.6029    4.6994    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -8.5029    4.2435    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  1  0      
  5  6  1  0      
  6  7  1  0      
  7  8  1  0      
  8  9  1  0      
  9 10  1  0      
 10 11  1  0      
 11 12  1  0      
  2 13  1  1      
  4 14  1  6      
  9 15  1  1      
 10 16  2  0      
 14 17  1  0      
 17 18  2  0      
 13 19  1  0      
 19 20  2  0      
  1 21  1  0      
 21 22  1  0      
 22 23  1  0      
 23 24  1  0      
 24 25  1  0      
 25 26  1  0      
 26 27  1  0      
 27 28  1  0      
 28 29  1  0      
 29 30  1  0      
 30 31  1  0      
 31 32  1  0      
 32 33  1  0      
 33 34  1  0      
 34 35  1  0      
 35 36  1  0      
 36 37  1  0      
 37 38  1  0      
 38 39  1  0      
 39 40  1  0      
 M  END

ENDDIMB20TEMPLATE

  my($DIMB22TemplateString)=<<ENDDIMB22TEMPLATE;
DIMB22 template structure
  LipdMAPS01101112452D

 42 41  0  0  0  0  0  0  0  0999 V2000
    7.7665    2.7442    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.3758    2.0676    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.7665    1.3912    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.3758    0.7147    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.7665    0.0383    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.3758   -0.6382    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.7665   -1.3147    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.3758   -1.9912    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.7665   -2.6677    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.3758   -3.3442    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.7665   -4.0206    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.3758   -4.6972    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.6381    2.0717    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    6.5514    0.7047    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    8.4825   -2.6804    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.7033   -3.3530    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.7702    0.3140    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.7702   -0.4020    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.8457    2.4439    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.8401    3.1712    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.3758    3.3518    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.8098    4.0245    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.2457    4.6537    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.5311    4.2415    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.6480    4.6972    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.7648    4.2415    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.8817    4.6972    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.9984    4.2415    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.1155    4.6972    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2322    4.2415    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3489    4.6972    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5340    4.2415    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4173    4.6972    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.3004    4.2415    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1836    4.6972    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.0667    4.2415    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.9498    4.6972    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.8330    4.2415    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.7163    4.6972    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.5992    4.2415    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -8.4825    4.6972    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -9.3825    4.2415    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  1  0      
  5  6  1  0      
  6  7  1  0      
  7  8  1  0      
  8  9  1  0      
  9 10  1  0      
 10 11  1  0      
 11 12  1  0      
  2 13  1  1      
  4 14  1  6      
  9 15  1  1      
 10 16  2  0      
 14 17  1  0      
 17 18  2  0      
 13 19  1  0      
 19 20  2  0      
  1 21  1  0      
 21 22  1  0      
 22 23  1  0      
 23 24  1  0      
 24 25  1  0      
 25 26  1  0      
 26 27  1  0      
 27 28  1  0      
 28 29  1  0      
 29 30  1  0      
 30 31  1  0      
 31 32  1  0      
 32 33  1  0      
 33 34  1  0      
 34 35  1  0      
 35 36  1  0      
 36 37  1  0      
 37 38  1  0      
 38 39  1  0      
 39 40  1  0      
 40 41  1  0      
 41 42  1  0      
 M  END

ENDDIMB22TEMPLATE

  my($PIM2Sn1Sn2TemplateString)=<<ENDPIM2TEMPLATE;
PIM2 sn1 acyl and sn2 acyl template structure
  LipdMAPS01111109592D

 50 52  0  0  0  0  0  0  0  0999 V2000
   -3.0504    1.0947    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.7627    1.5047    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.4754    1.0947    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.1874    1.5047    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.1874    2.3282    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6385    0.3824    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -3.4621    0.3824    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.8998    1.0947    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.3377    1.5060    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6254    1.0947    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6116    1.3830    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9745    0.7538    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6116    2.1332    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.2084   -0.0387    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.2084   -0.8625    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.9209    0.3729    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.4475    0.8166    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.0431   -0.2192    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.1447    0.8166    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.2938    0.4854    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.7846    0.0692    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.8954    0.1091    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.7374   -0.2192    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.3273    0.5809    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.3330   -2.9707    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.7374   -1.9348    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.6359   -2.9707    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.4868   -2.6395    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.3434   -2.5465    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.8452   -2.1164    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.8852   -3.0869    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.7374   -1.1111    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.1709   -3.1000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.8852   -2.2631    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.0431   -1.9348    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.8487   -3.7671    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.3253   -0.6354    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.8895    3.0776    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.2938    2.0418    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.1922    3.0776    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.0431    2.7464    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.8998    2.6535    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.4015    2.2235    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.4416    3.1938    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.7344    1.0591    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.7274    3.2070    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.4416    2.3701    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.5995    2.0418    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.4368    3.7671    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.1029    0.9705    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  2  0      
  4  8  1  0      
  1  7  1  6      
  1  6  1  1      
  9  1  1  0      
 10  9  1  0      
 11 12  1  0      
 11 13  2  0      
 14 15  2  0      
 14 16  1  0      
 14  7  1  0      
 11 10  1  0      
 17 18  1  0      
 18 22  1  1      
 19 23  1  1      
 19 20  1  0      
 17 20  1  0      
 18 21  1  0      
 23 22  1  1      
 20 45  1  0      
 17 24  1  0      
 22 37  1  0      
 23 32  1  0      
 25 26  1  0      
 26 34  1  1      
 27 35  1  1      
 27 28  1  0      
 25 28  1  0      
 27 29  1  0      
 30 35  1  0      
 31 34  1  0      
 26 32  1  0      
 28 33  1  0      
 34 35  1  1      
 33 36  1  0      
 38 39  1  0      
 39 47  1  1      
 40 48  1  1      
 40 41  1  0      
 38 41  1  0      
 40 42  1  0      
 43 48  1  0      
 44 47  1  0      
 39 45  1  0      
 41 46  1  0      
 47 48  1  1      
 46 49  1  0      
 11 50  1  0      
 50 19  1  0      
M  END

ENDPIM2TEMPLATE

  my($PIM1Sn1Sn2TemplateString)=<<ENDPIM1TEMPLATE;
PIM1 sn1 acyl and sn2 acyl structure
  LipdMAPS01111109592D

 39 40  0  0  0  0  0  0  0  0999 V2000
   -2.7722    1.8141    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.4845    2.2242    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.1972    1.8141    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.9093    2.2242    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.9093    3.0477    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.3603    1.1018    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1839    1.1018    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.6216    1.8141    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0595    2.2255    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3472    1.8141    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3334    2.1024    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6963    1.4732    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3334    2.8527    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.9302    0.6807    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.9302   -0.1430    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.6427    1.0923    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.7257    1.5360    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.3213    0.5002    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4229    1.5360    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.5720    1.2048    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.0627    0.7886    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.1736    0.8285    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0156    0.5002    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.6055    1.3003    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.6112   -2.2513    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.0156   -1.2154    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.9141   -2.2513    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.7650   -1.9201    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.6216   -1.8271    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.1233   -1.3970    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.1634   -2.3675    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.0156   -0.3916    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.4491   -2.3806    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.1634   -1.5437    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.3213   -1.2154    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.1269   -3.0477    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.6035    0.0840    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.0126    1.7785    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.3811    1.6899    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  2  0      
  4  8  1  0      
  1  7  1  6      
  1  6  1  1      
  9  1  1  0      
 10  9  1  0      
 11 12  1  0      
 11 13  2  0      
 14 15  2  0      
 14 16  1  0      
 14  7  1  0      
 11 10  1  0      
 17 18  1  0      
 18 22  1  1      
 19 23  1  1      
 19 20  1  0      
 17 20  1  0      
 18 21  1  0      
 23 22  1  1      
 20 38  1  0      
 17 24  1  0      
 22 37  1  0      
 23 32  1  0      
 25 26  1  0      
 26 34  1  1      
 27 35  1  1      
 27 28  1  0      
 25 28  1  0      
 27 29  1  0      
 30 35  1  0      
 31 34  1  0      
 26 32  1  0      
 28 33  1  0      
 34 35  1  1      
 33 36  1  0      
 11 39  1  0      
 39 19  1  0      
M  END

ENDPIM1TEMPLATE

  my($PIM4Sn1Sn2TemplateString)=<<ENDPIM4TEMPLATE;
PIM4 sn1 acyl and sn2 acyl template structure
  LipdMAPS01111109592D

 72 76  0  0  0  0  0  0  0  0999 V2000
   -5.4428   -1.0932    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.1551   -0.6832    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.8676   -1.0932    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -7.5797   -0.6832    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.5797    0.1403    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.0309   -1.8055    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -5.8545   -1.8055    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -8.2921   -1.0932    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.7301   -0.6819    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.0177   -1.0932    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0039   -0.8050    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -3.3668   -1.4341    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0039   -0.0547    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -6.6006   -2.2266    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.6006   -3.0503    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -7.3131   -1.8150    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0550   -1.3713    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.6506   -2.4071    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2477   -1.3713    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0985   -1.7025    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.3921   -2.1187    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.5030   -2.0788    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6549   -2.4071    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.9348   -1.6070    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0593   -5.1586    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6549   -4.1227    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.2434   -5.1586    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0943   -4.8274    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.9510   -4.7344    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.4527   -4.3043    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.4928   -5.2748    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6549   -3.2989    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.7784   -5.2879    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4928   -4.4510    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.6506   -4.1227    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4562   -5.9549    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.9328   -2.8232    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.4146    0.5802    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1810   -0.4554    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.7173    0.5802    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5681    0.2492    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.4249    0.1563    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.9266   -0.2737    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.9666    0.6964    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.3421   -1.1288    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.2524    0.7096    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.9666   -0.1271    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.1245   -0.4554    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.8481    2.9805    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.2524    1.9447    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.1508    2.9805    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.0017    2.6494    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.8584    2.5564    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.3601    2.1263    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.4001    3.0967    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.2524    1.2859    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.6859    3.1099    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.4001    2.2730    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.5581    1.9447    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.6858    3.7068    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2895   -1.2174    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.2819    5.4349    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.6862    4.3991    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.5845    5.4349    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.4354    5.1038    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.2921    5.0108    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.7938    4.5807    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.8339    5.5511    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.1196    5.5643    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.8339    4.7274    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.9918    4.3991    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.8514    5.9549    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  2  0      
  4  8  1  0      
  1  7  1  6      
  1  6  1  1      
  9  1  1  0      
 10  9  1  0      
 11 12  1  0      
 11 13  2  0      
 14 15  2  0      
 14 16  1  0      
 14  7  1  0      
 11 10  1  0      
 17 18  1  0      
 18 22  1  1      
 19 23  1  1      
 19 20  1  0      
 17 20  1  0      
 18 21  1  0      
 23 22  1  1      
 20 45  1  0      
 17 24  1  0      
 22 37  1  0      
 23 32  1  0      
 25 26  1  0      
 26 34  1  1      
 27 35  1  1      
 27 28  1  0      
 25 28  1  0      
 27 29  1  0      
 30 35  1  0      
 31 34  1  0      
 26 32  1  0      
 28 33  1  0      
 34 35  1  1      
 33 36  1  0      
 38 39  1  0      
 39 47  1  1      
 40 48  1  1      
 40 41  1  0      
 38 41  1  0      
 40 42  1  0      
 43 48  1  0      
 44 47  1  0      
 39 45  1  0      
 41 46  1  0      
 47 48  1  1      
 46 56  1  0      
 49 50  1  0      
 50 58  1  1      
 51 59  1  1      
 51 52  1  0      
 49 52  1  0      
 51 53  1  0      
 54 59  1  0      
 55 58  1  0      
 50 56  1  0      
 52 57  1  0      
 58 59  1  1      
 57 60  1  0      
 11 61  1  0      
 61 19  1  0      
 62 63  1  0      
 63 70  1  1      
 64 71  1  1      
 64 65  1  0      
 62 65  1  0      
 64 66  1  0      
 67 71  1  0      
 68 70  1  0      
 65 69  1  0      
 70 71  1  1      
 69 72  1  0      
 60 63  1  0      
M  END

ENDPIM4TEMPLATE

  my($PIM5Sn1Sn2TemplateString)=<<ENDPIM5TEMPLATE;
PIM5 sn1 acyl and sn2 acyl structure
  LipdMAPS01111109592D

 83 88  0  0  0  0  0  0  0  0999 V2000
   -6.6596   -2.2482    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.3719   -1.8382    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -8.0845   -2.2482    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -8.7966   -1.8382    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -8.7966   -1.0147    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -6.2477   -2.9605    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -7.0713   -2.9605    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -9.5090   -2.2482    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.9469   -1.8369    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.2346   -2.2482    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.2208   -1.9600    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -4.5837   -2.5891    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.2208   -1.2097    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -7.8175   -3.3816    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.8175   -4.2053    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -8.5300   -2.9700    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1618   -2.5263    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4338   -3.5621    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4646   -2.5263    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3154   -2.8575    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.1752   -3.2737    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7139   -3.2338    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8718   -3.5621    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7180   -2.7620    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2762   -6.3136    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8718   -5.2777    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0266   -6.3136    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1225   -5.9824    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.7341   -5.8894    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.2358   -5.4593    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7241   -6.4298    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8718   -4.4539    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.5616   -6.4429    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7241   -5.6060    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4338   -5.2777    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.2394   -7.1099    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2840   -3.9782    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8023   -0.5748    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3979   -1.6104    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5004   -0.5748    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3513   -0.9058    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.2080   -0.9987    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.7097   -1.4287    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2502   -0.4586    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8748   -2.2838    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.0355   -0.4454    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2502   -1.2821    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.9077   -1.6104    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.6312    1.8255    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.0355    0.7897    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.9339    1.8255    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.7848    1.4944    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.6415    1.4014    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.1432    0.9713    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.1832    1.9417    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.0355    0.1309    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.4690    1.9549    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.1832    1.1180    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.3412    0.7897    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.4689    2.5518    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.5064   -2.3724    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.0650    4.2799    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.4693    3.2441    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.3677    4.2799    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.2186    3.9488    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.0752    3.8558    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    6.5770    3.4257    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.6170    4.3961    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.9028    4.4093    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.6170    3.5724    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.7750    3.2441    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.9027    4.9236    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    6.4987    6.5899    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.9030    5.5541    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.8014    6.5899    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.6523    6.2588    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    9.5090    6.1658    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    9.0107    5.7357    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.0507    6.7061    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    8.3365    6.7193    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.0507    5.8824    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.2087    5.5541    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.0683    7.1099    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  2  0      
  4  8  1  0      
  1  7  1  6      
  1  6  1  1      
  9  1  1  0      
 10  9  1  0      
 11 12  1  0      
 11 13  2  0      
 14 15  2  0      
 14 16  1  0      
 14  7  1  0      
 11 10  1  0      
 17 18  1  0      
 18 22  1  1      
 19 23  1  1      
 19 20  1  0      
 17 20  1  0      
 18 21  1  0      
 23 22  1  1      
 20 45  1  0      
 17 24  1  0      
 22 37  1  0      
 23 32  1  0      
 25 26  1  0      
 26 34  1  1      
 27 35  1  1      
 27 28  1  0      
 25 28  1  0      
 27 29  1  0      
 30 35  1  0      
 31 34  1  0      
 26 32  1  0      
 28 33  1  0      
 34 35  1  1      
 33 36  1  0      
 38 39  1  0      
 39 47  1  1      
 40 48  1  1      
 40 41  1  0      
 38 41  1  0      
 40 42  1  0      
 43 48  1  0      
 44 47  1  0      
 39 45  1  0      
 41 46  1  0      
 47 48  1  1      
 46 56  1  0      
 49 50  1  0      
 50 58  1  1      
 51 59  1  1      
 51 52  1  0      
 49 52  1  0      
 51 53  1  0      
 54 59  1  0      
 55 58  1  0      
 50 56  1  0      
 52 57  1  0      
 58 59  1  1      
 57 60  1  0      
 11 61  1  0      
 61 19  1  0      
 62 63  1  0      
 63 70  1  1      
 64 71  1  1      
 64 65  1  0      
 62 65  1  0      
 64 66  1  0      
 67 71  1  0      
 68 70  1  0      
 65 69  1  0      
 70 71  1  1      
 69 72  1  0      
 60 63  1  0      
 73 74  1  0      
 74 81  1  1      
 75 82  1  1      
 75 76  1  0      
 73 76  1  0      
 75 77  1  0      
 78 82  1  0      
 79 81  1  0      
 76 80  1  0      
 81 82  1  1      
 80 83  1  0      
 74 72  1  0      
M  END

ENDPIM5TEMPLATE

  my($PIM6Sn1Sn2TemplateString)=<<ENDPIM6TEMPLATE;
PIM6 sn1 acyl and sn2 acyl structure
  LipdMAPS01111109592D

 94100  0  0  0  0  0  0  0  0999 V2000
   -7.8765   -3.4009    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -8.5888   -2.9908    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -9.3014   -3.4009    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  -10.0135   -2.9908    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  -10.0135   -2.1673    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -7.4646   -4.1131    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -8.2882   -4.1131    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  -10.7259   -3.4009    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.1638   -2.9895    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.4514   -3.4009    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.4376   -3.1126    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -5.8005   -3.7417    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.4376   -2.3623    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -9.0344   -4.5342    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -9.0344   -5.3579    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -9.7469   -4.1226    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3787   -3.6789    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7831   -4.7147    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.6815   -3.6789    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5323   -4.0101    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0416   -4.4263    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9308   -4.3864    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0887   -4.7147    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4989   -3.9146    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4931   -7.4662    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0887   -6.4303    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1903   -7.4662    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3394   -7.1350    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5173   -7.0420    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.0190   -6.6119    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9410   -7.5824    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0887   -5.6065    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6553   -7.5955    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9410   -6.7586    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7831   -6.4303    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9775   -8.2625    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5009   -5.1308    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0192   -1.7274    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6148   -2.7630    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.2836   -1.7274    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8656   -2.0584    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.9912   -2.1513    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.4929   -2.5813    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4671   -1.6112    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0917   -3.4364    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1813   -1.5980    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4671   -2.4347    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3092   -2.7630    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4144    0.6729    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1813   -0.3629    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.7170    0.6729    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5679    0.3418    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.4246    0.2488    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.9263   -0.1813    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.9664    0.7891    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1813   -1.0217    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.2521    0.8023    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.9664   -0.0346    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.1243   -0.3629    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.2521    1.3991    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.7233   -3.5250    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.8481    3.1273    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.2524    2.0915    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.1508    3.1273    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.0017    2.7962    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.8584    2.7032    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.3601    2.2731    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.4001    3.2435    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.6859    3.2567    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.4001    2.4198    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.5581    2.0915    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.6858    3.7710    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.2819    5.4373    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.6862    4.4015    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.5845    5.4373    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.4354    5.1062    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.2921    5.0132    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.7938    4.5831    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.8339    5.5535    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.1196    5.5667    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.8339    4.7298    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.9918    4.4015    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.1196    6.0604    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.7156    7.7425    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.1199    6.7067    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.0183    7.7425    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.8692    7.4114    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.7259    7.3184    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   10.2276    6.8883    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    8.2676    7.8587    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    9.5534    7.8719    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.2676    7.0350    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    9.4256    6.7067    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    9.2852    8.2625    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  2  0      
  4  8  1  0      
  1  7  1  6      
  1  6  1  1      
  9  1  1  0      
 10  9  1  0      
 11 12  1  0      
 11 13  2  0      
 14 15  2  0      
 14 16  1  0      
 14  7  1  0      
 11 10  1  0      
 17 18  1  0      
 18 22  1  1      
 19 23  1  1      
 19 20  1  0      
 17 20  1  0      
 18 21  1  0      
 23 22  1  1      
 20 45  1  0      
 17 24  1  0      
 22 37  1  0      
 23 32  1  0      
 25 26  1  0      
 26 34  1  1      
 27 35  1  1      
 27 28  1  0      
 25 28  1  0      
 27 29  1  0      
 30 35  1  0      
 31 34  1  0      
 26 32  1  0      
 28 33  1  0      
 34 35  1  1      
 33 36  1  0      
 38 39  1  0      
 39 47  1  1      
 40 48  1  1      
 40 41  1  0      
 38 41  1  0      
 40 42  1  0      
 43 48  1  0      
 44 47  1  0      
 39 45  1  0      
 41 46  1  0      
 47 48  1  1      
 46 56  1  0      
 49 50  1  0      
 50 58  1  1      
 51 59  1  1      
 51 52  1  0      
 49 52  1  0      
 51 53  1  0      
 54 59  1  0      
 55 58  1  0      
 50 56  1  0      
 52 57  1  0      
 58 59  1  1      
 57 60  1  0      
 11 61  1  0      
 61 19  1  0      
 62 63  1  0      
 63 70  1  1      
 64 71  1  1      
 64 65  1  0      
 62 65  1  0      
 64 66  1  0      
 67 71  1  0      
 68 70  1  0      
 65 69  1  0      
 70 71  1  1      
 69 72  1  0      
 60 63  1  0      
 73 74  1  0      
 74 81  1  1      
 75 82  1  1      
 75 76  1  0      
 73 76  1  0      
 75 77  1  0      
 78 82  1  0      
 79 81  1  0      
 76 80  1  0      
 81 82  1  1      
 80 83  1  0      
 74 72  1  0      
 84 85  1  0      
 85 92  1  1      
 86 93  1  1      
 86 87  1  0      
 84 87  1  0      
 86 88  1  0      
 89 93  1  0      
 90 92  1  0      
 87 91  1  0      
 92 93  1  1      
 91 94  1  0      
 83 85  1  0      
M  END

ENDPIM6TEMPLATE


  my($PIM3LYSOSn1TemplateString)=<<ENDPIM3LYSOTEMPLATE;
PIM3.mol
  ChemDraw01111120592D

 58 61  0  0  0  0  0  0  0  0999 V2000
   -4.2671   -0.2270    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.9794    0.1831    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.6920   -0.2270    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -6.4041    0.1831    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.4041    1.0066    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.8552   -0.9393    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -4.6788   -0.9393    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -7.1165   -0.2270    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.5544    0.1844    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.8421   -0.2270    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8283    0.0613    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1912   -0.5679    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8283    0.8116    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.2307   -0.5051    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.8263   -1.5409    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0721   -0.5051    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0771   -0.8363    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.5677   -1.2525    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.6786   -1.2126    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5207   -1.5409    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.1105   -0.7408    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.1163   -4.2923    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.5207   -3.2565    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.4191   -4.2923    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.2700   -3.9611    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.1266   -3.8681    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.6283   -3.4381    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.6684   -4.4085    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.5207   -2.4327    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.9541   -4.4216    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.6684   -3.5847    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.8263   -3.2565    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.6319   -5.0886    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.1085   -1.9570    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.6727    1.7559    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.0771    0.7202    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.9754    1.7559    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.8263    1.4248    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.6830    1.3319    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.1847    0.9019    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.2248    1.8721    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.5177   -0.2626    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.5105    1.8853    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.2248    1.0485    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.3827    0.7202    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.1062    4.5686    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.5105    3.5328    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.4089    4.5686    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.2598    4.2375    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.1165    4.1445    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    6.6182    3.7144    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.6582    4.6848    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.5105    2.7091    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.9440    4.6980    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.6582    3.8611    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.8162    3.5328    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.6758    5.0886    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1139   -0.3512    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  2  0      
  4  8  1  0      
  1  7  1  6      
  1  6  1  1      
  9  1  1  0      
 10  9  1  0      
 11 12  1  0      
 11 13  2  0      
 11 10  1  0      
 14 15  1  0      
 15 19  1  1      
 16 20  1  1      
 16 17  1  0      
 14 17  1  0      
 15 18  1  0      
 20 19  1  1      
 17 42  1  0      
 14 21  1  0      
 19 34  1  0      
 20 29  1  0      
 22 23  1  0      
 23 31  1  1      
 24 32  1  1      
 24 25  1  0      
 22 25  1  0      
 24 26  1  0      
 27 32  1  0      
 28 31  1  0      
 23 29  1  0      
 25 30  1  0      
 31 32  1  1      
 30 33  1  0      
 35 36  1  0      
 36 44  1  1      
 37 45  1  1      
 37 38  1  0      
 35 38  1  0      
 37 39  1  0      
 40 45  1  0      
 41 44  1  0      
 36 42  1  0      
 38 43  1  0      
 44 45  1  1      
 43 53  1  0      
 46 47  1  0      
 47 55  1  1      
 48 56  1  1      
 48 49  1  0      
 46 49  1  0      
 48 50  1  0      
 51 56  1  0      
 52 55  1  0      
 47 53  1  0      
 49 54  1  0      
 55 56  1  1      
 54 57  1  0      
 11 58  1  0      
 58 16  1  0      
M  END

ENDPIM3LYSOTEMPLATE

  my($PIM1LYSOSn1TemplateString)=<<ENDPIM1LYSOTEMPLATE;
LYSOPIM1.mol
  ChemDraw01111121372D

 36 37  0  0  0  0  0  0  0  0999 V2000
   -2.7725    1.8143    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.4848    2.2244    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.1976    1.8143    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.9097    2.2244    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.9097    3.0480    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.3605    1.1019    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1842    1.1019    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.6221    1.8143    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0597    2.2257    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3473    1.8143    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3334    2.1026    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6964    1.4733    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3334    2.8530    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.7260    1.5361    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.3217    0.5002    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4230    1.5361    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.5722    1.2049    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.0632    0.7887    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.1739    0.8286    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0158    0.5002    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.6059    1.3004    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.6114   -2.2515    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.0158   -1.2155    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.9145   -2.2515    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.7653   -1.9203    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.6221   -1.8273    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.1238   -1.3971    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.1637   -2.3677    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.0158   -0.3916    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.4495   -2.3808    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.1637   -1.5438    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.3217   -1.2155    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.1273   -3.0480    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.6038    0.0840    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.0129    1.7787    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.3811    1.6901    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  2  0      
  4  8  1  0      
  1  7  1  6      
  1  6  1  1      
  9  1  1  0      
 10  9  1  0      
 11 12  1  0      
 11 13  2  0      
 11 10  1  0      
 14 15  1  0      
 15 19  1  1      
 16 20  1  1      
 16 17  1  0      
 14 17  1  0      
 15 18  1  0      
 20 19  1  1      
 17 35  1  0      
 14 21  1  0      
 19 34  1  0      
 20 29  1  0      
 22 23  1  0      
 23 31  1  1      
 24 32  1  1      
 24 25  1  0      
 22 25  1  0      
 24 26  1  0      
 27 32  1  0      
 28 31  1  0      
 23 29  1  0      
 25 30  1  0      
 31 32  1  1      
 30 33  1  0      
 11 36  1  0      
 36 16  1  0      
M  END

ENDPIM1LYSOTEMPLATE

  my($PIM2LYSOSn1TemplateString)=<<ENDPIM2LYSOTEMPLATE;
LYSOPIM2.mol
  ChemDraw01111121442D

 47 49  0  0  0  0  0  0  0  0999 V2000
   -3.0503    1.0947    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.7626    1.5047    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.4753    1.0947    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.1872    1.5047    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.1872    2.3281    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6384    0.3824    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -3.4620    0.3824    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.8996    1.0947    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.3376    1.5060    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6253    1.0947    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6116    1.3830    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9745    0.7538    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6116    2.1331    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.4474    0.8166    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.0430   -0.2192    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.1447    0.8166    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.2937    0.4854    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.7845    0.0692    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.8953    0.1091    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.7373   -0.2192    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.3272    0.5809    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.3329   -2.9706    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.7373   -1.9347    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.6358   -2.9706    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.4867   -2.6394    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.3432   -2.5464    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.8450   -2.1163    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.8851   -3.0868    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.7373   -1.1111    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.1708   -3.0999    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.8851   -2.2630    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.0430   -1.9347    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.8486   -3.7670    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.3252   -0.6354    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.8894    3.0775    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.2937    2.0417    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.1920    3.0775    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.0430    2.7463    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.8996    2.6534    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.4013    2.2234    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.4415    3.1937    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.7343    1.0591    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.7273    3.2069    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.4415    2.3700    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.5994    2.0417    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.4367    3.7670    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.1029    0.9705    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  2  0      
  4  8  1  0      
  1  7  1  6      
  1  6  1  1      
  9  1  1  0      
 10  9  1  0      
 11 12  1  0      
 11 13  2  0      
 11 10  1  0      
 14 15  1  0      
 15 19  1  1      
 16 20  1  1      
 16 17  1  0      
 14 17  1  0      
 15 18  1  0      
 20 19  1  1      
 17 42  1  0      
 14 21  1  0      
 19 34  1  0      
 20 29  1  0      
 22 23  1  0      
 23 31  1  1      
 24 32  1  1      
 24 25  1  0      
 22 25  1  0      
 24 26  1  0      
 27 32  1  0      
 28 31  1  0      
 23 29  1  0      
 25 30  1  0      
 31 32  1  1      
 30 33  1  0      
 35 36  1  0      
 36 44  1  1      
 37 45  1  1      
 37 38  1  0      
 35 38  1  0      
 37 39  1  0      
 40 45  1  0      
 41 44  1  0      
 36 42  1  0      
 38 43  1  0      
 44 45  1  1      
 43 46  1  0      
 11 47  1  0      
 47 16  1  0      
M  END

ENDPIM2LYSOTEMPLATE

  my($PIM4LYSOSn1TemplateString)=<<ENDPIM4LYSOTEMPLATE;
LYSOPIM4.mol
  ChemDraw01111121472D

 69 73  0  0  0  0  0  0  0  0999 V2000
   -5.4427   -1.0932    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.1550   -0.6832    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.8675   -1.0932    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -7.5796   -0.6832    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.5796    0.1403    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.0308   -1.8055    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -5.8544   -1.8055    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -8.2920   -1.0932    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.7300   -0.6819    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.0176   -1.0932    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0039   -0.8050    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -3.3668   -1.4341    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0039   -0.0547    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.0550   -1.3713    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.6506   -2.4071    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2477   -1.3713    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0985   -1.7025    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.3921   -2.1187    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.5030   -2.0788    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6549   -2.4071    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.9348   -1.6070    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0593   -5.1585    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6549   -4.1226    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.2434   -5.1585    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0943   -4.8273    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.9510   -4.7343    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.4527   -4.3042    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.4928   -5.2747    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6549   -3.2989    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.7784   -5.2878    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4928   -4.4509    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.6506   -4.1226    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4562   -5.9548    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.9328   -2.8232    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.4146    0.5802    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1810   -0.4554    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.7173    0.5802    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5681    0.2492    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.4249    0.1563    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.9266   -0.2737    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.9666    0.6964    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.3421   -1.1288    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.2524    0.7096    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.9666   -0.1271    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.1245   -0.4554    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.8481    2.9805    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.2524    1.9447    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.1507    2.9805    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.0016    2.6494    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.8583    2.5564    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.3600    2.1263    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.4001    3.0967    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.2524    1.2859    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.6858    3.1099    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.4001    2.2730    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.5580    1.9447    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.6857    3.7068    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2895   -1.2174    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.2818    5.4348    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.6861    4.3990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.5844    5.4348    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.4353    5.1037    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.2920    5.0107    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.7937    4.5806    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.8338    5.5510    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.1195    5.5642    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.8338    4.7273    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.9917    4.3990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.8513    5.9548    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  2  0      
  4  8  1  0      
  1  7  1  6      
  1  6  1  1      
  9  1  1  0      
 10  9  1  0      
 11 12  1  0      
 11 13  2  0      
 11 10  1  0      
 14 15  1  0      
 15 19  1  1      
 16 20  1  1      
 16 17  1  0      
 14 17  1  0      
 15 18  1  0      
 20 19  1  1      
 17 42  1  0      
 14 21  1  0      
 19 34  1  0      
 20 29  1  0      
 22 23  1  0      
 23 31  1  1      
 24 32  1  1      
 24 25  1  0      
 22 25  1  0      
 24 26  1  0      
 27 32  1  0      
 28 31  1  0      
 23 29  1  0      
 25 30  1  0      
 31 32  1  1      
 30 33  1  0      
 35 36  1  0      
 36 44  1  1      
 37 45  1  1      
 37 38  1  0      
 35 38  1  0      
 37 39  1  0      
 40 45  1  0      
 41 44  1  0      
 36 42  1  0      
 38 43  1  0      
 44 45  1  1      
 43 53  1  0      
 46 47  1  0      
 47 55  1  1      
 48 56  1  1      
 48 49  1  0      
 46 49  1  0      
 48 50  1  0      
 51 56  1  0      
 52 55  1  0      
 47 53  1  0      
 49 54  1  0      
 55 56  1  1      
 54 57  1  0      
 11 58  1  0      
 58 16  1  0      
 59 60  1  0      
 60 67  1  1      
 61 68  1  1      
 61 62  1  0      
 59 62  1  0      
 61 63  1  0      
 64 68  1  0      
 65 67  1  0      
 62 66  1  0      
 67 68  1  1      
 66 69  1  0      
 57 60  1  0      
M  END

ENDPIM4LYSOTEMPLATE

  my($PIM5LYSOSn1TemplateString)=<<ENDPIM5LYSOTEMPLATE;
LYSOPIM5.mol
  ChemDraw01111121502D

 80 85  0  0  0  0  0  0  0  0999 V2000
   -6.6599   -2.2483    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.3723   -1.8383    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -8.0849   -2.2483    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -8.7970   -1.8383    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -8.7970   -1.0147    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -6.2480   -2.9606    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -7.0716   -2.9606    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -9.5095   -2.2483    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.9472   -1.8370    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.2349   -2.2483    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.2210   -1.9601    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -4.5839   -2.5892    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.2210   -1.2098    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1618   -2.5264    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4338   -3.5623    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4647   -2.5264    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3155   -2.8576    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.1753   -3.2739    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7139   -3.2340    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8719   -3.5623    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7180   -2.7621    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2763   -6.3139    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8719   -5.2780    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0266   -6.3139    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1225   -5.9827    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.7342   -5.8897    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.2359   -5.4596    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7241   -6.4301    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8719   -4.4541    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.5616   -6.4432    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7241   -5.6063    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4338   -5.2780    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.2394   -7.1102    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2840   -3.9784    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8023   -0.5748    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3980   -1.6105    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5005   -0.5748    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3513   -0.9058    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.2081   -0.9987    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.7098   -1.4288    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2502   -0.4586    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8748   -2.2839    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.0356   -0.4454    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2502   -1.2822    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.9077   -1.6105    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.6313    1.8256    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.0356    0.7897    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.9341    1.8256    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.7849    1.4945    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.6417    1.4015    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.1434    0.9713    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.1833    1.9418    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.0356    0.1309    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.4692    1.9550    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.1833    1.1181    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.3414    0.7897    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.4691    2.5519    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.5066   -2.3725    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.0652    4.2801    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.4695    3.2443    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.3680    4.2801    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.2189    3.9490    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.0755    3.8560    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    6.5773    3.4259    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.6172    4.3963    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.9031    4.4095    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.6172    3.5726    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.7753    3.2443    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.9030    4.9238    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    6.4990    6.5902    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.9033    5.5544    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.8018    6.5902    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.6527    6.2591    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    9.5095    6.1661    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    9.0111    5.7360    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.0510    6.7064    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    8.3369    6.7196    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.0510    5.8827    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.2091    5.5544    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.0687    7.1102    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  2  0      
  4  8  1  0      
  1  7  1  6      
  1  6  1  1      
  9  1  1  0      
 10  9  1  0      
 11 12  1  0      
 11 13  2  0      
 11 10  1  0      
 14 15  1  0      
 15 19  1  1      
 16 20  1  1      
 16 17  1  0      
 14 17  1  0      
 15 18  1  0      
 20 19  1  1      
 17 42  1  0      
 14 21  1  0      
 19 34  1  0      
 20 29  1  0      
 22 23  1  0      
 23 31  1  1      
 24 32  1  1      
 24 25  1  0      
 22 25  1  0      
 24 26  1  0      
 27 32  1  0      
 28 31  1  0      
 23 29  1  0      
 25 30  1  0      
 31 32  1  1      
 30 33  1  0      
 35 36  1  0      
 36 44  1  1      
 37 45  1  1      
 37 38  1  0      
 35 38  1  0      
 37 39  1  0      
 40 45  1  0      
 41 44  1  0      
 36 42  1  0      
 38 43  1  0      
 44 45  1  1      
 43 53  1  0      
 46 47  1  0      
 47 55  1  1      
 48 56  1  1      
 48 49  1  0      
 46 49  1  0      
 48 50  1  0      
 51 56  1  0      
 52 55  1  0      
 47 53  1  0      
 49 54  1  0      
 55 56  1  1      
 54 57  1  0      
 11 58  1  0      
 58 16  1  0      
 59 60  1  0      
 60 67  1  1      
 61 68  1  1      
 61 62  1  0      
 59 62  1  0      
 61 63  1  0      
 64 68  1  0      
 65 67  1  0      
 62 66  1  0      
 67 68  1  1      
 66 69  1  0      
 57 60  1  0      
 70 71  1  0      
 71 78  1  1      
 72 79  1  1      
 72 73  1  0      
 70 73  1  0      
 72 74  1  0      
 75 79  1  0      
 76 78  1  0      
 73 77  1  0      
 78 79  1  1      
 77 80  1  0      
 71 69  1  0      
M  END

ENDPIM5LYSOTEMPLATE

  my($PIM6LYSOSn1TemplateString)=<<ENDPIM6LYSOTEMPLATE;
LYSOPIM6.mol
  ChemDraw01111121522D

 91 97  0  0  0  0  0  0  0  0999 V2000
   -7.8764   -3.4009    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -8.5887   -2.9908    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -9.3013   -3.4009    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  -10.0134   -2.9908    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  -10.0134   -2.1673    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -7.4645   -4.1130    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -8.2881   -4.1130    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  -10.7257   -3.4009    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.1637   -2.9895    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.4513   -3.4009    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.4375   -3.1126    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -5.8004   -3.7416    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.4375   -2.3623    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3787   -3.6788    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7831   -4.7146    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.6814   -3.6788    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5323   -4.0100    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0416   -4.4262    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9308   -4.3863    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0887   -4.7146    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4989   -3.9145    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4931   -7.4661    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0887   -6.4302    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1903   -7.4661    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3394   -7.1349    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5173   -7.0419    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.0190   -6.6118    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9410   -7.5823    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0887   -5.6064    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6553   -7.5954    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9410   -6.7585    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7831   -6.4302    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9775   -8.2624    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5009   -5.1307    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0192   -1.7274    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6148   -2.7630    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.2836   -1.7274    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8656   -2.0584    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.9912   -2.1513    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.4929   -2.5813    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4671   -1.6112    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0917   -3.4364    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1813   -1.5980    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4671   -2.4347    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3092   -2.7630    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4144    0.6729    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1813   -0.3629    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.7170    0.6729    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5679    0.3418    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.4246    0.2488    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.9263   -0.1813    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.9664    0.7891    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1813   -1.0217    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.2521    0.8023    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.9664   -0.0346    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.1243   -0.3629    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.2521    1.3991    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.7232   -3.5249    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.8481    3.1273    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.2524    2.0915    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.1507    3.1273    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.0016    2.7962    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.8583    2.7032    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.3600    2.2731    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.4001    3.2435    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.6858    3.2567    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.4001    2.4198    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.5580    2.0915    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.6857    3.7709    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.2818    5.4372    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.6861    4.4014    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.5844    5.4372    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.4353    5.1061    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.2920    5.0131    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.7937    4.5830    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.8338    5.5534    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.1195    5.5666    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.8338    4.7297    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.9917    4.4014    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.1195    6.0603    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.7155    7.7424    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.1198    6.7066    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.0182    7.7424    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.8691    7.4113    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.7257    7.3183    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   10.2275    6.8882    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    8.2675    7.8586    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    9.5533    7.8718    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.2675    7.0349    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    9.4255    6.7066    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    9.2851    8.2624    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  2  0      
  4  8  1  0      
  1  7  1  6      
  1  6  1  1      
  9  1  1  0      
 10  9  1  0      
 11 12  1  0      
 11 13  2  0      
 11 10  1  0      
 14 15  1  0      
 15 19  1  1      
 16 20  1  1      
 16 17  1  0      
 14 17  1  0      
 15 18  1  0      
 20 19  1  1      
 17 42  1  0      
 14 21  1  0      
 19 34  1  0      
 20 29  1  0      
 22 23  1  0      
 23 31  1  1      
 24 32  1  1      
 24 25  1  0      
 22 25  1  0      
 24 26  1  0      
 27 32  1  0      
 28 31  1  0      
 23 29  1  0      
 25 30  1  0      
 31 32  1  1      
 30 33  1  0      
 35 36  1  0      
 36 44  1  1      
 37 45  1  1      
 37 38  1  0      
 35 38  1  0      
 37 39  1  0      
 40 45  1  0      
 41 44  1  0      
 36 42  1  0      
 38 43  1  0      
 44 45  1  1      
 43 53  1  0      
 46 47  1  0      
 47 55  1  1      
 48 56  1  1      
 48 49  1  0      
 46 49  1  0      
 48 50  1  0      
 51 56  1  0      
 52 55  1  0      
 47 53  1  0      
 49 54  1  0      
 55 56  1  1      
 54 57  1  0      
 11 58  1  0      
 58 16  1  0      
 59 60  1  0      
 60 67  1  1      
 61 68  1  1      
 61 62  1  0      
 59 62  1  0      
 61 63  1  0      
 64 68  1  0      
 65 67  1  0      
 62 66  1  0      
 67 68  1  1      
 66 69  1  0      
 57 60  1  0      
 70 71  1  0      
 71 78  1  1      
 72 79  1  1      
 72 73  1  0      
 70 73  1  0      
 72 74  1  0      
 75 79  1  0      
 76 78  1  0      
 73 77  1  0      
 78 79  1  1      
 77 80  1  0      
 71 69  1  0      
 81 82  1  0      
 82 89  1  1      
 83 90  1  1      
 83 84  1  0      
 81 84  1  0      
 83 85  1  0      
 86 90  1  0      
 87 89  1  0      
 84 88  1  0      
 89 90  1  1      
 88 91  1  0      
 80 82  1  0      
M  END


ENDPIM6LYSOTEMPLATE

  my($DIMBG15TemplateString)=<<ENDDIMBG15TEMPLATE;
GlycDIMB15Et.mol
  ChemDraw01121109052D

 76 79  0  0  0  0  0  0  0  0999 V2000
   10.4533    1.9870    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.0627    1.3103    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.4533    0.6338    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.0627   -0.0427    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.4533   -0.7191    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.0627   -1.3957    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.4533   -2.0722    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.0627   -2.7489    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.4533   -3.4253    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.0627   -4.1019    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.4533   -4.7783    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.0627   -5.4550    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    9.3250    1.3144    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    9.2381   -0.0527    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   11.1695   -3.4380    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    9.3902   -4.1107    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    8.4569   -0.4434    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.4569   -1.1595    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    8.5324    1.6866    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.5268    2.4140    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   10.0627    2.5947    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.4967    3.2673    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.0357    3.8348    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    9.1840    4.9075    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.3009    5.3633    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.4174    4.9075    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.5344    5.3633    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.6512    4.9075    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.7677    5.3633    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.8848    4.9075    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0015    5.3633    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.1183    4.9075    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2350    5.3633    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3520    4.9075    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5923    5.4550    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5386    4.9075    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5386    3.8160    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5923    3.2686    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3520    3.8160    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.3232    3.4215    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.9725    4.2112    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.7722    4.0084    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.0564    3.5981    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.6622    2.8681    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.8625    3.0708    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.5782    3.4812    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.0767    4.5526    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.5454    4.6205    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.0428    2.4976    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.8587    3.4060    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -6.5419    4.1820    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.3416    3.9792    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.6258    3.5689    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.2316    2.8389    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.4319    3.0416    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.1476    3.4520    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -7.6461    4.5234    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -7.1973    4.5946    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.6123    2.4684    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -7.4282    3.3768    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -9.0712    4.1159    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -9.8709    3.9131    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -9.1551    3.5028    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -8.7609    2.7729    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.9612    2.9756    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -8.6769    3.3859    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -9.6764    4.4040    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -9.6941    2.8782    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  -10.5337    3.4219    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -8.3955    2.3238    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   10.5452    4.5475    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    9.9883    5.2694    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  -11.1695    3.4728    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  -10.2606    2.8800    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.7660    2.5330    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.7091    2.5395    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  1  0      
  5  6  1  0      
  6  7  1  0      
  7  8  1  0      
  8  9  1  0      
  9 10  1  0      
 10 11  1  0      
 11 12  1  0      
  2 13  1  1      
  4 14  1  6      
  9 15  1  1      
 10 16  2  0      
 14 17  1  0      
 17 18  2  0      
 13 19  1  0      
 19 20  2  0      
  1 21  1  0      
 21 22  1  0      
 22 23  1  0      
 24 25  1  0      
 25 26  1  0      
 26 27  1  0      
 27 28  1  0      
 28 29  1  0      
 29 30  1  0      
 30 31  1  0      
 31 32  1  0      
 32 33  1  0      
 33 34  1  0      
 34 35  2  0      
 35 36  1  0      
 36 37  2  0      
 37 38  1  0      
 38 39  2  0      
 39 34  1  0      
 40 45  1  0      
 41 42  1  0      
 43 42  1  1      
 43 44  1  1      
 45 44  1  1      
 45 46  1  0      
 46 41  1  0      
 42 47  1  0      
 41 48  1  0      
 44 49  1  0      
 43 50  1  0      
 51 52  1  0      
 52 53  1  1      
 53 54  1  1      
 55 54  1  1      
 55 56  1  0      
 56 51  1  0      
 52 57  1  0      
 51 58  1  0      
 54 59  1  0      
 53 60  1  0      
 50 55  1  0      
 61 62  1  0      
 62 63  1  1      
 63 64  1  1      
 65 64  1  1      
 65 66  1  0      
 66 61  1  0      
 61 67  1  0      
 60 65  1  0      
 63 68  1  0      
 62 69  1  0      
 64 70  1  0      
 37 40  1  0      
 23 71  1  0      
 71 72  1  0      
 24 72  1  0      
 69 73  1  0      
 68 74  1  0      
 70 75  1  0      
 49 76  1  0      
M  END

ENDDIMBG15TEMPLATE

  my($DIMBG15MeTemplateString)=<<ENDDIMBG15METEMPLATE;
GlycDIMB16Et.mol
  ChemDraw01121109192D

 75 78  0  0  0  0  0  0  0  0999 V2000
   10.4534    1.6487    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.0628    0.9720    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.4534    0.2955    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.0628   -0.3811    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.4534   -1.0575    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.0628   -1.7341    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.4534   -2.4106    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.0628   -3.0873    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.4534   -3.7637    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.0628   -4.4403    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.4534   -5.1167    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    9.3251    0.9761    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    9.2382   -0.3911    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   11.1696   -3.7764    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    9.3903   -4.4491    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    8.4570   -0.7818    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.4570   -1.4979    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    8.5325    1.3483    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.5269    2.0757    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   10.0628    2.2564    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.4968    2.9290    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.0358    3.4965    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    9.1841    4.5692    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.3010    5.0250    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.4175    4.5692    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.5345    5.0250    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.6513    4.5692    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.7678    5.0250    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.8849    4.5692    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0015    5.0250    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.1183    4.5692    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2350    5.0250    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3520    4.5692    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5923    5.1167    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5386    4.5692    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5386    3.4777    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5923    2.9303    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3520    3.4777    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.3232    3.0832    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.9726    3.8729    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.7723    3.6701    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.0565    3.2598    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.6622    2.5298    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.8625    2.7325    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.5782    3.1429    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.0768    4.2143    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.5455    4.2822    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.0429    2.1593    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.8588    3.0677    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -6.5420    3.8437    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.3417    3.6409    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.6259    3.2306    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.2317    2.5006    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.4320    2.7033    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.1477    3.1137    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -7.6462    4.1851    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -7.1974    4.2563    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.6124    2.1301    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -7.4283    3.0385    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -9.0713    3.7776    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -9.8710    3.5748    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -9.1552    3.1645    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -8.7610    2.4346    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.9613    2.6373    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -8.6770    3.0476    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -9.6765    4.0657    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -9.6942    2.5399    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  -10.5338    3.0836    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -8.3956    1.9855    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   10.5453    4.2092    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    9.9884    4.9311    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  -11.1696    3.1345    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  -10.2607    2.5417    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.7661    2.1947    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.7092    2.2012    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  1  0      
  5  6  1  0      
  6  7  1  0      
  7  8  1  0      
  8  9  1  0      
  9 10  1  0      
 10 11  1  0      
  2 12  1  1      
  4 13  1  6      
  9 14  1  1      
 10 15  2  0      
 13 16  1  0      
 16 17  2  0      
 12 18  1  0      
 18 19  2  0      
  1 20  1  0      
 20 21  1  0      
 21 22  1  0      
 23 24  1  0      
 24 25  1  0      
 25 26  1  0      
 26 27  1  0      
 27 28  1  0      
 28 29  1  0      
 29 30  1  0      
 30 31  1  0      
 31 32  1  0      
 32 33  1  0      
 33 34  2  0      
 34 35  1  0      
 35 36  2  0      
 36 37  1  0      
 37 38  2  0      
 38 33  1  0      
 39 44  1  0      
 40 41  1  0      
 42 41  1  1      
 42 43  1  1      
 44 43  1  1      
 44 45  1  0      
 45 40  1  0      
 41 46  1  0      
 40 47  1  0      
 43 48  1  0      
 42 49  1  0      
 50 51  1  0      
 51 52  1  1      
 52 53  1  1      
 54 53  1  1      
 54 55  1  0      
 55 50  1  0      
 51 56  1  0      
 50 57  1  0      
 53 58  1  0      
 52 59  1  0      
 49 54  1  0      
 60 61  1  0      
 61 62  1  1      
 62 63  1  1      
 64 63  1  1      
 64 65  1  0      
 65 60  1  0      
 60 66  1  0      
 59 64  1  0      
 62 67  1  0      
 61 68  1  0      
 63 69  1  0      
 36 39  1  0      
 22 70  1  0      
 70 71  1  0      
 23 71  1  0      
 68 72  1  0      
 67 73  1  0      
 69 74  1  0      
 48 75  1  0      
M  END

ENDDIMBG15METEMPLATE

  my($DIMBG17TemplateString)=<<ENDDIMBG17TEMPLATE;
GlycDIMB16Et.mol
  ChemDraw01121109182D

 78 81  0  0  0  0  0  0  0  0999 V2000
   11.1958    1.9664    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.8052    1.2896    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   11.1958    0.6132    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.8052   -0.0633    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   11.1958   -0.7398    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.8052   -1.4163    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   11.1958   -2.0928    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.8052   -2.7695    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   11.1958   -3.4459    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.8052   -4.1225    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   11.1958   -4.7989    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.8052   -5.4756    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.0675    1.2937    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    9.9806   -0.0733    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   11.9120   -3.4586    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.1327   -4.1313    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    9.1994   -0.4640    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    9.1994   -1.1801    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    9.2749    1.6660    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    9.2693    2.3934    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   10.8052    2.5740    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   11.2392    3.2467    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.7782    3.8142    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    9.9265    4.8869    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    9.0434    5.3426    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.1599    4.8869    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.2769    5.3426    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.3937    4.8869    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.5102    5.3426    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.6273    4.8869    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.7440    5.3426    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.8608    4.8869    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.9775    5.3426    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3905    4.9282    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3348    5.4756    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2811    4.9282    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2811    3.8367    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3348    3.2892    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3905    3.8367    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0657    3.4421    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.7150    4.2318    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.5147    4.0290    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.7989    3.6187    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.4047    2.8887    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.6050    3.0914    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.3207    3.5018    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.8192    4.5732    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.2879    4.6411    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.7853    2.5182    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.6012    3.4266    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -7.2844    4.2026    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -8.0841    3.9998    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.3683    3.5895    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.9741    2.8596    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.1744    3.0623    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.8901    3.4726    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -8.3886    4.5441    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -7.9398    4.6152    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.3548    2.4890    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -8.1707    3.3974    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -9.8137    4.1365    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  -10.6134    3.9337    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -9.8976    3.5235    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -9.5034    2.7935    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -8.7037    2.9962    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -9.4194    3.4066    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  -10.4189    4.4246    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  -10.4366    2.8989    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  -11.2762    3.4425    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -9.1380    2.3444    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   11.2877    4.5269    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.7308    5.2487    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  -11.9120    3.4934    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  -11.0031    2.9006    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -8.5085    2.5536    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.4516    2.5601    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2227    4.9187    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4389    5.2900    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  1  0      
  5  6  1  0      
  6  7  1  0      
  7  8  1  0      
  8  9  1  0      
  9 10  1  0      
 10 11  1  0      
 11 12  1  0      
  2 13  1  1      
  4 14  1  6      
  9 15  1  1      
 10 16  2  0      
 14 17  1  0      
 17 18  2  0      
 13 19  1  0      
 19 20  2  0      
  1 21  1  0      
 21 22  1  0      
 22 23  1  0      
 24 25  1  0      
 25 26  1  0      
 26 27  1  0      
 27 28  1  0      
 28 29  1  0      
 29 30  1  0      
 30 31  1  0      
 31 32  1  0      
 32 33  1  0      
 34 35  2  0      
 35 36  1  0      
 36 37  2  0      
 37 38  1  0      
 38 39  2  0      
 39 34  1  0      
 40 45  1  0      
 41 42  1  0      
 43 42  1  1      
 43 44  1  1      
 45 44  1  1      
 45 46  1  0      
 46 41  1  0      
 42 47  1  0      
 41 48  1  0      
 44 49  1  0      
 43 50  1  0      
 51 52  1  0      
 52 53  1  1      
 53 54  1  1      
 55 54  1  1      
 55 56  1  0      
 56 51  1  0      
 52 57  1  0      
 51 58  1  0      
 54 59  1  0      
 53 60  1  0      
 50 55  1  0      
 61 62  1  0      
 62 63  1  1      
 63 64  1  1      
 65 64  1  1      
 65 66  1  0      
 66 61  1  0      
 61 67  1  0      
 60 65  1  0      
 63 68  1  0      
 62 69  1  0      
 64 70  1  0      
 37 40  1  0      
 23 71  1  0      
 71 72  1  0      
 24 72  1  0      
 69 73  1  0      
 68 74  1  0      
 70 75  1  0      
 49 76  1  0      
 33 77  1  0      
 77 78  1  0      
 78 34  1  0      
M  END

ENDDIMBG17TEMPLATE

  my($DIMBG17MeTemplateString)=<<ENDDIMBG17METEMPLATE;
GlycDIMB18Et.mol
  ChemDraw01121109192D

 77 80  0  0  0  0  0  0  0  0999 V2000
   11.1958    1.6280    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.8052    0.9513    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   11.1958    0.2748    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.8052   -0.4017    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   11.1958   -1.0781    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.8052   -1.7547    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   11.1958   -2.4312    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.8052   -3.1078    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   11.1958   -3.7843    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.8052   -4.4609    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   11.1958   -5.1373    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.0675    0.9554    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    9.9806   -0.4117    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   11.9120   -3.7970    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.1327   -4.4697    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    9.1994   -0.8023    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    9.1994   -1.5185    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    9.2749    1.3276    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    9.2693    2.0550    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   10.8052    2.2357    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   11.2392    2.9083    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.7782    3.4758    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    9.9265    4.5485    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    9.0434    5.0043    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.1599    4.5485    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.2769    5.0043    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.3937    4.5485    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.5102    5.0043    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.6273    4.5485    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.7440    5.0043    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.8608    4.5485    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.9775    5.0043    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3905    4.5898    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3348    5.1373    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2811    4.5898    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2811    3.4983    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3348    2.9508    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3905    3.4983    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0657    3.1037    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.7150    3.8934    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.5147    3.6906    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.7989    3.2804    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.4047    2.5504    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.6050    2.7531    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.3207    3.1634    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.8192    4.2349    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.2879    4.3027    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.7853    2.1799    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.6012    3.0882    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -7.2844    3.8642    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -8.0841    3.6615    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.3683    3.2512    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.9741    2.5212    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.1744    2.7239    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.8901    3.1343    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -8.3886    4.2057    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -7.9398    4.2769    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.3548    2.1507    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -8.1707    3.0590    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -9.8137    3.7982    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  -10.6134    3.5954    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -9.8976    3.1851    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -9.5034    2.4551    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -8.7037    2.6578    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -9.4194    3.0682    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  -10.4189    4.0862    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  -10.4366    2.5605    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  -11.2762    3.1041    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -9.1380    2.0061    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   11.2877    4.1885    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.7308    4.9104    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  -11.9120    3.1550    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  -11.0031    2.5623    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -8.5085    2.2152    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.4516    2.2218    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2227    4.5804    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4389    4.9516    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  1  0      
  5  6  1  0      
  6  7  1  0      
  7  8  1  0      
  8  9  1  0      
  9 10  1  0      
 10 11  1  0      
  2 12  1  1      
  4 13  1  6      
  9 14  1  1      
 10 15  2  0      
 13 16  1  0      
 16 17  2  0      
 12 18  1  0      
 18 19  2  0      
  1 20  1  0      
 20 21  1  0      
 21 22  1  0      
 23 24  1  0      
 24 25  1  0      
 25 26  1  0      
 26 27  1  0      
 27 28  1  0      
 28 29  1  0      
 29 30  1  0      
 30 31  1  0      
 31 32  1  0      
 33 34  2  0      
 34 35  1  0      
 35 36  2  0      
 36 37  1  0      
 37 38  2  0      
 38 33  1  0      
 39 44  1  0      
 40 41  1  0      
 42 41  1  1      
 42 43  1  1      
 44 43  1  1      
 44 45  1  0      
 45 40  1  0      
 41 46  1  0      
 40 47  1  0      
 43 48  1  0      
 42 49  1  0      
 50 51  1  0      
 51 52  1  1      
 52 53  1  1      
 54 53  1  1      
 54 55  1  0      
 55 50  1  0      
 51 56  1  0      
 50 57  1  0      
 53 58  1  0      
 52 59  1  0      
 49 54  1  0      
 60 61  1  0      
 61 62  1  1      
 62 63  1  1      
 64 63  1  1      
 64 65  1  0      
 65 60  1  0      
 60 66  1  0      
 59 64  1  0      
 62 67  1  0      
 61 68  1  0      
 63 69  1  0      
 36 39  1  0      
 22 70  1  0      
 70 71  1  0      
 23 71  1  0      
 68 72  1  0      
 67 73  1  0      
 69 74  1  0      
 48 75  1  0      
 32 76  1  0      
 76 77  1  0      
 77 33  1  0      
M  END

ENDDIMBG17METEMPLATE

  my($AC2SGLSn1Sn2TemplateString)=<<ENDAC2SGLSN1SN2TEMPLATE;
Ac2SGL.mol
  ChemDraw01121111282D

 32 33  0  0  0  0  0  0  0  0999 V2000
   -0.6654    2.5899    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0697    1.5540    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.9681    2.5899    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8190    2.2587    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.6757    2.1657    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1775    1.7357    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4006    1.1461    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5033    2.7193    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2174    1.8824    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.3755    1.5540    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5033    3.3070    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.5259   -1.3563    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0697   -0.3204    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.8286   -1.3563    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.6796   -1.0250    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.5362   -0.9321    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.0380   -0.5019    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0697    0.7302    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.0780   -0.6488    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.2360   -0.3204    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2914   -1.2531    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.6757   -0.1338    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0669   -0.8123    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0669   -0.1530    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7699   -1.2181    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.2192   -2.1067    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.5692   -2.6041    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5692   -3.3070    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8131    0.8367    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0
   -2.3287    1.2286    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2356    0.3417    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2669    0.4654    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  9  1  1      
  3 10  1  1      
  3  4  1  0      
  1  4  1  0      
  3  5  1  0      
  6 10  1  0      
  2 18  1  0      
  7  9  1  0      
  4  8  1  0      
  9 10  1  1      
  8 11  1  0      
 13 12  1  1      
 13 19  1  0      
 14 20  1  0      
 14 15  1  1      
 12 15  1  1      
 14 16  1  0      
 17 20  1  0      
 13 18  1  0      
 19 20  1  0      
 12 21  1  0      
 15 26  1  0      
 17 22  1  0      
 21 23  1  0      
 23 24  2  0      
 23 25  1  0      
 26 27  1  0      
 27 28  2  0      
  7 29  1  0      
 29 30  2  0      
 29 31  2  0      
 29 32  1  0      
M  END

ENDAC2SGLSN1SN2TEMPLATE

  my($DIMAG15TemplateString)=<<ENDDIMAG15TEMPLATE;
GlycDIMA15Et.mol
  ChemDraw01121112232D

 77 80  0  0  0  0  0  0  0  0999 V2000
   10.4534    1.9870    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.0628    1.3103    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.4534    0.6338    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.0628   -0.0427    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.4534   -0.7191    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.0628   -1.3957    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.4534   -2.0722    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.0628   -2.7489    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.4534   -3.4253    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.0628   -4.1020    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.4534   -4.7784    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.0628   -5.4551    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    9.3251    1.3144    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    9.2382   -0.0527    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   11.1696   -3.4380    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    9.3903   -4.1108    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    8.4570   -0.4434    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.4570   -1.1595    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    8.5325    1.6866    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.5269    2.4140    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   10.0628    2.5947    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.4968    3.2673    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.0358    3.8349    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    9.1841    4.9076    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.3010    5.3634    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.4175    4.9076    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.5345    5.3634    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.6513    4.9076    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.7678    5.3634    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.8849    4.9076    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0015    5.3634    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.1183    4.9076    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2350    5.3634    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3520    4.9076    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5923    5.4551    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5386    4.9076    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5386    3.8161    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5923    3.2686    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3520    3.8161    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.3232    3.4215    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.9726    4.2113    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.7723    4.0085    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.0565    3.5981    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.6622    2.8681    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.8625    3.0708    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.5782    3.4812    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.0768    4.5527    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.5455    4.6206    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.0429    2.4976    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.8588    3.4060    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -6.5420    4.1821    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.3417    3.9793    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.6259    3.5689    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.2317    2.8389    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.4320    3.0416    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.1477    3.4520    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -7.6462    4.5235    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -7.1974    4.5947    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.6124    2.4684    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -7.4283    3.3768    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -9.0713    4.1160    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -9.8710    3.9132    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -9.1552    3.5028    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -8.7610    2.7729    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.9613    2.9756    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -8.6770    3.3859    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -9.6765    4.4041    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -9.6942    2.8782    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  -10.5338    3.4219    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -8.3956    2.3238    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   10.5453    4.5476    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    9.9884    5.2695    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  -11.1696    3.4728    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  -10.2607    2.8800    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.7661    2.5330    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.7092    2.5395    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.9383   -4.7392    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  1  0      
  5  6  1  0      
  6  7  1  0      
  7  8  1  0      
  8  9  1  0      
  9 10  1  0      
 10 11  1  0      
 11 12  1  0      
  2 13  1  1      
  4 14  1  6      
  9 15  1  1      
 10 16  1  0      
 14 17  1  0      
 17 18  2  0      
 13 19  1  0      
 19 20  2  0      
  1 21  1  0      
 21 22  1  0      
 22 23  1  0      
 24 25  1  0      
 25 26  1  0      
 26 27  1  0      
 27 28  1  0      
 28 29  1  0      
 29 30  1  0      
 30 31  1  0      
 31 32  1  0      
 32 33  1  0      
 33 34  1  0      
 34 35  2  0      
 35 36  1  0      
 36 37  2  0      
 37 38  1  0      
 38 39  2  0      
 39 34  1  0      
 40 45  1  0      
 41 42  1  0      
 43 42  1  1      
 43 44  1  1      
 45 44  1  1      
 45 46  1  0      
 46 41  1  0      
 42 47  1  0      
 41 48  1  0      
 44 49  1  0      
 43 50  1  0      
 51 52  1  0      
 52 53  1  1      
 53 54  1  1      
 55 54  1  1      
 55 56  1  0      
 56 51  1  0      
 52 57  1  0      
 51 58  1  0      
 54 59  1  0      
 53 60  1  0      
 50 55  1  0      
 61 62  1  0      
 62 63  1  1      
 63 64  1  1      
 65 64  1  1      
 65 66  1  0      
 66 61  1  0      
 61 67  1  0      
 60 65  1  0      
 63 68  1  0      
 62 69  1  0      
 64 70  1  0      
 37 40  1  0      
 23 71  1  0      
 71 72  1  0      
 24 72  1  0      
 69 73  1  0      
 68 74  1  0      
 70 75  1  0      
 49 76  1  0      
 16 77  1  0      
M  END

ENDDIMAG15TEMPLATE

  my($DIMAG15MeTemplateString)=<<ENDDIMAG15METEMPLATE;
GlycDIMA15Et.mol
  ChemDraw01121112242D

 76 79  0  0  0  0  0  0  0  0999 V2000
   10.4534    1.6487    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.0628    0.9720    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.4534    0.2955    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.0628   -0.3811    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.4534   -1.0575    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.0628   -1.7341    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.4534   -2.4106    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.0628   -3.0873    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.4534   -3.7637    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.0628   -4.4403    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.4534   -5.1167    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    9.3251    0.9761    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    9.2382   -0.3911    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   11.1696   -3.7764    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    9.3903   -4.4491    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    8.4570   -0.7818    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.4570   -1.4979    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    8.5325    1.3483    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.5269    2.0757    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   10.0628    2.2564    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.4968    2.9290    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.0358    3.4965    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    9.1841    4.5692    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.3010    5.0250    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.4175    4.5692    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.5345    5.0250    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.6513    4.5692    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.7678    5.0250    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.8849    4.5692    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0015    5.0250    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.1183    4.5692    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2350    5.0250    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3520    4.5692    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5923    5.1167    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5386    4.5692    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5386    3.4777    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5923    2.9303    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3520    3.4777    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.3232    3.0832    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.9726    3.8729    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.7723    3.6701    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.0565    3.2598    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.6622    2.5298    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.8625    2.7325    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.5782    3.1429    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.0768    4.2143    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.5455    4.2822    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.0429    2.1593    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.8588    3.0677    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -6.5420    3.8437    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.3417    3.6409    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.6259    3.2306    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.2317    2.5006    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.4320    2.7033    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.1477    3.1137    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -7.6462    4.1851    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -7.1974    4.2563    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.6124    2.1301    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -7.4283    3.0385    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -9.0713    3.7776    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -9.8710    3.5748    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -9.1552    3.1645    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -8.7610    2.4346    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.9613    2.6373    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -8.6770    3.0476    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -9.6765    4.0657    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -9.6942    2.5399    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  -10.5338    3.0836    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -8.3956    1.9855    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   10.5453    4.2092    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    9.9884    4.9311    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  -11.1696    3.1345    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  -10.2607    2.5417    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.7661    2.1947    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.7092    2.2012    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.9383   -5.0776    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  1  0      
  5  6  1  0      
  6  7  1  0      
  7  8  1  0      
  8  9  1  0      
  9 10  1  0      
 10 11  1  0      
  2 12  1  1      
  4 13  1  6      
  9 14  1  1      
 10 15  1  0      
 13 16  1  0      
 16 17  2  0      
 12 18  1  0      
 18 19  2  0      
  1 20  1  0      
 20 21  1  0      
 21 22  1  0      
 23 24  1  0      
 24 25  1  0      
 25 26  1  0      
 26 27  1  0      
 27 28  1  0      
 28 29  1  0      
 29 30  1  0      
 30 31  1  0      
 31 32  1  0      
 32 33  1  0      
 33 34  2  0      
 34 35  1  0      
 35 36  2  0      
 36 37  1  0      
 37 38  2  0      
 38 33  1  0      
 39 44  1  0      
 40 41  1  0      
 42 41  1  1      
 42 43  1  1      
 44 43  1  1      
 44 45  1  0      
 45 40  1  0      
 41 46  1  0      
 40 47  1  0      
 43 48  1  0      
 42 49  1  0      
 50 51  1  0      
 51 52  1  1      
 52 53  1  1      
 54 53  1  1      
 54 55  1  0      
 55 50  1  0      
 51 56  1  0      
 50 57  1  0      
 53 58  1  0      
 52 59  1  0      
 49 54  1  0      
 60 61  1  0      
 61 62  1  1      
 62 63  1  1      
 64 63  1  1      
 64 65  1  0      
 65 60  1  0      
 60 66  1  0      
 59 64  1  0      
 62 67  1  0      
 61 68  1  0      
 63 69  1  0      
 36 39  1  0      
 22 70  1  0      
 70 71  1  0      
 23 71  1  0      
 68 72  1  0      
 67 73  1  0      
 69 74  1  0      
 48 75  1  0      
 15 76  1  0      
M  END

ENDDIMAG15METEMPLATE

  my($DIMAG17TemplateString)=<<ENDDIMAG17TEMPLATE;
GlycDIMB17Et.mol
  ChemDraw01121112252D

 79 82  0  0  0  0  0  0  0  0999 V2000
   11.1959    1.9664    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.8053    1.2896    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   11.1959    0.6132    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.8053   -0.0633    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   11.1959   -0.7398    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.8053   -1.4163    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   11.1959   -2.0928    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.8053   -2.7695    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   11.1959   -3.4459    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.8053   -4.1226    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   11.1959   -4.7990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.8053   -5.4757    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.0676    1.2937    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    9.9807   -0.0733    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   11.9121   -3.4586    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.1328   -4.1314    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    9.1995   -0.4640    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    9.1995   -1.1801    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    9.2750    1.6660    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    9.2694    2.3934    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   10.8053    2.5740    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   11.2393    3.2467    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.7783    3.8142    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    9.9266    4.8870    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    9.0435    5.3427    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.1600    4.8870    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.2770    5.3427    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.3938    4.8870    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.5103    5.3427    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.6274    4.8870    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.7440    5.3427    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.8608    4.8870    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.9775    5.3427    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3905    4.9283    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3348    5.4757    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2811    4.9283    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2811    3.8367    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3348    3.2892    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3905    3.8367    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0657    3.4421    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.7151    4.2319    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.5148    4.0290    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.7990    3.6187    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.4048    2.8887    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.6050    3.0914    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.3208    3.5018    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.8193    4.5733    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.2880    4.6412    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.7854    2.5182    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.6013    3.4266    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -7.2845    4.2027    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -8.0842    3.9998    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.3684    3.5895    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.9742    2.8596    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.1745    3.0623    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.8902    3.4726    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -8.3887    4.5442    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -7.9399    4.6153    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.3549    2.4890    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -8.1708    3.3974    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -9.8138    4.1366    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  -10.6135    3.9337    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -9.8977    3.5235    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -9.5035    2.7935    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -8.7038    2.9962    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -9.4195    3.4066    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  -10.4190    4.4247    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  -10.4367    2.8989    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  -11.2763    3.4425    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -9.1381    2.3444    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   11.2878    4.5270    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.7309    5.2488    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  -11.9121    3.4934    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  -11.0032    2.9006    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -8.5086    2.5536    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.4517    2.5601    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2227    4.9188    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4389    5.2901    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    9.7220   -4.7496    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  1  0      
  5  6  1  0      
  6  7  1  0      
  7  8  1  0      
  8  9  1  0      
  9 10  1  0      
 10 11  1  0      
 11 12  1  0      
  2 13  1  1      
  4 14  1  6      
  9 15  1  1      
 10 16  1  0      
 14 17  1  0      
 17 18  2  0      
 13 19  1  0      
 19 20  2  0      
  1 21  1  0      
 21 22  1  0      
 22 23  1  0      
 24 25  1  0      
 25 26  1  0      
 26 27  1  0      
 27 28  1  0      
 28 29  1  0      
 29 30  1  0      
 30 31  1  0      
 31 32  1  0      
 32 33  1  0      
 34 35  2  0      
 35 36  1  0      
 36 37  2  0      
 37 38  1  0      
 38 39  2  0      
 39 34  1  0      
 40 45  1  0      
 41 42  1  0      
 43 42  1  1      
 43 44  1  1      
 45 44  1  1      
 45 46  1  0      
 46 41  1  0      
 42 47  1  0      
 41 48  1  0      
 44 49  1  0      
 43 50  1  0      
 51 52  1  0      
 52 53  1  1      
 53 54  1  1      
 55 54  1  1      
 55 56  1  0      
 56 51  1  0      
 52 57  1  0      
 51 58  1  0      
 54 59  1  0      
 53 60  1  0      
 50 55  1  0      
 61 62  1  0      
 62 63  1  1      
 63 64  1  1      
 65 64  1  1      
 65 66  1  0      
 66 61  1  0      
 61 67  1  0      
 60 65  1  0      
 63 68  1  0      
 62 69  1  0      
 64 70  1  0      
 37 40  1  0      
 23 71  1  0      
 71 72  1  0      
 24 72  1  0      
 69 73  1  0      
 68 74  1  0      
 70 75  1  0      
 49 76  1  0      
 33 77  1  0      
 77 78  1  0      
 78 34  1  0      
 16 79  1  0      
M  END

ENDDIMAG17TEMPLATE

  my($DIMAG17MeTemplateString)=<<ENDDIMAG17METEMPLATE;
GlycDIMA17Et.mol
  ChemDraw01121112252D

 78 81  0  0  0  0  0  0  0  0999 V2000
   11.1959    1.6281    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.8053    0.9513    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   11.1959    0.2749    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.8053   -0.4017    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   11.1959   -1.0782    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.8053   -1.7547    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   11.1959   -2.4312    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.8053   -3.1079    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   11.1959   -3.7843    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.8053   -4.4609    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   11.1959   -5.1373    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.0676    0.9554    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    9.9807   -0.4117    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   11.9121   -3.7970    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.1328   -4.4697    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    9.1995   -0.8024    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    9.1995   -1.5185    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    9.2750    1.3277    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    9.2694    2.0551    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   10.8053    2.2357    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   11.2393    2.9084    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.7783    3.4759    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    9.9266    4.5486    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    9.0435    5.0043    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.1600    4.5486    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.2770    5.0043    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.3938    4.5486    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.5103    5.0043    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.6274    4.5486    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.7440    5.0043    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.8608    4.5486    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.9775    5.0043    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3905    4.5899    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3348    5.1373    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2811    4.5899    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2811    3.4984    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3348    2.9509    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3905    3.4984    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0657    3.1038    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.7151    3.8935    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.5148    3.6907    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.7990    3.2804    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.4048    2.5504    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.6050    2.7531    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.3208    3.1635    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.8193    4.2349    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.2880    4.3028    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.7854    2.1799    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.6013    3.0883    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -7.2845    3.8643    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -8.0842    3.6615    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.3684    3.2512    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.9742    2.5213    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.1745    2.7240    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.8902    3.1343    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -8.3887    4.2058    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -7.9399    4.2769    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.3549    2.1507    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -8.1708    3.0591    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -9.8138    3.7982    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  -10.6135    3.5954    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -9.8977    3.1852    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -9.5035    2.4552    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -8.7038    2.6579    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -9.4195    3.0683    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  -10.4190    4.0863    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  -10.4367    2.5606    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  -11.2763    3.1042    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -9.1381    2.0061    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   11.2878    4.1886    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.7309    4.9104    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  -11.9121    3.1551    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  -11.0032    2.5623    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -8.5086    2.2153    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.4517    2.2218    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2227    4.5804    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4389    4.9517    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    9.7220   -5.0879    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  1  0      
  5  6  1  0      
  6  7  1  0      
  7  8  1  0      
  8  9  1  0      
  9 10  1  0      
 10 11  1  0      
  2 12  1  1      
  4 13  1  6      
  9 14  1  1      
 10 15  1  0      
 13 16  1  0      
 16 17  2  0      
 12 18  1  0      
 18 19  2  0      
  1 20  1  0      
 20 21  1  0      
 21 22  1  0      
 23 24  1  0      
 24 25  1  0      
 25 26  1  0      
 26 27  1  0      
 27 28  1  0      
 28 29  1  0      
 29 30  1  0      
 30 31  1  0      
 31 32  1  0      
 33 34  2  0      
 34 35  1  0      
 35 36  2  0      
 36 37  1  0      
 37 38  2  0      
 38 33  1  0      
 39 44  1  0      
 40 41  1  0      
 42 41  1  1      
 42 43  1  1      
 44 43  1  1      
 44 45  1  0      
 45 40  1  0      
 41 46  1  0      
 40 47  1  0      
 43 48  1  0      
 42 49  1  0      
 50 51  1  0      
 51 52  1  1      
 52 53  1  1      
 54 53  1  1      
 54 55  1  0      
 55 50  1  0      
 51 56  1  0      
 50 57  1  0      
 53 58  1  0      
 52 59  1  0      
 49 54  1  0      
 60 61  1  0      
 61 62  1  1      
 62 63  1  1      
 64 63  1  1      
 64 65  1  0      
 65 60  1  0      
 60 66  1  0      
 59 64  1  0      
 62 67  1  0      
 61 68  1  0      
 63 69  1  0      
 36 39  1  0      
 22 70  1  0      
 70 71  1  0      
 23 71  1  0      
 68 72  1  0      
 67 73  1  0      
 69 74  1  0      
 48 75  1  0      
 32 76  1  0      
 76 77  1  0      
 77 33  1  0      
 15 78  1  0      
M  END

ENDDIMAG17METEMPLATE

  my($Ac1PIM1Sn1Sn2Sn3TemplateString)=<<ENDAC1PIM1TEMPLATE;
Ac1PIM1.mol
  ChemDraw01191113062D

 42 43  0  0  0  0  0  0  0  0999 V2000
   -2.7725    2.3263    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.4848    2.7364    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.1976    2.3263    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.9097    2.7364    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.9097    3.5600    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.3605    1.6139    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1842    1.6139    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.6221    2.3263    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0597    2.7377    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3473    2.3263    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3334    2.6146    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6964    1.9854    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3334    3.3650    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.9306    1.1928    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.9306    0.3690    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.6431    1.6044    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.7260    2.0482    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.3217    1.0123    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4230    2.0482    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.5722    1.7170    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.0632    1.3007    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.1739    1.3406    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0158    1.0123    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.6059    1.8125    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.6114   -1.7395    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.0158   -0.7035    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.9145   -1.7395    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.7653   -1.4082    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.6221   -1.3152    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.1238   -0.8851    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.1637   -1.8557    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.0158    0.1204    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.4495   -1.8688    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.1637   -1.0318    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.3217   -0.7035    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.1273   -2.5359    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.6038    0.5961    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.0129    2.2907    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.3811    2.2021    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.6325   -2.8588    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.6325   -3.5600    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.9932   -2.4463    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  2  0      
  4  8  1  0      
  1  7  1  6      
  1  6  1  1      
  9  1  1  0      
 10  9  1  0      
 11 12  1  0      
 11 13  2  0      
 14 15  2  0      
 14 16  1  0      
 14  7  1  0      
 11 10  1  0      
 17 18  1  0      
 18 22  1  1      
 19 23  1  1      
 19 20  1  0      
 17 20  1  0      
 18 21  1  0      
 23 22  1  1      
 20 38  1  0      
 17 24  1  0      
 22 37  1  0      
 23 32  1  0      
 25 26  1  0      
 26 34  1  1      
 27 35  1  1      
 27 28  1  0      
 25 28  1  0      
 27 29  1  0      
 30 35  1  0      
 31 34  1  0      
 26 32  1  0      
 28 33  1  0      
 34 35  1  1      
 33 36  1  0      
 11 39  1  0      
 39 19  1  0      
 36 40  1  0      
 40 41  2  0      
 40 42  1  0      
M  END

ENDAC1PIM1TEMPLATE

  my($Ac1PIM2Sn1Sn2Sn3TemplateString)=<<ENDAC1PIM2TEMPLATE;
PIM2.mol
  ChemDraw01191113192D

 53 55  0  0  0  0  0  0  0  0999 V2000
   -3.0503    1.6288    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.7626    2.0388    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.4753    1.6288    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.1872    2.0388    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.1872    2.8622    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6384    0.9165    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -3.4620    0.9165    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.8996    1.6288    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.3376    2.0401    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6254    1.6288    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6116    1.9171    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9745    1.2879    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6116    2.6672    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.2083    0.4954    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.2083   -0.3284    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.9207    0.9070    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.4474    1.3507    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.0430    0.3149    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.1447    1.3507    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.2937    1.0195    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.7845    0.6033    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.8953    0.6432    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.7373    0.3149    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.3272    1.1150    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.3329   -2.4365    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.7373   -1.4006    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.6358   -2.4365    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.4867   -2.1053    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.3432   -2.0123    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.8450   -1.5822    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.8851   -2.5527    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.7373   -0.5770    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.1708   -2.5658    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.8851   -1.7289    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.0430   -1.4006    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.8486   -3.2329    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.3252   -0.1013    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.8894    3.6116    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.2937    2.5758    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.1920    3.6116    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.0430    3.2804    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.8996    3.1875    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.4013    2.7575    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.4415    3.7278    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.7343    1.5932    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.7273    3.7410    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.4415    2.9041    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.5994    2.5758    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.4367    4.3011    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.1029    1.5046    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.1994   -3.6617    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.1994   -4.3011    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.5394   -3.2080    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  2  0      
  4  8  1  0      
  1  7  1  6      
  1  6  1  1      
  9  1  1  0      
 10  9  1  0      
 11 12  1  0      
 11 13  2  0      
 14 15  2  0      
 14 16  1  0      
 14  7  1  0      
 11 10  1  0      
 17 18  1  0      
 18 22  1  1      
 19 23  1  1      
 19 20  1  0      
 17 20  1  0      
 18 21  1  0      
 23 22  1  1      
 20 45  1  0      
 17 24  1  0      
 22 37  1  0      
 23 32  1  0      
 25 26  1  0      
 26 34  1  1      
 27 35  1  1      
 27 28  1  0      
 25 28  1  0      
 27 29  1  0      
 30 35  1  0      
 31 34  1  0      
 26 32  1  0      
 28 33  1  0      
 34 35  1  1      
 33 36  1  0      
 38 39  1  0      
 39 47  1  1      
 40 48  1  1      
 40 41  1  0      
 38 41  1  0      
 40 42  1  0      
 43 48  1  0      
 44 47  1  0      
 39 45  1  0      
 41 46  1  0      
 47 48  1  1      
 46 49  1  0      
 11 50  1  0      
 50 19  1  0      
 36 51  1  0      
 51 52  2  0      
 51 53  1  0      
M  END

ENDAC1PIM2TEMPLATE

  my($Ac1PIM3Sn1Sn2Sn3TemplateString)=<<ENDAC1PIM3TEMPLATE;
PIM3.mol
  ChemDraw01191113192D
  
 64 67  0  0  0  0  0  0  0  0999 V2000
   -4.2671    0.3269    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.9794    0.7370    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.6920    0.3269    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -6.4041    0.7370    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.4041    1.5605    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.8552   -0.3854    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -4.6788   -0.3854    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -7.1165    0.3269    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.5544    0.7383    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.8421    0.3269    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8283    0.6152    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1912   -0.0140    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8283    1.3655    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.4250   -0.8065    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.4250   -1.6302    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -6.1375   -0.3949    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.2307    0.0488    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.8263   -0.9870    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0721    0.0488    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0771   -0.2824    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.5677   -0.6986    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.6786   -0.6587    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5207   -0.9870    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.1105   -0.1869    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.1163   -3.7384    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.5207   -2.7025    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.4191   -3.7384    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.2700   -3.4072    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.1266   -3.3142    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.6283   -2.8841    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.6684   -3.8546    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.5207   -1.8788    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.9541   -3.8677    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.6684   -3.0308    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.8263   -2.7025    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.6319   -4.5347    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.1085   -1.4031    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.6727    2.3098    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.0771    1.2741    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.9754    2.3098    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.8263    1.9787    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.6830    1.8858    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.1847    1.4558    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.2248    2.4260    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.5177    0.2913    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.5105    2.4392    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.2248    1.6024    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.3827    1.2741    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.1062    5.1225    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.5105    4.0868    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.4089    5.1225    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.2598    4.7914    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.1165    4.6984    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    6.6182    4.2684    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.6582    5.2387    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.5105    3.2630    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.9440    5.2519    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.6582    4.4151    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.8162    4.0868    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.6758    5.6425    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1139    0.2027    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.0311   -4.9825    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0311   -5.6425    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.3917   -4.5494    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  2  0      
  4  8  1  0      
  1  7  1  6      
  1  6  1  1      
  9  1  1  0      
 10  9  1  0      
 11 12  1  0      
 11 13  2  0      
 14 15  2  0      
 14 16  1  0      
 14  7  1  0      
 11 10  1  0      
 17 18  1  0      
 18 22  1  1      
 19 23  1  1      
 19 20  1  0      
 17 20  1  0      
 18 21  1  0      
 23 22  1  1      
 20 45  1  0      
 17 24  1  0      
 22 37  1  0      
 23 32  1  0      
 25 26  1  0      
 26 34  1  1      
 27 35  1  1      
 27 28  1  0      
 25 28  1  0      
 27 29  1  0      
 30 35  1  0      
 31 34  1  0      
 26 32  1  0      
 28 33  1  0      
 34 35  1  1      
 33 36  1  0      
 38 39  1  0      
 39 47  1  1      
 40 48  1  1      
 40 41  1  0      
 38 41  1  0      
 40 42  1  0      
 43 48  1  0      
 44 47  1  0      
 39 45  1  0      
 41 46  1  0      
 47 48  1  1      
 46 56  1  0      
 49 50  1  0      
 50 58  1  1      
 51 59  1  1      
 51 52  1  0      
 49 52  1  0      
 51 53  1  0      
 54 59  1  0      
 55 58  1  0      
 50 56  1  0      
 52 57  1  0      
 58 59  1  1      
 57 60  1  0      
 11 61  1  0      
 61 19  1  0      
 36 62  1  0      
 62 63  2  0      
 62 64  1  0      
M  END

ENDAC1PIM3TEMPLATE

  my($Ac1PIM4Sn1Sn2Sn3TemplateString)=<<ENDAC1PIM4TEMPLATE;
PIM4.mol
  ChemDraw01191113462D

 75 79  0  0  0  0  0  0  0  0999 V2000
   -5.4427   -0.5599    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.1550   -0.1499    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.8675   -0.5599    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -7.5796   -0.1499    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.5796    0.6736    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.0308   -1.2722    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -5.8544   -1.2722    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -8.2920   -0.5599    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.7300   -0.1486    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.0176   -0.5599    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0039   -0.2717    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -3.3668   -0.9008    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0039    0.4786    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -6.6005   -1.6933    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.6005   -2.5170    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -7.3130   -1.2817    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0550   -0.8380    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.6506   -1.8738    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2477   -0.8380    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0985   -1.1692    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.3921   -1.5854    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.5030   -1.5455    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6549   -1.8738    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.9348   -1.0737    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0593   -4.6252    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6549   -3.5893    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.2434   -4.6252    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0943   -4.2940    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.9510   -4.2010    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.4527   -3.7709    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.4928   -4.7414    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6549   -2.7655    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.7784   -4.7545    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4928   -3.9176    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.6506   -3.5893    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4562   -5.4215    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.9328   -2.2899    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.4146    1.1135    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1810    0.0779    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.7173    1.1135    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5681    0.7825    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.4249    0.6896    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.9266    0.2596    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.9666    1.2297    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.3421   -0.5955    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.2524    1.2429    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.9666    0.4062    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.1245    0.0779    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.8481    3.5138    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.2524    2.4780    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.1507    3.5138    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.0016    3.1827    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.8583    3.0897    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.3600    2.6596    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.4001    3.6300    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.2524    1.8192    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.6858    3.6432    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.4001    2.8063    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.5580    2.4780    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.6857    4.2401    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2895   -0.6841    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.2818    5.9681    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.6861    4.9324    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.5844    5.9681    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.4353    5.6370    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.2920    5.5440    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.7937    5.1139    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.8338    6.0843    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.1195    6.0975    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.8338    5.2606    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.9917    4.9324    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.8513    6.4881    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.8556   -5.8488    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.8556   -6.4881    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.2163   -5.3950    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  2  0      
  4  8  1  0      
  1  7  1  6      
  1  6  1  1      
  9  1  1  0      
 10  9  1  0      
 11 12  1  0      
 11 13  2  0      
 14 15  2  0      
 14 16  1  0      
 14  7  1  0      
 11 10  1  0      
 17 18  1  0      
 18 22  1  1      
 19 23  1  1      
 19 20  1  0      
 17 20  1  0      
 18 21  1  0      
 23 22  1  1      
 20 45  1  0      
 17 24  1  0      
 22 37  1  0      
 23 32  1  0      
 25 26  1  0      
 26 34  1  1      
 27 35  1  1      
 27 28  1  0      
 25 28  1  0      
 27 29  1  0      
 30 35  1  0      
 31 34  1  0      
 26 32  1  0      
 28 33  1  0      
 34 35  1  1      
 33 36  1  0      
 38 39  1  0      
 39 47  1  1      
 40 48  1  1      
 40 41  1  0      
 38 41  1  0      
 40 42  1  0      
 43 48  1  0      
 44 47  1  0      
 39 45  1  0      
 41 46  1  0      
 47 48  1  1      
 46 56  1  0      
 49 50  1  0      
 50 58  1  1      
 51 59  1  1      
 51 52  1  0      
 49 52  1  0      
 51 53  1  0      
 54 59  1  0      
 55 58  1  0      
 50 56  1  0      
 52 57  1  0      
 58 59  1  1      
 57 60  1  0      
 11 61  1  0      
 61 19  1  0      
 62 63  1  0      
 63 70  1  1      
 64 71  1  1      
 64 65  1  0      
 62 65  1  0      
 64 66  1  0      
 67 71  1  0      
 68 70  1  0      
 65 69  1  0      
 70 71  1  1      
 69 72  1  0      
 60 63  1  0      
 36 73  1  0      
 73 74  2  0      
 73 75  1  0      
M  END

ENDAC1PIM4TEMPLATE

  my($Ac1PIM5Sn1Sn2Sn3TemplateString)=<<ENDAC1PIM5TEMPLATE;
PIM5.mol
  ChemDraw01191113542D

 86 91  0  0  0  0  0  0  0  0999 V2000
   -6.6599   -1.6806    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.3723   -1.2705    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -8.0849   -1.6806    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -8.7970   -1.2705    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -8.7970   -0.4470    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -6.2480   -2.3929    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -7.0716   -2.3929    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -9.5095   -1.6806    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.9472   -1.2692    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.2349   -1.6806    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.2210   -1.3923    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -4.5839   -2.0215    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.2210   -0.6420    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -7.8179   -2.8140    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.8179   -3.6378    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -8.5304   -2.4024    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1618   -1.9587    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4338   -2.9945    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4647   -1.9587    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3155   -2.2899    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.1753   -2.7061    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7139   -2.6662    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8719   -2.9945    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7180   -2.1944    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2763   -5.7462    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8719   -4.7102    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0266   -5.7462    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1225   -5.4149    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.7342   -5.3219    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.2359   -4.8918    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7241   -5.8624    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8719   -3.8864    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.5616   -5.8755    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7241   -5.0385    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4338   -4.7102    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.2394   -6.5425    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2840   -3.4106    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8023   -0.0071    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3980   -1.0427    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5005   -0.0071    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3513   -0.3381    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.2081   -0.4310    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.7098   -0.8610    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2502    0.1091    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8748   -1.7162    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.0356    0.1223    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2502   -0.7144    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.9077   -1.0427    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.6313    2.3933    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.0356    1.3575    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.9341    2.3933    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.7849    2.0622    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.6417    1.9692    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.1434    1.5391    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.1833    2.5095    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.0356    0.6987    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.4692    2.5227    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.1833    1.6858    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.3414    1.3575    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.4691    3.1197    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.5066   -1.8048    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.0652    4.8479    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.4695    3.8120    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.3680    4.8479    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.2189    4.5167    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.0755    4.4237    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    6.5773    3.9936    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.6172    4.9641    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.9031    4.9773    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.6172    4.1403    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.7753    3.8120    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.9030    5.4916    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    6.4990    7.1580    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.9033    6.1221    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.8018    7.1580    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.6527    6.8269    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    9.5095    6.7339    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    9.0111    6.3037    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.0510    7.2742    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    8.3369    7.2874    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.0510    6.4504    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.2091    6.1221    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.0687    7.6780    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4856   -7.0180    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4856   -7.6780    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1456   -6.5849    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  2  0      
  4  8  1  0      
  1  7  1  6      
  1  6  1  1      
  9  1  1  0      
 10  9  1  0      
 11 12  1  0      
 11 13  2  0      
 14 15  2  0      
 14 16  1  0      
 14  7  1  0      
 11 10  1  0      
 17 18  1  0      
 18 22  1  1      
 19 23  1  1      
 19 20  1  0      
 17 20  1  0      
 18 21  1  0      
 23 22  1  1      
 20 45  1  0      
 17 24  1  0      
 22 37  1  0      
 23 32  1  0      
 25 26  1  0      
 26 34  1  1      
 27 35  1  1      
 27 28  1  0      
 25 28  1  0      
 27 29  1  0      
 30 35  1  0      
 31 34  1  0      
 26 32  1  0      
 28 33  1  0      
 34 35  1  1      
 33 36  1  0      
 38 39  1  0      
 39 47  1  1      
 40 48  1  1      
 40 41  1  0      
 38 41  1  0      
 40 42  1  0      
 43 48  1  0      
 44 47  1  0      
 39 45  1  0      
 41 46  1  0      
 47 48  1  1      
 46 56  1  0      
 49 50  1  0      
 50 58  1  1      
 51 59  1  1      
 51 52  1  0      
 49 52  1  0      
 51 53  1  0      
 54 59  1  0      
 55 58  1  0      
 50 56  1  0      
 52 57  1  0      
 58 59  1  1      
 57 60  1  0      
 11 61  1  0      
 61 19  1  0      
 62 63  1  0      
 63 70  1  1      
 64 71  1  1      
 64 65  1  0      
 62 65  1  0      
 64 66  1  0      
 67 71  1  0      
 68 70  1  0      
 65 69  1  0      
 70 71  1  1      
 69 72  1  0      
 60 63  1  0      
 73 74  1  0      
 74 81  1  1      
 75 82  1  1      
 75 76  1  0      
 73 76  1  0      
 75 77  1  0      
 78 82  1  0      
 79 81  1  0      
 76 80  1  0      
 81 82  1  1      
 80 83  1  0      
 74 72  1  0      
 36 84  1  0      
 84 85  2  0      
 84 86  1  0      
M  END

ENDAC1PIM5TEMPLATE

  my($Ac1PIM6Sn1Sn2Sn3TemplateString)=<<ENDAC1PIM6TEMPLATE;
PIM6.mol
  ChemDraw01191113582D

 97103  0  0  0  0  0  0  0  0999 V2000
   -7.8764   -2.7890    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -8.5887   -2.3789    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -9.3013   -2.7890    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  -10.0134   -2.3789    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  -10.0134   -1.5554    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -7.4645   -3.5012    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -8.2881   -3.5012    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  -10.7257   -2.7890    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.1637   -2.3776    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.4513   -2.7890    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.4375   -2.5007    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -5.8004   -3.1298    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.4375   -1.7504    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -9.0343   -3.9223    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -9.0343   -4.7460    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -9.7468   -3.5107    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3787   -3.0670    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7831   -4.1028    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.6814   -3.0670    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5323   -3.3982    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0416   -3.8144    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9308   -3.7745    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0887   -4.1028    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4989   -3.3027    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4931   -6.8542    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0887   -5.8183    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1903   -6.8542    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3394   -6.5230    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5173   -6.4300    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.0190   -5.9999    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9410   -6.9704    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0887   -4.9945    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6553   -6.9835    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9410   -6.1466    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7831   -5.8183    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9775   -7.6505    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5009   -4.5189    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0192   -1.1155    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6148   -2.1511    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.2836   -1.1155    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8656   -1.4465    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.9912   -1.5394    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.4929   -1.9694    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4671   -0.9993    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0917   -2.8245    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1813   -0.9861    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4671   -1.8228    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3092   -2.1511    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4144    1.2848    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1813    0.2490    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.7170    1.2848    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5679    0.9537    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.4246    0.8607    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.9263    0.4306    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.9664    1.4010    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1813   -0.4098    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.2521    1.4142    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.9664    0.5773    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.1243    0.2490    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.2521    2.0109    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.7232   -2.9131    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.8481    3.7391    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.2524    2.7033    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.1507    3.7391    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.0016    3.4080    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.8583    3.3150    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.3600    2.8849    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.4001    3.8553    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.6858    3.8685    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.4001    3.0316    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.5580    2.7033    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.6857    4.3828    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.2818    6.0491    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.6861    5.0133    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.5844    6.0491    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.4353    5.7180    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.2920    5.6250    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.7937    5.1949    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.8338    6.1653    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.1195    6.1785    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.8338    5.3416    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.9917    5.0133    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.1195    6.6722    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.7155    8.3543    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.1198    7.3185    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.0182    8.3543    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.8691    8.0232    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.7257    7.9302    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   10.2275    7.5001    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    8.2675    8.4705    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    9.5533    8.4837    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.2675    7.6468    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    9.4255    7.3185    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    9.2851    8.8743    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5781   -8.1318    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5781   -8.8743    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.3206   -7.6368    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  2  0      
  4  8  1  0      
  1  7  1  6      
  1  6  1  1      
  9  1  1  0      
 10  9  1  0      
 11 12  1  0      
 11 13  2  0      
 14 15  2  0      
 14 16  1  0      
 14  7  1  0      
 11 10  1  0      
 17 18  1  0      
 18 22  1  1      
 19 23  1  1      
 19 20  1  0      
 17 20  1  0      
 18 21  1  0      
 23 22  1  1      
 20 45  1  0      
 17 24  1  0      
 22 37  1  0      
 23 32  1  0      
 25 26  1  0      
 26 34  1  1      
 27 35  1  1      
 27 28  1  0      
 25 28  1  0      
 27 29  1  0      
 30 35  1  0      
 31 34  1  0      
 26 32  1  0      
 28 33  1  0      
 34 35  1  1      
 33 36  1  0      
 38 39  1  0      
 39 47  1  1      
 40 48  1  1      
 40 41  1  0      
 38 41  1  0      
 40 42  1  0      
 43 48  1  0      
 44 47  1  0      
 39 45  1  0      
 41 46  1  0      
 47 48  1  1      
 46 56  1  0      
 49 50  1  0      
 50 58  1  1      
 51 59  1  1      
 51 52  1  0      
 49 52  1  0      
 51 53  1  0      
 54 59  1  0      
 55 58  1  0      
 50 56  1  0      
 52 57  1  0      
 58 59  1  1      
 57 60  1  0      
 11 61  1  0      
 61 19  1  0      
 62 63  1  0      
 63 70  1  1      
 64 71  1  1      
 64 65  1  0      
 62 65  1  0      
 64 66  1  0      
 67 71  1  0      
 68 70  1  0      
 65 69  1  0      
 70 71  1  1      
 69 72  1  0      
 60 63  1  0      
 73 74  1  0      
 74 81  1  1      
 75 82  1  1      
 75 76  1  0      
 73 76  1  0      
 75 77  1  0      
 78 82  1  0      
 79 81  1  0      
 76 80  1  0      
 81 82  1  1      
 80 83  1  0      
 74 72  1  0      
 84 85  1  0      
 85 92  1  1      
 86 93  1  1      
 86 87  1  0      
 84 87  1  0      
 86 88  1  0      
 89 93  1  0      
 90 92  1  0      
 87 91  1  0      
 92 93  1  1      
 91 94  1  0      
 83 85  1  0      
 36 95  1  0      
 95 96  2  0      
 95 97  1  0      
M  END

ENDAC1PIM6TEMPLATE

  my($CLSn1Sn2Sn3Sn4TemplateString)=<<ENDTEMPLATE;
CL sn1 acyl, sn2, sn3 and sn4 acyl template structure
  LipdMAPS06110711152D

 39 38  0  0  0  0  0  0  0  0999 V2000
    2.5696    2.2606    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    2.5696    3.0855    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.3568    1.4635    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.3667    2.0477    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.9378    0.8810    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.7267    0.0855    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.1391   -0.6280    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8620    1.3504    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1371    0.1178    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1371    0.9379    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4170    1.3504    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.8867    2.0639    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1664    3.2997    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1664    2.4763    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4365    2.0639    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0093    2.0639    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7213    2.4763    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4112    1.3504    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.5238    0.2984    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.9297    0.2984    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.6442   -1.5531    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.8683   -1.1975    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    1.8683   -0.3891    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.6938    2.4294    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5773   -1.9413    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.3927    2.0258    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.7630   -1.3477    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.0486   -0.9352    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6659   -1.3477    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3804   -0.9352    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0949   -1.3477    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.8093   -0.9352    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.5238   -1.3477    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.8093   -0.1102    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0784   -2.0622    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7929   -2.4747    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5074   -2.0622    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7929   -3.2997    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2534   -2.0622    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0      
  1  3  1  0      
  1  4  1  0      
  3  5  1  0      
  5  6  1  0      
  6  7  1  0      
  9 10  2  0      
  8 10  1  0      
 10 11  1  0      
 13 14  2  0      
 12 14  1  0      
 14 15  1  0      
 15 17  1  0      
 16 17  1  0      
  6 19  1  1      
  6 20  1  1      
  7 21  1  0      
 21 22  1  0      
 22 23  2  0      
 22 25  1  0      
 24 16  1  0      
 24 26  1  0      
 26  1  1  0      
 27 22  1  0      
 16 11  1  6      
 16 18  1  1      
 27 28  1  0      
 29 28  1  0      
 30 29  1  0      
 31 30  1  0      
 32 31  1  0      
 33 32  1  0      
 32 34  2  0      
 29 35  1  6      
 36 35  1  0      
 37 36  1  0      
 36 38  2  0      
 29 39  1  1      
M  END

ENDTEMPLATE

  my($AC2PIM2ReverseSn2TemplateString)=<<ENDAC2PIM2TEMPLATE;
Ac2PIM2.mol
  ChemDraw01251109492D

 56 58  0  0  0  0  0  0  0  0999 V2000
   -3.0508    1.6806    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.7632    2.0907    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.4760    1.6806    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.1880    2.0907    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.1880    2.9142    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6388    0.9682    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -3.4625    0.9682    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.9005    1.6806    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.3380    2.0920    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6256    1.6806    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6117    1.9690    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9746    1.3397    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6117    2.7192    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.2089    0.5470    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.2089   -0.2769    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.9214    0.9587    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.4479    1.4025    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.0436    0.3665    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.1449    1.4025    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.2940    1.0712    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.7852    0.6550    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.8957    0.6949    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.7376    0.3665    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.3279    1.1667    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.3333   -2.3853    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.7376   -1.3492    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.6365   -2.3853    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.4872   -2.0541    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.3440   -1.9610    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.8457   -1.5309    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.8855   -2.5015    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.7376   -0.5255    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.1714   -2.5146    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.8855   -1.6776    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.0436   -1.3492    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.8492   -3.1818    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.3257   -0.0498    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.8898    3.6637    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.2940    2.6278    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.1928    3.6637    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.0436    3.3325    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.9005    3.2395    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.4021    2.8095    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.4420    3.7799    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.7347    1.6450    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.7280    3.7931    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.4420    2.9561    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.6001    2.6278    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.4374    4.3533    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.1029    1.5564    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.4885   -3.5874    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.4993   -4.3533    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.1385   -3.2188    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.0141   -0.3286    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.0141   -1.0299    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.5297    0.0632    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  2  0      
  4  8  1  0      
  1  7  1  6      
  1  6  1  1      
  9  1  1  0      
 10  9  1  0      
 11 12  1  0      
 11 13  2  0      
 14 15  2  0      
 14 16  1  0      
 14  7  1  0      
 11 10  1  0      
 17 18  1  0      
 18 22  1  1      
 19 23  1  1      
 19 20  1  0      
 17 20  1  0      
 18 21  1  0      
 23 22  1  1      
 20 45  1  0      
 17 24  1  0      
 22 37  1  0      
 23 32  1  0      
 25 26  1  0      
 26 34  1  1      
 27 35  1  1      
 27 28  1  0      
 25 28  1  0      
 27 29  1  0      
 30 35  1  0      
 31 34  1  0      
 26 32  1  0      
 28 33  1  0      
 34 35  1  1      
 33 36  1  0      
 38 39  1  0      
 39 47  1  1      
 40 48  1  1      
 40 41  1  0      
 38 41  1  0      
 40 42  1  0      
 43 48  1  0      
 44 47  1  0      
 39 45  1  0      
 41 46  1  0      
 47 48  1  1      
 46 49  1  0      
 11 50  1  0      
 50 19  1  0      
 36 51  1  0      
 51 52  2  0      
 51 53  1  0      
 37 54  1  0      
 54 55  2  0      
 54 56  1  0      
M  END

ENDAC2PIM2TEMPLATE

  my($AC2PIM3ReverseSn2TemplateString)=<<ENDAC2PIM3TEMPLATE;
Ac2PIM2.mol
  ChemDraw01251113122D

 67 70  0  0  0  0  0  0  0  0999 V2000
   -4.2590    0.4428    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.9714    0.8529    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.6842    0.4428    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -6.3962    0.8529    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.3962    1.6764    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.8471   -0.2696    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -4.6708   -0.2696    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -7.1088    0.4428    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.5462    0.8542    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.8339    0.4428    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8200    0.7311    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1829    0.1018    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8200    1.4813    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.4172   -0.6908    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.4172   -1.5147    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -6.1297   -0.2791    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.2397    0.1646    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.8353   -0.8713    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0634    0.1646    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0858   -0.1666    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.5770   -0.5829    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.6875   -0.5430    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5293   -0.8713    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.1196   -0.0711    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.1250   -3.6231    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.5293   -2.5871    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.4282   -3.6231    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.2790   -3.2919    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.1357   -3.1989    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.6375   -2.7687    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.6773   -3.7393    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.5293   -1.7633    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.9632   -3.7524    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.6773   -2.9154    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.8353   -2.5871    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.6409   -4.4196    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.1174   -1.2876    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.6816    2.4259    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.0858    1.3899    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.9845    2.4259    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.8353    2.0946    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.6922    2.0017    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.1939    1.5717    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.2338    2.5421    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.5264    0.4072    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.5198    2.5553    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.2338    1.7183    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.3918    1.3899    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.4972    3.2805    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1053    0.3186    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.2802   -4.8252    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.2910   -5.5911    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.9303   -4.4566    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.8058   -1.5664    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.8058   -2.2677    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.3215   -1.1746    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.0982    5.0711    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.5024    4.0353    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.4011    5.0711    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.2519    4.7400    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.1088    4.6469    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    6.6104    4.2169    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.6502    5.1873    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.9362    5.2005    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.6502    4.3636    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.8083    4.0353    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.6679    5.5911    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  2  0      
  4  8  1  0      
  1  7  1  6      
  1  6  1  1      
  9  1  1  0      
 10  9  1  0      
 11 12  1  0      
 11 13  2  0      
 14 15  2  0      
 14 16  1  0      
 14  7  1  0      
 11 10  1  0      
 17 18  1  0      
 18 22  1  1      
 19 23  1  1      
 19 20  1  0      
 17 20  1  0      
 18 21  1  0      
 23 22  1  1      
 20 45  1  0      
 17 24  1  0      
 22 37  1  0      
 23 32  1  0      
 25 26  1  0      
 26 34  1  1      
 27 35  1  1      
 27 28  1  0      
 25 28  1  0      
 27 29  1  0      
 30 35  1  0      
 31 34  1  0      
 26 32  1  0      
 28 33  1  0      
 34 35  1  1      
 33 36  1  0      
 38 39  1  0      
 39 47  1  1      
 40 48  1  1      
 40 41  1  0      
 38 41  1  0      
 40 42  1  0      
 43 48  1  0      
 44 47  1  0      
 39 45  1  0      
 41 46  1  0      
 47 48  1  1      
 46 49  1  0      
 11 50  1  0      
 50 19  1  0      
 36 51  1  0      
 51 52  2  0      
 51 53  1  0      
 37 54  1  0      
 54 55  2  0      
 54 56  1  0      
 57 58  1  0      
 58 65  1  1      
 59 66  1  1      
 59 60  1  0      
 57 60  1  0      
 59 61  1  0      
 62 66  1  0      
 63 65  1  0      
 60 64  1  0      
 65 66  1  1      
 64 67  1  0      
 49 58  1  0      
M  END

ENDAC2PIM3TEMPLATE

  my($AC2PIM4ReverseSn2TemplateString)=<<ENDAC2PIM4TEMPLATE;
Ac2PIM3.mol
  ChemDraw01251113242D

 78 82  0  0  0  0  0  0  0  0999 V2000
   -5.4774   -0.7637    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.1898   -0.3536    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.9026   -0.7637    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -7.6147   -0.3536    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.6147    0.4699    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.0655   -1.4761    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -5.8892   -1.4761    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -8.3272   -0.7637    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.7646   -0.3523    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.0523   -0.7637    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0384   -0.4753    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -3.4013   -1.1046    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0384    0.2749    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -6.6356   -1.8972    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.6356   -2.7212    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -7.3481   -1.4856    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0213   -1.0418    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.6169   -2.0778    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2818   -1.0418    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1326   -1.3731    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.3586   -1.7893    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.4691   -1.7494    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6891   -2.0778    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.9012   -1.2775    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0934   -4.8296    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6891   -3.7935    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.2098   -4.8296    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0606   -4.4983    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.9173   -4.4053    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.4191   -3.9751    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.4589   -4.9458    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6891   -2.9698    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.7448   -4.9589    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4589   -4.1219    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.6169   -3.7935    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4225   -5.6261    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.8990   -2.4940    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.4632    1.2194    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1326    0.1835    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.7661    1.2194    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.6169    0.8882    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.4738    0.7953    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.9754    0.3652    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.0154    1.3357    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.3080   -0.7993    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.3013    1.3489    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0154    0.5118    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.1734    0.1835    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.2788    2.0740    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.3238   -0.8879    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.0618   -6.0316    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0726   -6.7976    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.7119   -5.6631    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5874   -2.7729    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5874   -3.4741    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.1031   -2.3810    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.8798    3.8646    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.2840    2.8288    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.1827    3.8646    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.0335    3.5335    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.8904    3.4405    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.3920    3.0105    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.4318    3.9808    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.7177    3.9940    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.4318    3.1572    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.5899    2.8288    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.7176    4.5497    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.3168    6.2776    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.7211    5.2418    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.6195    6.2776    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.4704    5.9464    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.3272    5.8534    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.8288    5.4233    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.8688    6.3938    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.1546    6.4070    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.8688    5.5700    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.0268    5.2418    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.8864    6.7976    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  2  0      
  4  8  1  0      
  1  7  1  6      
  1  6  1  1      
  9  1  1  0      
 10  9  1  0      
 11 12  1  0      
 11 13  2  0      
 14 15  2  0      
 14 16  1  0      
 14  7  1  0      
 11 10  1  0      
 17 18  1  0      
 18 22  1  1      
 19 23  1  1      
 19 20  1  0      
 17 20  1  0      
 18 21  1  0      
 23 22  1  1      
 20 45  1  0      
 17 24  1  0      
 22 37  1  0      
 23 32  1  0      
 25 26  1  0      
 26 34  1  1      
 27 35  1  1      
 27 28  1  0      
 25 28  1  0      
 27 29  1  0      
 30 35  1  0      
 31 34  1  0      
 26 32  1  0      
 28 33  1  0      
 34 35  1  1      
 33 36  1  0      
 38 39  1  0      
 39 47  1  1      
 40 48  1  1      
 40 41  1  0      
 38 41  1  0      
 40 42  1  0      
 43 48  1  0      
 44 47  1  0      
 39 45  1  0      
 41 46  1  0      
 47 48  1  1      
 46 49  1  0      
 11 50  1  0      
 50 19  1  0      
 36 51  1  0      
 51 52  2  0      
 51 53  1  0      
 37 54  1  0      
 54 55  2  0      
 54 56  1  0      
 57 58  1  0      
 58 65  1  1      
 59 66  1  1      
 59 60  1  0      
 57 60  1  0      
 59 61  1  0      
 62 66  1  0      
 63 65  1  0      
 60 64  1  0      
 65 66  1  1      
 64 67  1  0      
 49 58  1  0      
 68 69  1  0      
 69 76  1  1      
 70 77  1  1      
 70 71  1  0      
 68 71  1  0      
 70 72  1  0      
 73 77  1  0      
 74 76  1  0      
 71 75  1  0      
 76 77  1  1      
 75 78  1  0      
 67 69  1  0      
M  END

ENDAC2PIM4TEMPLATE

  my($AC2PIM5ReverseSn2TemplateString)=<<ENDAC2PIM5TEMPLATE;
Ac2PIM4.mol
  ChemDraw01251113352D

 89 94  0  0  0  0  0  0  0  0999 V2000
   -6.6943   -1.9599    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.4067   -1.5498    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -8.1195   -1.9599    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -8.8315   -1.5498    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -8.8315   -0.7263    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -6.2823   -2.6723    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -7.1061   -2.6723    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -9.5440   -1.9599    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.9815   -1.5485    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.2692   -1.9599    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.2552   -1.6716    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -4.6182   -2.3009    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.2552   -0.9214    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -7.8525   -3.0935    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.8525   -3.9174    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -8.5650   -2.6818    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1956   -2.2381    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4001   -3.2740    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4987   -2.2381    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3495   -2.5693    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.1417   -2.9856    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7478   -2.9457    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9060   -3.2740    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6843   -2.4738    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3103   -6.0258    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9060   -4.9898    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.9930   -6.0258    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1563   -5.6946    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.7005   -5.6016    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.2022   -5.1714    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7580   -6.1420    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9060   -4.1660    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.5279   -6.1551    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7580   -5.3181    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4001   -4.9898    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.2056   -6.8223    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3178   -3.6903    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7537    0.0232    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3495   -1.0128    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5492    0.0232    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4001   -0.3081    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.2569   -0.4010    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.7586   -0.8310    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2015    0.1394    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9088   -1.9955    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.0845    0.1526    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2015   -0.6844    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.9566   -1.0128    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0620    0.8778    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.5406   -2.0841    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.8449   -7.2279    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.8557   -7.9938    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.4950   -6.8593    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3706   -3.9691    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3706   -4.6704    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.8862   -3.5773    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.6629    2.6684    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.0672    1.6326    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.9658    2.6684    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.8166    2.3373    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.6735    2.2442    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.1751    1.8142    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.2150    2.7846    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.5009    2.7978    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.2150    1.9609    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.3731    1.6326    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.5008    3.3534    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.0999    5.0813    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.5042    4.0456    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.4027    5.0813    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.2535    4.7502    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.1103    4.6572    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    6.6120    4.2271    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.6520    5.1975    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.9377    5.2107    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.6520    4.3738    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.8099    4.0456    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.9376    5.7663    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    6.5337    7.4738    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.9380    6.4381    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.8364    7.4738    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.6872    7.1427    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    9.5440    7.0497    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    9.0457    6.6196    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.0857    7.5900    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    8.3715    7.6032    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.0857    6.7663    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.2437    6.4381    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.1033    7.9938    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  2  0      
  4  8  1  0      
  1  7  1  6      
  1  6  1  1      
  9  1  1  0      
 10  9  1  0      
 11 12  1  0      
 11 13  2  0      
 14 15  2  0      
 14 16  1  0      
 14  7  1  0      
 11 10  1  0      
 17 18  1  0      
 18 22  1  1      
 19 23  1  1      
 19 20  1  0      
 17 20  1  0      
 18 21  1  0      
 23 22  1  1      
 20 45  1  0      
 17 24  1  0      
 22 37  1  0      
 23 32  1  0      
 25 26  1  0      
 26 34  1  1      
 27 35  1  1      
 27 28  1  0      
 25 28  1  0      
 27 29  1  0      
 30 35  1  0      
 31 34  1  0      
 26 32  1  0      
 28 33  1  0      
 34 35  1  1      
 33 36  1  0      
 38 39  1  0      
 39 47  1  1      
 40 48  1  1      
 40 41  1  0      
 38 41  1  0      
 40 42  1  0      
 43 48  1  0      
 44 47  1  0      
 39 45  1  0      
 41 46  1  0      
 47 48  1  1      
 46 49  1  0      
 11 50  1  0      
 50 19  1  0      
 36 51  1  0      
 51 52  2  0      
 51 53  1  0      
 37 54  1  0      
 54 55  2  0      
 54 56  1  0      
 57 58  1  0      
 58 65  1  1      
 59 66  1  1      
 59 60  1  0      
 57 60  1  0      
 59 61  1  0      
 62 66  1  0      
 63 65  1  0      
 60 64  1  0      
 65 66  1  1      
 64 67  1  0      
 49 58  1  0      
 68 69  1  0      
 69 76  1  1      
 70 77  1  1      
 70 71  1  0      
 68 71  1  0      
 70 72  1  0      
 73 77  1  0      
 74 76  1  0      
 71 75  1  0      
 76 77  1  1      
 75 78  1  0      
 67 69  1  0      
 79 80  1  0      
 80 87  1  1      
 81 88  1  1      
 81 82  1  0      
 79 82  1  0      
 81 83  1  0      
 84 88  1  0      
 85 87  1  0      
 82 86  1  0      
 87 88  1  1      
 86 89  1  0      
 78 80  1  0      
M  END

ENDAC2PIM5TEMPLATE

  my($AC2PIM6ReverseSn2TemplateString)=<<ENDAC2PIM6TEMPLATE;
Ac2PIM6.mol
  ChemDraw01251113542D

100106  0  0  0  0  0  0  0  0999 V2000
   -7.9112   -3.0892    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -8.6236   -2.6791    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -9.3364   -3.0892    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  -10.0484   -2.6791    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  -10.0484   -1.8556    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -7.4992   -3.8016    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -8.3229   -3.8016    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  -10.7609   -3.0892    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.1984   -2.6778    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.4861   -3.0892    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.4721   -2.8009    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -5.8351   -3.4302    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.4721   -2.0506    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -9.0694   -4.2228    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -9.0694   -5.0467    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -9.7819   -3.8111    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4125   -3.3673    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8168   -4.4033    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.7155   -3.3673    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5664   -3.6986    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0752   -4.1149    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9647   -4.0750    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1229   -4.4033    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5326   -3.6031    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5272   -7.2995    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1229   -6.2634    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2239   -7.2995    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3732   -6.9682    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4836   -6.8752    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0147   -6.4451    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9749   -7.4157    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1229   -5.2953    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6890   -7.4288    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9749   -6.5918    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8168   -6.2634    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0112   -8.0960    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5347   -4.8196    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9706   -1.1061    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5664   -2.1421    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3324   -1.1061    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8168   -1.4374    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0401   -1.5303    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.5417   -1.9603    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4184   -0.9899    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1257   -3.1248    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1324   -0.9767    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4184   -1.8137    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2603   -2.1421    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1549   -0.2515    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.7575   -3.2134    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3719   -8.5016    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3612   -9.2675    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.2781   -8.1330    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8463   -5.0984    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8463   -5.7997    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2069   -4.7478    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4460    1.5391    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1497    0.5033    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.7489    1.5391    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5997    1.2080    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.4566    1.1150    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.9583    0.6849    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.9981    1.6553    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.2840    1.6685    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.9981    0.8316    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.1562    0.5033    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.2839    2.2241    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.8831    3.9520    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.2873    2.9163    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.1858    3.9520    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.0366    3.6209    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.8934    3.5279    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.3951    3.0978    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.4351    4.0682    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.7209    4.0814    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.4351    3.2445    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.5930    2.9163    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.7208    4.6370    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.3168    6.3445    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.7211    5.3088    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.6195    6.3445    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.4704    6.0134    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.3272    5.9204    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.8288    5.4903    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.8688    6.4607    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.1546    6.4739    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.8688    5.6370    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.0268    5.3088    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.1545    7.0502    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.7506    8.7475    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.1548    7.7117    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.0533    8.7475    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.9041    8.4164    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.7609    8.3234    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   10.2626    7.8932    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    8.3026    8.8637    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    9.5884    8.8769    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.3026    8.0399    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    9.4605    7.7117    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    9.3201    9.2675    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  2  0      
  4  8  1  0      
  1  7  1  6      
  1  6  1  1      
  9  1  1  0      
 10  9  1  0      
 11 12  1  0      
 11 13  2  0      
 14 15  2  0      
 14 16  1  0      
 14  7  1  0      
 11 10  1  0      
 17 18  1  0      
 18 22  1  1      
 19 23  1  1      
 19 20  1  0      
 17 20  1  0      
 18 21  1  0      
 23 22  1  1      
 20 45  1  0      
 17 24  1  0      
 22 37  1  0      
 23 32  1  0      
 25 26  1  0      
 26 34  1  1      
 27 35  1  1      
 27 28  1  0      
 25 28  1  0      
 27 29  1  0      
 30 35  1  0      
 31 34  1  0      
 26 32  1  0      
 28 33  1  0      
 34 35  1  1      
 33 36  1  0      
 38 39  1  0      
 39 47  1  1      
 40 48  1  1      
 40 41  1  0      
 38 41  1  0      
 40 42  1  0      
 43 48  1  0      
 44 47  1  0      
 39 45  1  0      
 41 46  1  0      
 47 48  1  1      
 46 49  1  0      
 11 50  1  0      
 50 19  1  0      
 36 51  1  0      
 51 52  2  0      
 51 53  1  0      
 37 54  1  0      
 54 55  2  0      
 54 56  1  0      
 57 58  1  0      
 58 65  1  1      
 59 66  1  1      
 59 60  1  0      
 57 60  1  0      
 59 61  1  0      
 62 66  1  0      
 63 65  1  0      
 60 64  1  0      
 65 66  1  1      
 64 67  1  0      
 49 58  1  0      
 68 69  1  0      
 69 76  1  1      
 70 77  1  1      
 70 71  1  0      
 68 71  1  0      
 70 72  1  0      
 73 77  1  0      
 74 76  1  0      
 71 75  1  0      
 76 77  1  1      
 75 78  1  0      
 67 69  1  0      
 79 80  1  0      
 80 87  1  1      
 81 88  1  1      
 81 82  1  0      
 79 82  1  0      
 81 83  1  0      
 84 88  1  0      
 85 87  1  0      
 82 86  1  0      
 87 88  1  1      
 86 89  1  0      
 78 80  1  0      
 90 91  1  0      
 91 98  1  1      
 92 99  1  1      
 92 93  1  0      
 90 93  1  0      
 92 94  1  0      
 95 99  1  0      
 96 98  1  0      
 93 97  1  0      
 98 99  1  1      
 97100  1  0      
 89 91  1  0      
M  END

ENDAC2PIM6TEMPLATE

  my($SL1TemplateString)=<<SL1TEMPLATE;
SL1.mol
  ChemDraw01251116442D

 35 36  0  0  0  0  0  0  0  0999 V2000
   -0.9684    1.8976    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3727    0.8617    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.2711    1.8976    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1220    1.5664    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.9787    1.4734    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.4804    1.0434    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8000    0.6337    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.8062    2.0270    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5203    1.1901    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6784    0.8617    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.8268    2.8622    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.2230   -1.8217    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3727   -0.7858    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.5257   -1.8217    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.3766   -1.4904    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.2333   -1.3975    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.7351   -0.9674    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3727    0.0379    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.7751   -1.1142    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.9331   -0.7858    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5943   -1.7186    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.3522   -0.6198    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3698   -1.2777    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3698   -0.5566    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.2051   -2.8196    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2295    3.2897    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2295    3.9909    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.9787   -0.0928    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.9787    0.6291    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.9680   -3.3309    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.9680   -3.9909    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2501    0.1959    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8995   -0.1341    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.7657   -0.1547    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6626    0.5878    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  9  1  1      
  3 10  1  1      
  3  4  1  0      
  1  4  1  0      
  3  5  1  0      
  6 10  1  0      
  2 18  1  0      
  7  9  1  0      
  4  8  1  0      
  9 10  1  1      
  8 11  1  0      
 13 12  1  1      
 13 19  1  0      
 14 20  1  0      
 14 15  1  1      
 12 15  1  1      
 14 16  1  0      
 17 20  1  0      
 13 18  1  0      
 19 20  1  0      
 12 21  1  0      
 15 25  1  0      
 17 22  1  0      
 21 23  1  0      
 23 24  2  0      
 11 26  1  0      
 26 27  2  0      
 22 28  1  0      
 28 29  2  0      
 25 30  1  0      
 30 31  2  0      
  7 32  1  0      
 32 33  2  0      
 32 34  1  0      
 32 35  2  0      
M  END

SL1TEMPLATE

  my($TESTTemplateString)=<<TEST;
TEST.mol
  ChemDraw01251117112D

 20 19  0  0  0  0  0  0  0  0999 V2000
    0.2888    4.1261    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1444    3.3759    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.2888    2.6257    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1444    1.8755    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.2888    1.1253    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1444    0.3751    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.2888   -0.3751    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1444   -1.1253    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.2888   -1.8755    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1444   -2.6257    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.2888   -3.3759    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1444   -4.1261    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0106    3.0742    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0106    1.5479    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9281   -0.6383    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0313   -2.1233    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0313    3.8167    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.0313    2.2492    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9281    0.0423    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0313   -1.3808    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  1  0      
  5  6  1  0      
  6  7  1  0      
  7  8  1  0      
  8  9  1  0      
  9 10  1  0      
 10 11  1  0      
 11 12  1  0      
  3 13  1  0      
  5 14  1  0      
  8 15  1  0      
 10 16  1  0      
 13 17  2  0      
 14 18  2  0      
 15 19  2  0      
 16 20  2  0      
M  END

TEST

  my($SL2TemplateString)=<<SL2PRIMETEMPLATE;
SL2prime.mol
  ChemDraw01261108482D

 35 36  0  0  0  0  0  0  0  0999 V2000
   -0.9684    1.6185    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3727    0.5826    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.2711    1.6185    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1220    1.2873    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.9786    1.1943    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.4804    0.7643    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8000    0.3546    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.8062    1.7479    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5203    0.9110    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6784    0.5826    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.8268    2.5831    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.2230   -2.1008    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3727   -1.0649    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.5257   -2.1008    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.3766   -1.7695    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.7258   -1.1130    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3727   -0.2412    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.7751   -1.3933    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.9331   -1.0649    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5943   -1.9977    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.3315   -0.7339    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3698   -1.5568    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3698   -0.8357    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.9676   -2.5005    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2295    3.0106    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2295    3.7118    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.9786   -0.3719    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.9786    0.3500    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2501   -0.0832    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8995   -0.4132    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.7657   -0.4338    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6626    0.3087    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.2208   -2.4330    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.8396   -2.9486    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.8396   -3.7118    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  9  1  1      
  3 10  1  1      
  3  4  1  0      
  1  4  1  0      
  3  5  1  0      
  6 10  1  0      
  2 17  1  0      
  7  9  1  0      
  4  8  1  0      
  9 10  1  1      
  8 11  1  0      
 13 12  1  1      
 13 18  1  0      
 14 19  1  0      
 14 15  1  1      
 12 15  1  1      
 16 19  1  0      
 13 17  1  0      
 18 19  1  0      
 12 20  1  0      
 15 24  1  0      
 16 21  1  0      
 20 22  1  0      
 22 23  2  0      
 11 25  1  0      
 25 26  2  0      
 21 27  1  0      
 27 28  2  0      
  7 29  1  0      
 29 30  2  0      
 29 31  1  0      
 29 32  2  0      
 14 33  1  0      
 33 34  1  0      
 34 35  2  0      
M  END

SL2PRIMETEMPLATE

  my($SL3TemplateString)=<<SL3TEMPLATE;
SL3.mol
  ChemDraw01261108492D

 33 34  0  0  0  0  0  0  0  0999 V2000
   -0.9684    2.6269    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3727    1.5910    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.2711    2.6269    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1220    2.2957    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.9786    2.2027    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.4804    1.7727    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8000    1.3630    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.8062    2.7563    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5203    1.9194    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6784    1.5910    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.3111    3.2615    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.2230   -1.0923    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3727   -0.0564    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.5257   -1.0923    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.3766   -0.7610    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.2333   -0.6681    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.7351   -0.2380    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3727    0.7672    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.7751   -0.3848    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.9331   -0.0564    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5943   -0.9892    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.3522    0.1096    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3698   -0.5483    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3698    0.1727    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.2051   -2.0902    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.9786    0.6365    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.9786    1.3584    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.9680   -2.6015    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.9680   -3.2615    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2501    0.9252    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8995    0.5952    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.7657    0.5746    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6626    1.3171    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  9  1  1      
  3 10  1  1      
  3  4  1  0      
  1  4  1  0      
  3  5  1  0      
  6 10  1  0      
  2 18  1  0      
  7  9  1  0      
  4  8  1  0      
  9 10  1  1      
  8 11  1  0      
 13 12  1  1      
 13 19  1  0      
 14 20  1  0      
 14 15  1  1      
 12 15  1  1      
 14 16  1  0      
 17 20  1  0      
 13 18  1  0      
 19 20  1  0      
 12 21  1  0      
 15 25  1  0      
 17 22  1  0      
 21 23  1  0      
 23 24  2  0      
 22 26  1  0      
 26 27  2  0      
 25 28  1  0      
 28 29  2  0      
  7 30  1  0      
 30 31  2  0      
 30 32  1  0      
 30 33  2  0      
M  END

SL3TEMPLATE

  my($PAT18ReverseSn2TemplateString)=<<PAT18TEMPLATE;
PAT18.mol
  ChemDraw02011117242D

 50 51  0  0  0  0  0  0  0  0999 V2000
    3.4391    2.6019    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.9504    1.5201    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.1408    2.7872    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.2610    2.3641    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0626    2.9488    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.0903    2.0547    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.4678    0.6666    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.6151    2.8789    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.8313    1.9401    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.6489    1.7055    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.6624    3.4655    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.7319   -1.9168    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.1354   -0.8796    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.0378   -1.9168    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.8871   -1.5851    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.7464   -1.4921    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.2474   -1.0614    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.9290    0.3584    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.2848   -1.2085    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    6.4443   -0.8796    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.9134   -1.8136    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.8654   -0.3210    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.1369   -1.3721    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.1369   -0.7120    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.4330   -1.7785    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.8406   -3.0607    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.5220   -3.5314    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.5321   -4.3246    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.7489   -1.3835    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0648   -1.7785    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3806   -1.3835    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3036   -1.7785    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9877   -1.3835    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6719   -1.7785    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.3560   -1.3835    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0401   -1.7785    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.7243   -1.3835    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.4085   -1.7785    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.0926   -1.3835    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.7768   -1.7785    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.4609   -1.3835    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.1450   -1.7785    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.8292   -1.3835    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -8.5134   -1.7785    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6616    3.4491    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6631    4.3246    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.6456    0.2074    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6456   -0.6187    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    8.5134    0.1655    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.5134    0.9297    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  9  1  1      
  3 10  1  1      
  3  4  1  0      
  1  4  1  0      
  3  5  1  0      
  6 10  1  0      
  2 18  1  0      
  7  9  1  0      
  4  8  1  0      
  9 10  1  1      
  8 11  1  0      
 13 12  1  1      
 13 19  1  0      
 14 20  1  0      
 14 15  1  1      
 12 15  1  1      
 14 16  1  0      
 17 20  1  0      
 13 18  1  0      
 19 20  1  0      
 12 21  1  0      
 15 26  1  0      
 17 22  1  0      
 21 23  1  0      
 23 24  2  0      
 23 25  1  0      
 26 27  1  0      
 27 28  2  0      
 25 29  1  0      
 29 30  1  0      
 30 31  1  0      
 31 32  1  0      
 32 33  1  0      
 33 34  1  0      
 34 35  1  0      
 35 36  1  0      
 36 37  1  0      
 37 38  1  0      
 38 39  1  0      
 39 40  1  0      
 40 41  1  0      
 41 42  1  0      
 42 43  1  0      
 43 44  1  0      
  5 45  1  0      
 45 46  2  0      
  7 47  1  0      
 47 48  2  0      
 22 49  1  0      
 49 50  2  0      
M  END

PAT18TEMPLATE

  my($PAT16ReverseSn2TemplateString)=<<PAT16TEMPLATE;
PAT18.mol
  ChemDraw02011117242D

 48 49  0  0  0  0  0  0  0  0999 V2000
    2.7549    2.6019    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.2662    1.5201    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4566    2.7872    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5768    2.3641    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6216    2.9488    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.4061    2.0547    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.7836    0.6666    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.9309    2.8789    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.1471    1.9401    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.9648    1.7055    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.9782    3.4655    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.0478   -1.9168    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.4513   -0.8796    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.3536   -1.9168    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.2029   -1.5851    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.0622   -1.4921    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    6.5632   -1.0614    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.2448    0.3584    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.6006   -1.2085    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.7601   -0.8796    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.2293   -1.8136    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.1812   -0.3210    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.4527   -1.3721    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.4527   -0.7120    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.7488   -1.7785    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.1564   -3.0607    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    6.8379   -3.5314    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.8479   -4.3246    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.0647   -1.3835    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3806   -1.7785    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3036   -1.3835    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9878   -1.7785    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6719   -1.3835    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.3561   -1.7785    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0402   -1.3835    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.7243   -1.7785    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.4085   -1.3835    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.0927   -1.7785    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.7768   -1.3835    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.4610   -1.7785    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.1451   -1.3835    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.8292   -1.7785    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3457    3.4491    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3473    4.3246    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0386    0.2074    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0386   -0.6187    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.8292    0.1655    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.8292    0.9297    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  9  1  1      
  3 10  1  1      
  3  4  1  0      
  1  4  1  0      
  3  5  1  0      
  6 10  1  0      
  2 18  1  0      
  7  9  1  0      
  4  8  1  0      
  9 10  1  1      
  8 11  1  0      
 13 12  1  1      
 13 19  1  0      
 14 20  1  0      
 14 15  1  1      
 12 15  1  1      
 14 16  1  0      
 17 20  1  0      
 13 18  1  0      
 19 20  1  0      
 12 21  1  0      
 15 26  1  0      
 17 22  1  0      
 21 23  1  0      
 23 24  2  0      
 23 25  1  0      
 26 27  1  0      
 27 28  2  0      
 25 29  1  0      
 29 30  1  0      
 30 31  1  0      
 31 32  1  0      
 32 33  1  0      
 33 34  1  0      
 34 35  1  0      
 35 36  1  0      
 36 37  1  0      
 37 38  1  0      
 38 39  1  0      
 39 40  1  0      
 40 41  1  0      
 41 42  1  0      
  5 43  1  0      
 43 44  2  0      
  7 45  1  0      
 45 46  2  0      
 22 47  1  0      
 47 48  2  0      
M  END

PAT16TEMPLATE

  my($ASNSn1TemplateString)=<<ASNTEMPLATE;
ASN.mol
  ChemDraw03111113532D

 12 11  0  0  0  0  0  0  0  0999 V2000
    0.3505   -1.5352    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0767   -0.4869    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3039    0.5272    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.9483    1.6027    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.2966    2.2919    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.8480    1.6104    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7908   -0.0656    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4189    1.2694    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8480    0.6518    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8480    1.7148    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.0957   -2.2919    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6759   -1.5287    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  1  0      
  4  6  2  0      
  3  7  1  6      
  3  8  1  1      
  7  9  1  0      
  9 10  2  0      
  1 11  1  0      
  1 12  2  0      
M  END

ASNTEMPLATE

  my($HISSn1TemplateString)=<<HISTEMPLATE;
HIS.mol
  ChemDraw03111113532D

 14 14  0  0  0  0  0  0  0  0999 V2000
    0.3505   -0.9592    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0767    0.0891    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3039    1.1031    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.9483    2.1787    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.2966    2.8679    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.8480    2.1863    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7908    0.5104    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4189    1.8453    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8480    1.2278    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8480    2.2908    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8122   -1.1980    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9444   -2.3777    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.1367   -2.8679    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.9370   -1.9912    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  1  0      
  4  6  2  0      
  3  7  1  6      
  3  8  1  1      
  7  9  1  0      
  9 10  2  0      
  1 11  2  0      
 11 12  1  0      
 12 13  1  0      
 13 14  2  0      
 14  1  1  0
M  END

HISTEMPLATE

  my($ALASn1TemplateString)=<<ALATEMPLATE;
ALA.mol
  ChemDraw03111113532D

  9  8  0  0  0  0  0  0  0  0999 V2000
    0.7492   -0.8126    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.2014   -0.3568    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6284    0.3559    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.1966    0.8126    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.2247    0.3609    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5241   -0.7497    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2776    0.1350    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2247   -0.2743    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2247    0.4302    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  3  5  2  0      
  2  6  1  6      
  2  7  1  1      
  6  8  1  0      
  8  9  2  0  
M  END

ALATEMPLATE

  my($GLNSn1TemplateString)=<<GLNTEMPLATE;
GLN.mol
  ChemDraw03111113532D

 13 12  0  0  0  0  0  0  0  0999 V2000
    0.1872   -1.9153    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.5962   -1.3471    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.1924   -0.6209    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6736    0.0738    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.1616    0.7459    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5885    1.4586    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.1567    1.9153    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.1848    1.4636    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.2645   -1.3552    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5639    0.3530    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3174    1.2377    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2645    0.8284    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2645    1.5329    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  1  0      
  5  6  1  0      
  6  7  1  0      
  6  8  2  0      
  2  9  2  0      
  5 10  1  6      
  5 11  1  1      
 10 12  1  0      
 12 13  2  0 
M  END

GLNTEMPLATE


  my($GLUSn1TemplateString)=<<GLUTEMPLATE;
GLU.mol
  ChemDraw03111113532D

 13 12  0  0  0  0  0  0  0  0999 V2000
    0.1872   -1.9153    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.5962   -1.3471    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.1924   -0.6209    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6736    0.0738    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.1616    0.7459    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5885    1.4586    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.1567    1.9153    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.1848    1.4636    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.2645   -1.3552    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5639    0.3530    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3174    1.2377    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2645    0.8284    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2645    1.5329    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  1  0      
  5  6  1  0      
  6  7  1  0      
  6  8  2  0      
  2  9  2  0      
  5 10  1  6      
  5 11  1  1      
 10 12  1  0      
 12 13  2  0
M  END

GLUTEMPLATE

  my($PROSn1TemplateString)=<<PROTEMPLATE;
PRO.mol
  ChemDraw03111113532D

 11 11  0  0  0  0  0  0  0  0999 V2000
    1.0767   -0.8235    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3039    0.1905    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.9483    1.2661    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.2966    1.9552    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.8480    1.2737    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7908   -0.4023    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4189    0.9327    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8480    0.3152    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8480    1.3782    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.5529   -1.9552    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8122   -1.6200    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  3  5  2  0      
  2  6  1  0      
  2  7  1  1      
  6  8  1  0      
  8  9  2  0      
 10  1  1  0      
 10 11  1  0      
 11  6  1  0      
M  END

PROTEMPLATE

  my($PHESn1TemplateString)=<<PHETEMPLATE;
PHE.mol
  ChemDraw03111113532D

 15 15  0  0  0  0  0  0  0  0999 V2000
    0.2527   -0.6706    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7762    0.0851    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.2191    0.8162    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6836    1.5915    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.2138    2.0884    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.3322    1.5970    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5701    0.3889    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3020    1.3512    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3322    0.9060    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3322    1.6724    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5659   -0.5676    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0643   -1.2251    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7441   -1.9854    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0745   -2.0884    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5728   -1.4309    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  1  0      
  4  6  2  0      
  3  7  1  6      
  3  8  1  1      
  7  9  1  0      
  9 10  2  0      
  1 11  2  0      
 11 12  1  0      
 12 13  2  0      
 13 14  1  0      
 14 15  2  0      
 15  1  1  0
M  END

PHETEMPLATE

  my($TRPSn1TemplateString)=<<TRPTEMPLATE;
TRP.mol
  ChemDraw03111113532D

 18 19  0  0  0  0  0  0  0  0999 V2000
    1.3624   -0.7661    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6644    0.2917    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.1288    0.9945    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5754    1.7398    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.1238    2.2174    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.1989    1.7451    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6298    0.5837    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3720    1.5088    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3624    1.0808    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3624    1.8175    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.6792   -0.3037    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0283   -0.8106    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3092   -1.5863    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.1337   -1.5588    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7840   -0.6661    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3153   -1.2972    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0344   -2.0729    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2221   -2.2174    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  2  3  1  0      
  3  4  1  0      
  4  5  1  0      
  4  6  2  0      
  3  7  1  6      
  3  8  1  1      
  7  9  1  0      
  9 10  2  0      
  1 11  2  0      
 11 12  1  0      
 12 13  2  0      
 13 14  1  0      
 14  1  1  0      
 12 15  1  0      
 15 16  2  0      
 16 17  1  0      
 17 18  2  0      
 18 13  1  0      
 11  2  1  0 
M  END

TRPTEMPLATE

  my($TYRSn1TemplateString)=<<TYRTEMPLATE;
TYR.mol
  ChemDraw03111113532D

 16 16  0  0  0  0  0  0  0  0999 V2000
    0.2527   -0.4780    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7762    0.2778    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.2191    1.0088    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6836    1.7842    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.2138    2.2810    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.3322    1.7897    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5701    0.5815    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3020    1.5439    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3322    1.0987    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3322    1.8650    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5659   -0.3750    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0643   -1.0325    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7441   -1.7928    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0745   -1.8957    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5728   -1.2383    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1750   -2.2810    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  1  0      
  4  6  2  0      
  3  7  1  6      
  3  8  1  1      
  7  9  1  0      
  9 10  2  0      
  1 11  2  0      
 11 12  1  0      
 12 13  2  0      
 13 14  1  0      
 14 15  2  0      
 15  1  1  0      
 13 16  1  0      
M  END

TYRTEMPLATE


  my($SERSn1TemplateString)=<<SERTEMPLATE;
SER.mol
  ChemDraw03111113532D

 10  9  0  0  0  0  0  0  0  0999 V2000
    0.2323   -1.2681    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.7135   -0.5734    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.2014    0.0986    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6284    0.8113    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.1966    1.2681    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.2247    0.8165    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5241   -0.2942    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2776    0.5905    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2247    0.1812    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2247    0.8857    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  1  0      
  4  6  2  0      
  3  7  1  6      
  3  8  1  1      
  7  9  1  0      
  9 10  2  0      
M  END

SERTEMPLATE


  my($THRSn1TemplateString)=<<THRTEMPLATE;
THR.mol
  ChemDraw03111113532D

 11 10  0  0  0  0  0  0  0  0999 V2000
    0.1564   -1.2681    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6377   -0.5734    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.1255    0.0986    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5525    0.8113    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.1207    1.2681    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.1488    0.8165    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5999   -0.2942    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3535    0.5905    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3006    0.1812    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3006    0.8857    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.3006   -0.5774    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  1  0      
  4  6  2  0      
  3  7  1  6      
  3  8  1  1      
  7  9  1  0      
  9 10  2  0      
  2 11  1  6 
M  END

THRTEMPLATE

  my($METSn1TemplateString)=<<METTEMPLATE;
MET.mol
  ChemDraw03111113532D

 12 11  0  0  0  0  0  0  0  0999 V2000
    0.2672   -0.7132    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.8001   -0.1010    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.2317    0.5890    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7230    1.4090    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.2261    1.9344    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.4090    1.4148    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6029    0.1370    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3194    1.1548    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4090    0.6840    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4090    1.4944    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.8562   -1.3941    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0
    0.2328   -1.9344    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  1  0      
  4  6  2  0      
  3  7  1  6      
  3  8  1  1      
  7  9  1  0      
  9 10  2  0      
  1 11  1  0      
 11 12  1  0   
M  END

METTEMPLATE


  my($ILESn1TemplateString)=<<ILETEMPLATE;
ILE.mol
  ChemDraw03111113532D

 12 11  0  0  0  0  0  0  0  0999 V2000
    0.2043   -1.1915    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7679   -0.3779    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.1681    0.4092    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6682    1.2439    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.1625    1.7788    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.3665    1.2498    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6815   -0.0509    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3928    0.9852    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5019    0.5059    0.0000 A   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5019    1.3309    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.7827   -1.7788    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5019   -0.3826    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  1  0      
  4  6  2  0      
  3  7  1  6      
  3  8  1  1      
  7  9  1  0      
  9 10  2  0      
  1 11  1  0      
  2 12  1  1 
M  END

ILETEMPLATE


  my($VALSn1TemplateString)=<<VALTEMPLATE;
VAL.mol
  ChemDraw03111113532D

 11 10  0  0  0  0  0  0  0  0999 V2000
    0.3822   -1.2410    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5563   -0.6005    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0442    0.0715    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4712    0.7843    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0394    1.2410    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.0675    0.7893    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6813   -0.3213    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4347    0.5634    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3818    0.1541    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3818    0.8586    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.3818   -0.6045    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  1  0      
  4  6  2  0      
  3  7  1  6      
  3  8  1  1      
  7  9  1  0      
  9 10  2  0      
  2 11  1  0      
M  END

VALTEMPLATE


  my($GABASn1TemplateString)=<<GABATEMPLATE;
GABA.mol
  ChemDraw03111113532D

  9  8  0  0  0  0  0  0  0  0999 V2000
   -3.2151   -0.3506    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1434    0.2681    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0717   -0.3506    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0000    0.2681    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0717   -0.3506    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.1434    0.2681    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.2151   -0.3506    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.1474    1.2581    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1945   -1.2581    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  1  0      
  5  6  1  0      
  6  7  1  0      
  6  8  2  0      
  1  9  2  0      
M  END

GABATEMPLATE


  my($TAURSn1TemplateString)=<<TAURTEMPLATE;
TAUR.mol
  ChemDraw03111113532D

  9  8  0  0  0  0  0  0  0  0999 V2000
   -2.4873    0.0864    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4156    0.7051    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3439    0.0864    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7278    0.7051    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.7995    0.0864    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4667   -0.8211    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.1319   -0.6473    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.3292    0.8211    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.4873   -0.5795    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  1  0      
  1  6  2  0      
  5  7  2  0      
  5  8  2  0      
  5  9  1  0   
M  END

TAURTEMPLATE


  my($ARGSn1TemplateString)=<<ARGTEMPLATE;
ARG.mol
  ChemDraw03141111042D

 14 13  0  0  0  0  0  0  0  0999 V2000
    0.2436   -0.2538    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7483    0.4747    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.2112    1.1795    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6591    1.9270    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.2061    2.4060    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.2844    1.9323    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5496    0.7675    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2911    1.6953    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2844    1.2661    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2844    2.0049    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.0950   -1.5804    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6512   -0.9711    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5915   -1.1685    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.1098   -2.4060    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  1  0      
  4  6  2  0      
  3  7  1  6      
  3  8  1  1      
  7  9  1  0      
  9 10  2  0      
 11 12  1  0      
 12  1  1  0      
 11 13  1  0      
 11 14  2  0      
M  END

ARGTEMPLATE

  my($WESn1Sn2TemplateString)=<<WETEMPLATE;
WE.mol
  ChemDraw02011115312D

  5  4  0  0  0  0  0  0  0  0999 V2000
   -0.3300    1.2066    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3300    0.5466    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3300    0.0309    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.3300   -0.7322    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3094   -1.2066    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0      
  2  3  1  0      
  3  4  1  0      
  4  5  1  0      
M  END

WETEMPLATE

  my($ACSn1TemplateString)=<<ACTEMPLATE;
Ac.mol
  ChemDraw02011115312D

  5  4  0  0  0  0  0  0  0  0999 V2000
   -1.3718   -0.0743    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4573   -0.6023    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.4573   -0.0743    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.3718   -0.6023    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4362    0.6023    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  3  5  2  0      
M  END

ACTEMPLATE

  my($MYCSn1Sn2TemplateString)=<<MYCTEMPLATE;
Mycolic_core.mol
  ChemDraw02241113572D

  8  7  0  0  0  0  0  0  0  0999 V2000
    0.9413    0.1660    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.3206    0.5223    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2984    0.1632    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9189    0.5194    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3189    1.2625    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9413   -1.2625    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2973   -0.8887    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9214    1.2396    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  2  5  2  0      
  6  7  1  0      
  3  7  1  0      
  4  8  1  0      
M  END

MYCTEMPLATE

  my($ASn1TemplateString)=<<STERATEMPLATE;
Cholesterol core
  ChemDraw04081107022D

 35 38  0  0  0  0  0  0  0  0999 V2000
   -0.5145   -0.9385    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5191    0.7325    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2425    0.3108    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2352   -0.5244    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.2120    0.3190    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.2091   -0.5168    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.6560   -0.5219    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.6589    0.3140    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.9369    0.7370    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5169   -0.3458    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.2061    0.9409    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.9280    1.5497    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.6478    1.9584    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.9381    2.1838    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.3528    1.5479    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0607    1.9534    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.7658    1.5430    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.4738    1.9485    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.2082   -1.1228    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.7658    0.7180    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9542   -0.9385    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9542   -1.7697    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2352   -2.1838    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5145   -1.7697    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6686   -0.5260    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.3831   -0.9385    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.3831   -1.7635    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6686   -2.1759    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9748   -0.3611    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4435    0.9289    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.2911    1.9503    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2558   -1.1018    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -3.9513   -2.1701    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.4798   -1.7370    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.4798   -1.1457    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  4  1  0      
  1  6  1  0      
  5  2  1  0      
  2  3  1  0      
  3  4  1  0      
  5  6  1  0      
  7  8  1  0      
  8  9  1  0      
  9  5  1  0      
  6  7  1  0      
  1 10  1  1      
  5 11  1  1      
  9 12  1  0      
 12 13  1  0      
 12 14  1  1      
 13 15  1  0      
 15 16  1  0      
 16 17  1  0      
 17 18  1  0      
  6 19  1  6      
 17 20  1  0      
  4 21  1  0      
 21 22  1  0      
 22 23  2  0      
 23 24  1  0      
 24  1  1  0      
 21 25  1  0      
 25 26  1  0      
 26 27  1  0      
 27 28  1  0      
 28 22  1  0      
 21 29  1  1      
  9 30  1  6      
 12 31  1  6      
  4 32  1  6      
 27 33  1  1      
 33 34  1  0      
 34 35  2  0      
M  END

STERATEMPLATE

  my($BSn1TemplateString)=<<STERBTEMPLATE;
Cholesterol core reversed
  ChemDraw04081107022D

 35 38  0  0  0  0  0  0  0  0999 V2000
   -0.5527   -0.3850    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5574    1.2860    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2807    0.8642    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2735    0.0291    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.1737    0.8725    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.1709    0.0367    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.6177    0.0316    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.6206    0.8674    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.8986    1.2905    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5552    0.2077    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.1678    1.4944    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.8897    2.1031    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.6095    2.5118    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.8998    2.7372    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.3145    2.1014    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0224    2.5068    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.7275    2.0964    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.4355    2.5019    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.1700   -0.5693    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.7275    1.2714    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9924   -0.3850    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9924   -1.2162    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2735   -1.6303    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5527   -1.2162    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.7069    0.0275    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.4213   -0.3850    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.4213   -1.2100    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.7069   -1.6224    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0130    0.1923    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4052    1.4824    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.2528    2.5037    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2941   -0.5483    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -3.9895   -1.6166    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.4871   -2.0498    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.4871   -2.7372    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  4  1  0      
  1  6  1  0      
  5  2  1  0      
  2  3  1  0      
  3  4  1  0      
  5  6  1  0      
  7  8  1  0      
  8  9  1  0      
  9  5  1  0      
  6  7  1  0      
  1 10  1  1      
  5 11  1  1      
  9 12  1  0      
 12 13  1  0      
 12 14  1  1      
 13 15  1  0      
 15 16  1  0      
 16 17  1  0      
 17 18  1  0      
  6 19  1  6      
 17 20  1  0      
  4 21  1  0      
 21 22  1  0      
 22 23  2  0      
 23 24  1  0      
 24  1  1  0      
 21 25  1  0      
 25 26  1  0      
 26 27  1  0      
 27 28  1  0      
 28 22  1  0      
 21 29  1  1      
  9 30  1  6      
 12 31  1  6      
  4 32  1  6      
 27 33  1  1      
 33 34  1  0      
 34 35  2  0      
M  END

STERBTEMPLATE

  my($CARSn1TemplateString)=<<CARTEMPLATE;
Carnitine core
  ChemDraw04081107022D

 14 13  0  0  0  0  0  0  0  0999 V2000
   -1.5339    2.0103    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5338    1.2285    0.0000 C  0  0  0  0  0  0  0  0  0  0  0  0
   -0.8567    0.8377    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3168    0.0865    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3223    0.7255    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0804   -0.3543    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4469   -0.3543    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2106    0.0865    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2106    0.6806    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.9505   -0.3407    0.0000 O   0  5  0  0  0  0  0  0  0  0  0  0
   -1.0804   -1.1853    0.0000 N   0  3  0  0  0  0  0  0  0  0  0  0
   -1.9505   -1.4183    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3987   -1.5788    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0804   -2.0103    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  2  1  2  0      
  3  2  1  0      
  4  5  1  1      
  6  4  1  0      
  7  4  1  0      
  8  7  1  0      
  8  9  2  0      
  8 10  1  0      
  4  3  1  6      
  6 11  1  0      
 11 12  1  0      
 11 13  1  0      
 11 14  1  0      
M  CHG  2  10  -1  11   1
M  END

CARTEMPLATE

  my($MESn1TemplateString)=<<METEMPLATE;
Ac.mol
  ChemDraw02011115312D

  4  3  0  0  0  0  0  0  0  0999 V2000
   -0.7625   -0.1176    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0423   -0.6128    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.7686   -0.1506    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7686    0.6128    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  1  4  2  0      
M  END

METEMPLATE

  my($OPFBSn1TemplateString)=<<OPFBTEMPLATE;
Ac.mol
  ChemDraw02011115312D

 15 15  0  0  0  0  0  0  0  0999 V2000
   -2.1325    0.8233    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3278    0.3281    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6015    0.7903    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1386    1.5537    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.1083    0.3116    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.1248   -0.5468    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.8764   -0.9617    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.6115   -0.5182    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.6115    0.3402    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.8764    0.7550    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3374   -0.8934    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    0.8841   -1.5537    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    2.1221   -0.8769    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    2.1386    0.6417    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    0.8764    1.3350    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  1  4  2  0      
  3  5  1  0      
  5  6  2  0      
  6  7  1  0      
  7  8  2  0      
  8  9  1  0      
  9 10  2  0      
 10  5  1  0      
  6 11  1  0      
  7 12  1  0      
  8 13  1  0      
  9 14  1  0      
 10 15  1  0      
M  END

OPFBTEMPLATE

  my($MGDGSn1Sn2TemplateString)=<<MGDGTEMPLATE;
Ac.mol
  ChemDraw02011115312D

 22 22  0  0  0  0  0  0  0  0999 V2000
   -0.9221    0.6834    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6365    0.2721    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.3508    0.6834    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0652    0.2721    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.7795    0.6834    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.7795    1.5092    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2236   -0.4423    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0494   -0.4423    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.7637   -0.8551    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.7637   -1.6809    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.2737    0.9731    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.7177    0.0101    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.6419    0.2919    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.5726   -0.0137    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.1286    0.9494    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.2046    0.6676    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4985    0.7807    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.2692    0.1579    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.7795   -0.2521    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2364    0.2876    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.7353    1.1985    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.2737    1.6809    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  1  0      
  5  6  2  0      
  8  9  1  0      
  9 10  2  0      
 14 15  1  0      
 15 16  1  0      
 16 11  1  0      
 11 12  1  1      
 14 13  1  1      
 12 13  1  1      
 15 17  1  0      
 12 18  1  0      
 18 19  1  0      
  1 20  1  0      
 16 21  1  0      
 11 22  1  0      
 20 14  1  0      
  2  8  1  6      
  2  7  1  1      
M  END

MGDGTEMPLATE

  my($DGDGSn1Sn2TemplateString)=<<DGDGTEMPLATE;
Ac.mol
  ChemDraw02011115312D

 33 34  0  0  0  0  0  0  0  0999 V2000
   -2.7596    0.5322    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.4740    0.1209    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.1883    0.5322    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.9027    0.1209    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.6170    0.5322    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.6170    1.3580    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0611   -0.5935    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -3.8869   -0.5935    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.6012   -1.0062    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.6012   -1.8320    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.4362    0.8219    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.8801   -0.1411    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1956    0.1407    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2649   -0.1649    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7090    0.7983    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3670    0.5164    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3390    0.6295    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.4317    0.0067    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.9419   -0.4033    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0739    0.1364    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.8977    1.0473    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.4362    1.5297    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.1113    1.1243    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.5552    0.1613    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.4795    0.4431    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.4102    0.1375    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.9662    1.1006    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.0421    0.8188    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.3361    0.9318    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.1068    0.3091    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.6170   -0.1009    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.5729    1.3497    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.1113    1.8320    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  1  0      
  5  6  2  0      
  8  9  1  0      
  9 10  2  0      
 14 15  1  0      
 15 16  1  0      
 16 11  1  0      
 11 12  1  1      
 14 13  1  1      
 12 13  1  1      
 15 17  1  0      
 12 18  1  0      
 18 19  1  0      
  1 20  1  0      
 16 21  1  0      
 11 22  1  0      
 20 14  1  0      
  2  8  1  6      
  2  7  1  1      
 26 27  1  0      
 27 28  1  0      
 28 23  1  0      
 23 24  1  1      
 26 25  1  1      
 24 25  1  1      
 27 29  1  0      
 24 30  1  0      
 30 31  1  0      
 28 32  1  0      
 23 33  1  0      
 26 19  1  0           
M  END

DGDGTEMPLATE

  my($SQDGSn1Sn2TemplateString)=<<SQDGTEMPLATE;
Ac.mol
  ChemDraw02011115312D

 25 25  0  0  0  0  0  0  0  0999 V2000
   -1.2914    0.7692    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0058    0.3580    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.7201    0.7692    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.4345    0.3580    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.1488    0.7692    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.1488    1.5950    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5929   -0.3564    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4187   -0.3564    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1330   -0.7692    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1330   -1.5950    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.9044    1.0589    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.3483    0.0959    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2726    0.3777    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.2033    0.0721    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7592    1.0353    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.8352    0.7534    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.1292    0.8665    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.8999    0.2437    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.4101   -0.1663    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6057    0.3734    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.3659    1.2843    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.5324    0.6735    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.0556   -0.7414    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.0093   -0.5553    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.1488    0.0494    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  1  0      
  5  6  2  0      
  8  9  1  0      
  9 10  2  0      
 14 15  1  0      
 15 16  1  0      
 16 11  1  0      
 11 12  1  1      
 14 13  1  1      
 12 13  1  1      
 15 17  1  0      
 12 18  1  0      
 18 19  1  0      
  1 20  1  0      
 16 21  1  0      
 11 22  1  0      
 20 14  1  0      
  2  8  1  6      
  2  7  1  1      
 19 23  2  0      
 19 24  1  0      
 19 25  2  0           
M  END

SQDGTEMPLATE


  my($DFPU3DAGSn1Sn2TemplateString)=<<DFPU3DAGTEMPLATE;
Ac.mol
  ChemDraw02011115312D

 22 22  0  0  0  0  0  0  0  0999 V2000
   -0.9961    0.5327    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7110    0.1211    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4258    0.5327    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1408    0.1211    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.8556    0.5327    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.8556    1.3591    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2978   -0.5938    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1242   -0.5938    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.8391   -1.0069    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.8391   -1.8333    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3098    0.1365    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.4113    0.5171    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.1125    0.1046    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    1.7725    0.4758    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.5563    0.0427    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.3232    0.5049    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.3232    1.4002    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.5563    1.8333    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.7725    1.3711    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.5563   -0.5348    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    3.8556    1.7546    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    0.4113    1.0740    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  1  0      
  5  6  2  0      
  8  9  1  0      
  9 10  2  0      
  1 11  1  0      
  2  8  1  6      
  2  7  1  1      
 11 12  1  0      
 12 13  1  0      
 13 14  1  0      
 14 15  2  0      
 15 16  1  0      
 16 17  2  0      
 17 18  1  0      
 18 19  2  0      
 19 14  1  0      
 15 20  1  0      
 17 21  1  0      
 12 22  2  0                
M  END

DFPU3DAGTEMPLATE


  my($DFPU2DAGSn1Sn2TemplateString)=<<DFPU2DAGTEMPLATE;
Ac.mol
  ChemDraw02011115312D

 23 23  0  0  0  0  0  0  0  0999 V2000
    1.7912    0.4501    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0762    0.0386    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3614    0.4501    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3536    0.0386    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0684    0.4501    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0684    1.2765    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.4894   -0.6764    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.7661   -0.7382    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.2782   -1.1514    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.2719   -1.7306    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3356   -0.7799    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9670   -1.1073    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7360   -0.6847    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4866   -1.1282    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4746   -2.0055    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7120   -2.4282    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9549   -1.9847    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4762   -2.2969    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0541   -2.3323    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    1.8029    1.1082    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.4904    1.6857    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.4766    2.4282    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0541    1.3832    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  1  0      
  5  6  2  0      
  8  9  1  0      
  9 10  2  0      
  2  8  1  6      
  2  7  1  1      
  9 11  1  0      
 11 12  1  0      
 12 13  2  0      
 13 14  1  0      
 14 15  2  0      
 15 16  1  0      
 16 17  2  0      
 17 12  1  0      
 17 18  1  0      
 15 19  1  0      
  1 20  1  0      
 20 21  1  0      
 21 22  1  0      
 21 23  2  0                     
M  END

DFPU2DAGTEMPLATE

  my($DFPU23MAGSn1TemplateString)=<<DFPU23MAGTEMPLATE;
Ac.mol
  ChemDraw02011115312D

 31 32  0  0  0  0  0  0  0  0999 V2000
   -0.0032    0.7889    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7181    0.3773    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4329    0.7889    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1480    0.3773    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.8628    0.7889    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.8628    1.6153    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3050   -0.3376    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0283   -0.3995    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5162   -0.8126    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5225   -1.3918    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.6830    0.3927    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.4041    0.7732    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.1054    0.3607    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    2.7654    0.7320    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.5491    0.2989    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.3161    0.7611    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.3161    1.6564    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.5491    2.0895    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.7654    1.6273    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.5491   -0.2786    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    4.8485    2.0107    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    1.4041    1.3301    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1300   -0.4412    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -2.7614   -0.7686    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.5304   -0.3460    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.2810   -0.7894    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.2690   -1.6668    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.5064   -2.0895    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.7493   -1.6459    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2706   -1.9582    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
   -4.8485   -1.9936    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  1  0      
  5  6  2  0      
  8  9  1  0      
  9 10  2  0      
  1 11  1  0      
  2  8  1  6      
  2  7  1  1      
 11 12  1  0      
 12 13  1  0      
 13 14  1  0      
 14 15  2  0      
 15 16  1  0      
 16 17  2  0      
 17 18  1  0      
 18 19  2  0      
 19 14  1  0      
 15 20  1  0      
 17 21  1  0      
 12 22  2  0      
  9 23  1  0      
 23 24  1  0      
 24 25  2  0      
 25 26  1  0      
 26 27  2  0      
 27 28  1  0      
 28 29  2  0      
 29 24  1  0      
 29 30  1  0      
 27 31  1  0                     
M  END

DFPU23MAGTEMPLATE

# original: "WESn1Sn2" => "WE|wax ester|0.1632|0.5866|-1.2066|-0.7622|0|0|0|0|2|5|0|0|1|2|0|0|0|0|0|LM|||$WESn1Sn2TemplateString",
# Format: ID => AbbrevID|HeadGroupName|Y1Sn1|Y2Sn1|Y1Sn2|Y2Sn2|Y1Sn3|Y2Sn3|Y1Sn4|Y2Sn4|Sn1ChainAtomNum|Sn2ChainAtomNum|Sn3ChainAtomNum|Sn4ChainAtomNum|Sn1ChainCarbons|Sn2ChainCarbons|Sn3ChainCarbons|Sn4ChainCarbons|Sn2CAtomNum|Sn2OAtomNum|Sn2HAtomNum|LMCategory|LMMainClass|LMSubClass|TemplateString
  # or
  # Format: ID => AbbrevID|HeadGroupName|Y1Sn1|Y2Sn1|Y1Sn2|Y2Sn2|Y1Sn3|Y2Sn3|Y1Sn4|Y2Sn4|Sn1ChainAtomNum|Sn2ChainAtomNum|Sn3ChainAtomNum|Sn4ChainAtomNum|Sn1ChainCarbons|Sn2ChainCarbons|Sn3ChainCarbons|Sn4ChainCarbons|Sn2CAtomNum|Sn2OAtomNum|Sn2HAtomNum|LMCategory|LMMainClass|LMSubClass|TemplateString|Sn1ReverseChainGrowth|Sn2ReverseChainGrowth|Sn3ReverseChainGrowth|Sn4ReverseChainGrowth
  #
  # Note: Template structure for Both Alkyl and Alkenyl substituted GP is the same; Different entries are used to facilitate assigning LM identity
  %LMTemplatesDataMap = (
"DFPU23MAGSn1" => "DFPU23MAG|2,3-diDFPU-monoacylglycerol|0.3773|0.7889|0|0|0|0|0|0|5|0|0|0|1|0|0|0|0|0|0|LM|||$DFPU23MAGSn1TemplateString",
"DFPU2DAGSn1Sn2" => "DFPU2DAG|2-DFPU-diacylglycerol|0.0386|0.4501|2.4282|2.9238|0|0|0|0|5|22|0|0|1|2|0|0|0|0|0|LM|||$DFPU2DAGSn1Sn2TemplateString",
"DFPU3DAGSn1Sn2" => "DFPU3DAG|3-DFPU-diacylglycerol|0.1211|0.5327|-1.0069|-0.5938|0|0|0|0|5|9|0|0|1|1|0|0|0|0|0|LM|||$DFPU3DAGSn1Sn2TemplateString",
"SQDGSn1Sn2" => "SQDG|sulfoquinovosyldiacylglycerol|0.3580|0.7692|-0.7692|-0.3564|0|0|0|0|5|9|0|0|1|1|0|0|0|0|0|LM|||$SQDGSn1Sn2TemplateString",
"DGDGSn1Sn2" => "DGDG|digalactosyldiacylglycerol|0.1209|0.5322|-1.0062|-0.5935|0|0|0|0|5|9|0|0|1|1|0|0|0|0|0|LM|||$DGDGSn1Sn2TemplateString",
"MGDGSn1Sn2" => "MGDG|monogalactosyldiacylglycerol|0.2721|0.6834|-0.8551|-0.4423|0|0|0|0|5|9|0|0|1|1|0|0|0|0|0|LM|||$MGDGSn1Sn2TemplateString",
               "OPFBSn1" => "OPFB|Pentafluorobenzyl ester|0.3281|0.8233|0|0|0|0|0|0|1|0|0|0|1|0|0|0|0|0|0|LM|||$OPFBSn1TemplateString",
               "MESn1" => "ME|methyl ester|-0.6128|-0.1176|0|0|0|0|0|0|1|0|0|0|1|0|0|0|0|0|0|LM|||$MESn1TemplateString",
               "CARSn1" => "CAR|carnitine|0.8377|1.2285|0|0|0|0|0|0|2|0|0|0|1|0|0|0|0|0|0|LM|||$CARSn1TemplateString",
      "CESn1" => "CE|cholest-5-en-3beta-ol|-2.1701|-1.7370|0|0|0|0|0|0|34|0|0|0|1|0|0|0|0|0|0|LM|||$ASn1TemplateString",
"CErevSn1" => "CErev|cholest-5-en-3beta-ol|-2.0498|-1.6166|0|0|0|0|0|0|34|0|0|0|1|0|0|0|0|0|0|LM|||$BSn1TemplateString",

     "ACSn1" => "AC|acetate|-0.6023|-0.0743|0|0|0|0|0|0|1|0|0|0|1|0|0|0|0|0|0|LM|||$ACSn1TemplateString",
  
     "ARGSn1" => "ARG|arginine|0.7675|1.2661|0|0|0|0|0|0|9|0|0|0|1|0|0|0|0|0|0|LM|||$ARGSn1TemplateString",
   "ASXSn1" => "ASX|asparagine|-0.0656|0.6518|0|0|0|0|0|0|9|0|0|0|1|0|0|0|0|0|0|LM|||$ASNSn1TemplateString",
     "HISSn1" => "HIS|histidine|0.5104|1.2278|0|0|0|0|0|0|9|0|0|0|1|0|0|0|0|0|0|LM|||$HISSn1TemplateString",
     "ALASn1" => "ALA|alanine|-0.7497|-0.2743|0|0|0|0|0|0|8|0|0|0|1|0|0|0|0|0|0|LM|||$ALASn1TemplateString",
    "GLNSn1" => "GLN|glutamine|0.3530|0.8284|0|0|0|0|0|0|12|0|0|0|1|0|0|0|0|0|0|LM|||$GLNSn1TemplateString",
    "GLUSn1" => "GLU|glutamate|0.3530|0.8284|0|0|0|0|0|0|12|0|0|0|1|0|0|0|0|0|0|LM|||$GLUSn1TemplateString",
      "PROSn1" => "PRO|proline|-0.4023|0.3152|0|0|0|0|0|0|8|0|0|0|1|0|0|0|0|0|0|LM|||$PROSn1TemplateString",
 "PHESn1" => "PHE|phenylalanine|0.3889|0.9060|0|0|0|0|0|0|9|0|0|0|1|0|0|0|0|0|0|LM|||$PHESn1TemplateString",
 "TRPSn1" => "TRP|phenylalanine|0.5837|1.0808|0|0|0|0|0|0|9|0|0|0|1|0|0|0|0|0|0|LM|||$TRPSn1TemplateString",
      "TYRSn1" => "TYR|tyrosine|0.5815|1.0987|0|0|0|0|0|0|9|0|0|0|1|0|0|0|0|0|0|LM|||$TYRSn1TemplateString",
       "SERSn1" => "SER|serine|-0.2942|0.1812|0|0|0|0|0|0|9|0|0|0|1|0|0|0|0|0|0|LM|||$SERSn1TemplateString",
    "THRSn1" => "THR|threonine|-0.2942|0.1812|0|0|0|0|0|0|9|0|0|0|1|0|0|0|0|0|0|LM|||$THRSn1TemplateString",
    "METSn1" => "MET|methionine|0.1370|0.6840|0|0|0|0|0|0|9|0|0|0|1|0|0|0|0|0|0|LM|||$METSn1TemplateString",
   "ILESn1" => "ILE|isoleucine|-0.0509|0.5059|0|0|0|0|0|0|9|0|0|0|1|0|0|0|0|0|0|LM|||$ILESn1TemplateString",
       "VALSn1" => "VAL|valine|-0.3213|0.1541|0|0|0|0|0|0|9|0|0|0|1|0|0|0|0|0|0|LM|||$VALSn1TemplateString",
       "GABASn1" => "GABA|gamma-aminobutyric acid|-0.3506|0.2681|0|0|0|0|0|0|1|0|0|0|1|0|0|0|0|0|0|LM|||$GABASn1TemplateString",
       "TAURSn1" => "TAUR|taurine|0.0864|0.7051|0|0|0|0|0|0|1|0|0|0|1|0|0|0|0|0|0|LM|||$TAURSn1TemplateString",

"WESn1Sn2" => "WE|wax ester|-1.2066|-0.7622|0.1632|0.5866|0|0|0|0|5|2|0|0|2|1|0|0|0|0|0|LM|||$WESn1Sn2TemplateString",
"MYCSn1Sn2" => "MYC|Mycolic acid|0.1632|0.5194|-1.2625|-0.8887|0|0|0|0|4|6|0|0|0|2|0|0|0|0|0|LM|||$MYCSn1Sn2TemplateString",
"TEST2Sn1Sn2Sn3Sn4" => "TEST2|sulfolipid-1|-1.1253|-0.6383|-2.6257|-2.1233|1.1253|1.5479|2.6257|3.0742|15|16|14|13|1|1|1|1|0|0|0|LM|||$TESTTemplateString|0|0|1|1",
"PAT18Sn1Sn2Sn3Sn4" => "PAT18|polyacyltrehalose (R1=C18)|-3.5314|-3.0607|-0.3210|0.1655|0.2074|0.6666|2.9488|3.4491|27|49|47|45|1|1|1|1|0|0|0|LM|||$PAT18ReverseSn2TemplateString|1|1|0|0",
"PAT16Sn1Sn2Sn3Sn4" => "PAT16|polyacyltrehalose (R1=C16)|-3.5060|-3.0607|-0.3210|0.1655|0.2074|0.6666|2.9488|3.4491|27|47|45|43|1|1|1|1|0|0|0|LM|||$PAT16ReverseSn2TemplateString|1|1|0|0",
               "SL1Sn1Sn2Sn3Sn4" => "SL1|sulfolipid I|-1.7186|-1.2777|-3.3309|-2.8196|-0.6198|-0.0928|2.8622|3.2897|23|30|28|26|1|1|1|1|0|0|0|LM|||$SL1TemplateString|0|1|1|1",
               "SL2Sn1Sn2Sn3Sn4" => "SL2|sulfolipid II|-1.7186|-1.2777|-3.3309|-2.8196|-0.6198|-0.0928|2.8622|3.2897|23|30|28|26|1|1|1|1|0|0|0|LM|||$SL1TemplateString|0|1|1|1",

      "SL2pSn1Sn2Sn3Sn4" => "SL2p|sulfolipid-II'|-1.9977|-1.5568|-2.9486|-2.4330|-0.7339|-0.3719|2.5831|3.0106|22|34|27|25|1|1|1|1|0|0|0|LM|||$SL2TemplateString|0|1|1|1",
      "SL3Sn1Sn2Sn3" => "SL3|sulfolipid-III|-0.9892|-0.5483|-2.6015|-2.0902|0.1096|0.6365|0|0|23|28|26|0|1|1|1|0|0|0|0|LM|||$SL3TemplateString|0|1|1",

"Ac2PIM2Sn1Sn2Sn3Sn4" => "Ac2PIM2|glycerophosphoinositoldiacyldimannoside|1.6806|2.0907|0.5470|0.9587|-3.5874|-3.2188|-0.3286|0.0632|8|16|53|56|2|2|2|2|0|0|0|LM|||$AC2PIM2ReverseSn2TemplateString|0|0|1|1",
"Ac2PIM3Sn1Sn2Sn3Sn4" => "Ac2PIM3|glycerophosphoinositoldiacyltrimannoside|0.4428|0.8529|-0.6908|-0.2791|-4.8252|-4.4566|-1.5664|-1.1746|8|16|53|56|2|2|2|2|0|0|0|LM|||$AC2PIM3ReverseSn2TemplateString|0|0|1|1",
"Ac2PIM4Sn1Sn2Sn3Sn4" => "Ac2PIM4|glycerophosphoinositoldiacyltetramannoside|-0.7637|-0.3536|-1.8972|-1.4856|-6.0316|-5.6631|-2.7729|-2.3810|8|16|53|56|2|2|2|2|0|0|0|LM|||$AC2PIM4ReverseSn2TemplateString|0|0|1|1",
"Ac2PIM5Sn1Sn2Sn3Sn4" => "Ac2PIM5|glycerophosphoinositoldiacyltetramannoside|-1.9599|-1.5498|-3.0935|-2.6818|-7.2279|-6.8593|-3.9691|-3.5773|8|16|53|56|2|2|2|2|0|0|0|LM|||$AC2PIM5ReverseSn2TemplateString|0|0|1|1",
"Ac2PIM6Sn1Sn2Sn3Sn4" => "Ac2PIM6|glycerophosphoinositoldiacyltetramannoside|-3.0892|-2.6791|-4.2228|-3.8111|-8.5016|-8.1330|-5.0984|-4.7478|8|16|53|56|2|2|2|2|0|0|0|LM|||$AC2PIM6ReverseSn2TemplateString|0|0|1|1",

"DIMA20Sn1Sn2" => "DIMA20|phthiocerol dimycosocerosate (DIMA,n=20,R3=Et)|2.0717|2.4439|0.3140|0.7047|0|0|0|0|20|18|0|0|1|1|0|0|0|0|0|LM|||$DIMA20TemplateString",
"DIMA22Sn1Sn2" => "DIMA22|phthiocerol dimycosocerosate (DIMA,n=22,R3=Et)|2.0718|2.4926|0.3141|0.7047|0|0|0|0|20|18|0|0|1|1|0|0|0|0|0|LM|||$DIMA22TemplateString",
"DIMA20MeSn1Sn2" => "DIMA20Me|phthiocerol dimycosocerosate (DIMA,n=20,R3=Me)|1.7334|2.1057|-0.0243|0.3664|0|0|0|0|19|17|0|0|1|1|0|0|0|0|0|LM|||$DIMA20MeTemplateString",
"DIMA22MeSn1Sn2" => "DIMA22Me|phthiocerol dimycosocerosate (DIMA,n=22,R3=Me)|1.7334|2.1057|-0.0243|0.3664|0|0|0|0|19|17|0|0|1|1|0|0|0|0|0|LM|||$DIMA22MeTemplateString",

"DIMB20Sn1Sn2" => "DIMB20|phthiocerol dimycosocerosate (DIMB,n=20,R3=Et)|2.0717|2.4439|0.3140|0.7047|0|0|0|0|19|17|0|0|1|1|0|0|0|0|0|LM|||$DIMB20TemplateString",
"DIMB22Sn1Sn2" => "DIMB22|phthiocerol dimycosocerosate (DIMB,n=22,R3=Et)|2.0718|2.4926|0.3141|0.7047|0|0|0|0|19|17|0|0|1|1|0|0|0|0|0|LM|||$DIMB22TemplateString",
"DIMB20MeSn1Sn2" => "DIMB20Me|phthiocerol dimycosocerosate (DIMB,n=20,R3=Me)|1.7334|2.1057|-0.0243|0.3664|0|0|0|0|18|16|0|0|1|1|0|0|0|0|0|LM|||$DIMB20MeTemplateString",
"DIMB22MeSn1Sn2" => "DIMB22Me|phthiocerol dimycosocerosate (DIMB,n=22,R3=Me)|1.7334|2.1057|-0.0243|0.3664|0|0|0|0|18|16|0|0|1|1|0|0|0|0|0|LM|||$DIMB22MeTemplateString",

"DIMBG15Sn1Sn2" => "DIMBG15|phthiocerol dimycosocerosate (Glycosylated DIMB,n=15,R3=Et)|1.3144|1.6866|-0.4434|-0.0527|0|0|0|0|19|17|0|0|1|1|0|0|0|0|0|LM|||$DIMBG15TemplateString",
"DIMBG15MeSn1Sn2" => "DIMBG15Me|phthiocerol dimycosocerosate (Glycosylated DIMB,n=15,R3=Et)|0.9761|1.3483|-0.7818|-0.3911|0|0|0|0|18|16|0|0|1|1|0|0|0|0|0|LM|||$DIMBG15MeTemplateString",
"DIMBG17Sn1Sn2" => "DIMBG17|phthiocerol dimycosocerosate (Glycosylated DIMB,n=15,R3=Et)|1.2937|1.6660|-0.4640|-0.0733|0|0|0|0|19|17|0|0|1|1|0|0|0|0|0|LM|||$DIMBG17TemplateString",
"DIMBG17MeSn1Sn2" => "DIMBG17Me|phthiocerol dimycosocerosate (Glycosylated DIMB,n=15,R3=Et)|0.9554|1.3276|-0.8023|-0.4117|0|0|0|0|18|16|0|0|1|1|0|0|0|0|0|LM|||$DIMBG17MeTemplateString",

"DIMAG15Sn1Sn2" => "DIMAG15|phthiocerol dimycosocerosate (Glycosylated DIMA,n=15,R3=Et)|1.3144|1.6866|-0.4434|-0.0527|0|0|0|0|19|17|0|0|1|1|0|0|0|0|0|LM|||$DIMAG15TemplateString",
"DIMAG15MeSn1Sn2" => "DIMAG15Me|phthiocerol dimycosocerosate (Glycosylated DIMA,n=15,R3=Et)|0.9761|1.3483|-0.7818|-0.3911|0|0|0|0|18|16|0|0|1|1|0|0|0|0|0|LM|||$DIMAG15MeTemplateString",
"DIMAG17Sn1Sn2" => "DIMAG17|phthiocerol dimycosocerosate (Glycosylated DIMA,n=15,R3=Et)|1.2937|1.6660|-0.4640|-0.0733|0|0|0|0|19|17|0|0|1|1|0|0|0|0|0|LM|||$DIMAG17TemplateString",
"DIMAG17MeSn1Sn2" => "DIMAG17Me|phthiocerol dimycosocerosate (Glycosylated DIMA,n=15,R3=Et)|0.9554|1.3276|-0.8023|-0.4117|0|0|0|0|18|16|0|0|1|1|0|0|0|0|0|LM|||$DIMAG17MeTemplateString",

     "PIM1Sn1Sn2" => "PIM1|glycerophosphoinositolmonomannoside|1.8141|2.2242|0.6807|1.0923|0|0|0|0|8|16|0|0|2|2|0|0|0|0|0|LM|||$PIM1Sn1Sn2TemplateString",
      "PIM2Sn1Sn2" => "PIM2|glycerophosphoinositoldimannoside|1.0947|1.5047|-0.0387|0.3729|0|0|0|0|8|16|0|0|2|2|0|0|0|0|0|LM|||$PIM2Sn1Sn2TemplateString",
   "PIM3Sn1Sn2" => "PIM3|glycerophosphoinositoltrimannoside|-0.2270|0.1831|-1.3604|-0.9488|0|0|0|0|8|16|0|0|2|2|0|0|0|0|0|LM|||$PIM3Sn1Sn2TemplateString",
"PIM4Sn1Sn2" => "PIM4|glycerophosphoinositoltetramannoside|-1.0932|-0.6832|-2.2266|-1.8150|0|0|0|0|8|16|0|0|2|2|0|0|0|0|0|LM|||$PIM4Sn1Sn2TemplateString",
"PIM5Sn1Sn2" => "PIM5|glycerophosphoinositolpentamannoside|-2.2482|-1.8382|-3.3816|-2.9700|0|0|0|0|8|16|0|0|2|2|0|0|0|0|0|LM|||$PIM5Sn1Sn2TemplateString",
 "PIM6Sn1Sn2" => "PIM6|glycerophosphoinositolhexamannoside|-3.4009|-2.9908|-4.5342|-4.1226|0|0|0|0|8|16|0|0|2|2|0|0|0|0|0|LM|||$PIM6Sn1Sn2TemplateString",

"LPIM1Sn1" => "LPIM1|glycerophosphoinositolmonomannoside|1.8143|2.2244|0|0|0|0|0|0|8|0|0|0|2|0|0|0|0|0|0|LM|||$PIM1LYSOSn1TemplateString",
"LPIM2Sn1" => "LPIM2|glycerophosphoinositoldimannoside|1.0947|1.5047|0|0|0|0|0|0|8|0|0|0|2|0|0|0|0|0|0|LM|||$PIM2LYSOSn1TemplateString",
"LPIM3Sn1" => "LPIM3|glycerophosphoinositoltrimannoside|-0.2270|0.1831|0|0|0|0|0|0|8|0|0|0|2|0|0|0|0|0|0|LM|||$PIM3LYSOSn1TemplateString",
"LPIM4Sn1" => "LPIM4|glycerophosphoinositoltetramannoside|-1.0932|-0.6832|0|0|0|0|0|0|8|0|0|0|2|0|0|0|0|0|0|LM|||$PIM4LYSOSn1TemplateString",
"LPIM5Sn1" => "LPIM5|glycerophosphoinositolpentamannoside|-2.2483|-1.8383|0|0|0|0|0|0|8|0|0|0|2|0|0|0|0|0|0|LM|||$PIM5LYSOSn1TemplateString",
"LPIM6Sn1" => "LPIM6|glycerophosphoinositolhexamannoside|-3.4009|-2.9908|0|0|0|0|0|0|8|0|0|0|2|0|0|0|0|0|0|LM|||$PIM6LYSOSn1TemplateString",

"Ac1PIM1Sn1Sn2Sn3" => "Ac1PIM1|glycerophosphoinositolacylmonomannoside|2.3263|2.7364|1.1928|1.6044|-2.8588|-2.4463|0|0|8|16|42|0|2|2|2|0|0|0|0|LM|||$Ac1PIM1Sn1Sn2Sn3TemplateString",
"Ac1PIM2Sn1Sn2Sn3" => "Ac1PIM2|glycerophosphoinositolacyldimannoside|1.6288|2.0388|0.4954|0.9070|-3.6617|-3.2080|0|0|8|16|53|0|2|2|2|0|0|0|0|LM|||$Ac1PIM2Sn1Sn2Sn3TemplateString",

"Ac1PIM3Sn1Sn2Sn3" => "Ac1PIM3|glycerophosphoinositolacyltrimannoside|0.3269|0.7370|-0.8065|-0.3949|-4.9825|-4.5494|0|0|8|16|64|0|2|2|2|0|0|0|0|LM|||$Ac1PIM3Sn1Sn2Sn3TemplateString",
"Ac1PIM4Sn1Sn2Sn3" => "Ac1PIM4|glycerophosphoinositolacyltetramannoside|-0.5599|-0.1499|-1.6933|-1.2817|-5.8488|-5.3950|0|0|8|16|75|0|2|2|2|0|0|0|0|LM|||$Ac1PIM4Sn1Sn2Sn3TemplateString",

"Ac1PIM5Sn1Sn2Sn3" => "Ac1PIM5|glycerophosphoinositolacylpentamannoside|-1.6806|-1.2705|-2.8140|-2.4024|-7.0180|-6.5849|0|0|8|16|86|0|2|2|2|0|0|0|0|LM|||$Ac1PIM5Sn1Sn2Sn3TemplateString",
"Ac1PIM6Sn1Sn2Sn3" => "Ac1PIM6|glycerophosphoinositolacylhexamannoside|-2.7890|-2.3789|-3.9223|-3.5107|-8.1318|-7.6368|0|0|8|16|97|0|2|2|2|0|0|0|0|LM|||$Ac1PIM6Sn1Sn2Sn3TemplateString",

"DATSn1Sn2" => "DAT|di-O-acyltrehalose| -1.0216| -0.6463|-2.3034|-1.8434|0|0|0|0|25|27|0|0|2|1|0|0|0|0|0|LM|||$DATSn1Sn2TemplateString",

"CoASn1" => "CoA|CoA|3.5305|3.9404|0|0|0|0|0|0|52|0|0|0|2|0|0|0|0|0|0|LM|||$COASn1TemplateString",

"CLSn1Sn2Sn3Sn4" => "CL|glycerol|2.0639|2.4763|0.9379|1.3504|-1.3477|-0.9352|-2.4747|-2.0622|12|8|33|37|2|2|2|2|6|20|19|GP|12|01|$CLSn1Sn2Sn3Sn4TemplateString",

"AC2SGLSn1Sn2" => "AC2SGL|di-O-acyl-2'-sulfotrehalose|-1.1046|-0.6988|-2.4906|-1.9932|0|0|0|0|25|27|0|0|2|1|0|0|0|0|0|LM|0|0|$AC2SGLSn1Sn2TemplateString",

			);
}
# Format: ID => AbbrevID|HeadGroupName|Y1Sn1|Y2Sn1|Y1Sn2|Y2Sn2|Y1Sn3|Y2Sn3|Y1Sn4|Y2Sn4|Sn1ChainAtomNum|Sn2ChainAtomNum|Sn3ChainAtomNum|Sn4ChainAtomNum|Sn1ChainCarbons|Sn2ChainCarbons|Sn3ChainCarbons|Sn4ChainCarbons|Sn2CAtomNum|Sn2OAtomNum|Sn2HAtomNum|LMCategory|LMMainClass|LMSubClass|TemplateString
# or
# Format: ID => AbbrevID|HeadGroupName|Y1Sn1|Y2Sn1|Y1Sn2|Y2Sn2|Y1Sn3|Y2Sn3|Y1Sn4|Y2Sn4|Sn1ChainAtomNum|Sn2ChainAtomNum|Sn3ChainAtomNum|Sn4ChainAtomNum|Sn1ChainCarbons|Sn2ChainCarbons|Sn3ChainCarbons|Sn4ChainCarbons|Sn2CAtomNum|Sn2OAtomNum|Sn2HAtomNum|LMCategory|LMMainClass|LMSubClass|TemplateString|Sn1ReverseChainGrowth|Sn2ReverseChainGrowth|Sn3ReverseChainGrowth|Sn4ReverseChainGrowth

# Initialize supported head groups...
sub _InitializeSupportedHeadGroupsData {
  my($GPType, $GPHeadGroup);
  %LMSupportedHeadGroupMap = ();

  for $GPType (keys %LMTemplatesDataMap) {
    ($GPHeadGroup) = split /\|/, $LMTemplatesDataMap{$GPType};
    if (!(exists $LMSupportedHeadGroupMap{$GPHeadGroup})) {
      $LMSupportedHeadGroupMap{$GPHeadGroup} = $GPHeadGroup;
    }
  }
}


1;

__END__

=head1 NAME

LMStr - LIPID MAPS arbitrary structure generation methods

=head1 SYNOPSIS

use LMStr;

use LMStr qw(:all);

=head1 DESCRIPTION

LMStr module provides these methods:

    ExpandLMCmpdAbbrevs - Expand abbreviation
    GenerateCmpdOntologyData - Generate ontology data
    GenerateCmpdOntologySDDataLines - Generate ontology data lines for SD file
    GenerateLMChainStrData - Generate chain structure data
    GenerateSDFile - Generate SD file
    GetLMTemplatesData - Get templates data
    GetLMSupportedHeadGroupMap - Get supported headgroups data
    GetLMTemplateID - Get templates ID
    IsLMChainsAbbrevSupported - Is it a supported abbreviation
    ParseLMAbbrev - Parse abbreviation
    SetupLMCmpdAbbrevTemplateDataMap - Setup template structure data map
    ValidateLMAbbrev - Validate abbreviation

=head1 METHODS

=over 4

=item B<ExpandLMCmpdAbbrevs>

    $ExpandedAbbrevArrayRef = ExpandLMCmpdAbbrevs($CmpdAbbrev);

Return a reference to an array containing complete LM abbreviations. Wild card
characters in LM abbreviation name are expanded to generate fully qualified
LM abbreviations.

=item B<GenerateCmpdOntologyData>

    $DataHashRef = GenerateCmpdOntologyData($CmpDataRef);

Return a reference to a hash containing ontology data with hash keys and values
corresponding to property names and values.

=item B<GenerateCmpdOntologySDDataLines>

    $DataLinesArrayRef = GenerateCmpdOntologySDDataLines($CmpdDataRef);

Return a reference to an array containing ontology data lines suitable for
generate SD file data block.

=item B<GenerateLMChainStrData>

    ($AtomLinesArrayRef, $BondLinesArrayRef) =
       GenerateLMChainStrData($ChainType, $CmpdDataRef);

Return array references containing atom and bond data lines for SD file. Appropriate atom
and bond data lines are generated using chain type and abbreviation template data.

=item B<GenerateSDFile>

    GenerateSDFile($SDFileName, $CmdAbbrevsRef);

Generate a SD file for compound abbreviations. Structure data for specified abbreviation
is generated sequentially and written to SD file.

=item B<GetLMTemplatesData>

    $TemplatesDataRef = GetLMTemplatesData();

Return a reference to a hash containing LM templates data

=item B<GetLMSupportedHeadGroupMap>

    $SupportedHeadGroupDataRef = GetLMSupportedHeadGroupMap();

Return a reference to a hash containing supported head groups data.

=item B<GetLMTemplateID>

    $HeadGroupID = GetLMTemplateID($HeadGroupAbbrev, $ChainsAbbrev);

Return a supported template ID for compound abbreviation.

=item B<IsLMChainsAbbrevSupported>

    $Status = IsLMChainsAbbrevSupported($Abbrev, $PrintWarning);

Return 1 or 0 based on whether LM abbreviated is supported. For unsupported LM abbreviations,
a warning is printed unless PrintWarning flag is set.

=item B<ParseLMAbbrev>

    ($HeadGroup, $ChainsAbbrev, $AbbrevModifier) =
       ParseLMAbbrev($Abbrev);

Parse LM abbreviation and return these values: HeadGroup, ChainsAbbrev,
AbbrevModifier.

=item B<SetupLMCmpdAbbrevTemplateDataMap>

    $AbbrevTemplateDataMapRef =
       SetupLMCmpdAbbrevTemplateDataMap($Abbrev);

Return a reference to a hash containing template data for compound abbreviation. The
template data is used to generate SD file for compound abbreviation.

=item B<ValidateLMAbbrev>

    $Status = ValidateLMAbbrev($Abbrev);

Return 1 or 0 based on whether a LM abbreviation is valid.

=back

=head1 AUTHOR

Manish Sud

=head1 CONTRIBUTOR

Eoin Fahy

=head1 SEE ALSO

ChainAbbrev.pm, ChainStr.pm, LMAPSStr.pm

=head1 COPYRIGHT

Copyright (C) 2006-2012. The Regents of the University of California. All Rights Reserved.

=head1 LICENSE

Modified BSD License

=cut
