package GPStr;
#
# File: GPStr.pm
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
@EXPORT_OK = qw(GenerateCmpdOntologyData GenerateCmpdOntologySDDataLines GenerateGPChainStrData GenerateSDFile GetGPTemplatesData GetGPSupportedHeadGroupMap GetGPTemplateID IsGPChainsAbbrevSupported ParseGPAbbrev ProcessGPCmpdAbbrevs SetupGPCmpdAbbrevTemplateDataMap ValidateGPAbbrev);
%EXPORT_TAGS = (
		all  => [@EXPORT, @EXPORT_OK],
		);

# Initialize GP data...
my(%GPTemplatesDataMap, %GPSupportedHeadGroupMap);
_InitializeData();

# Process specified compound abbreviations containing any wild cards
# and count the number of valid abbreviations along with optional generation
# of a SD file.
#
sub ProcessGPCmpdAbbrevs {
  my($CmpdAbbrevsRef, $AllowArbitraryChainAbbrev, $WriteSDFile, $SDFileName) = @_;

  $AllowArbitraryChainAbbrev = defined $AllowArbitraryChainAbbrev ? $AllowArbitraryChainAbbrev : 0;

  $WriteSDFile = defined $WriteSDFile ? $WriteSDFile : 0;
  if ($WriteSDFile &&  IsEmpty($SDFileName)) {
    warn "Warning: GPStr::ProcessGPCmpdAbbrevs: No SD file name specified. Ingoring structure generation.\n";
    return;
  }

  if ($WriteSDFile) {
    print "Generating SD file $SDFileName...\n";

    open SDFILE, ">$SDFileName" or die "Error: Couldn't open $SDFileName: $! \n";

    _ProcessGPCmpdAbbrevs($CmpdAbbrevsRef, $AllowArbitraryChainAbbrev, $WriteSDFile, \*SDFILE);

    close SDFILE;
  }
  else {
    _ProcessGPCmpdAbbrevs($CmpdAbbrevsRef, $AllowArbitraryChainAbbrev, $WriteSDFile);
  }
}

# Process specified compound abbreviations containing any wild cards
# and count the number of valid abbreviations along with optional generation
# of a SD file.
#
sub _ProcessGPCmpdAbbrevs {
  my($CmpdAbbrevsRef, $AllowArbitraryChainAbbrev, $WriteSDFile, $SDFileRef) = @_;
  my($Abbrev, $AbbrevCount, $HeadGroupAbbrev, $ChainsAbbrev, $AbbrevModifier, $NewAbbrev, $Sn1Abbrev, $Sn2Abbrev, $AllowSubstituents, $AllowRings, @ChainsAbbrevList);

  $AbbrevCount = 0;

  ABBREV: for $Abbrev (@{$CmpdAbbrevsRef}) {
    $Sn1Abbrev = "0:0"; $Sn2Abbrev = "0:0";

    if (!($Abbrev =~ /\(/ && $Abbrev =~ /\)/)) {
      warn "Warning: Ignored GP compound abbreviation $Abbrev due to incorrect format. Valid format: HeadGroupAbbrev(Chain1Abbrev/Chain2Abbrev)\n";
      next ABBREV;
    }
    if (!ValidateGPAbbrev($Abbrev)) {
      next ABBREV;
    }

    ($HeadGroupAbbrev, $ChainsAbbrev, $AbbrevModifier) = ParseGPAbbrev($Abbrev);

    @ChainsAbbrevList = ();
    @ChainsAbbrevList = split /\//, $ChainsAbbrev;
    if (@ChainsAbbrevList != 2) {
      warn "Warning: Ignored GP compound Abbreviation $Abbrev due to incorrect format: Must contain chain abbreviations for both sn1 and sn2 chains only.\n";
      next ABBREV;
    }
    ($Sn1Abbrev, $Sn2Abbrev) = @ChainsAbbrevList;

    ($AllowSubstituents, $AllowRings) = (1, 0);
    if (!(ChainAbbrev::IsChainAbbrevOkay($Sn1Abbrev, $AllowSubstituents, $AllowRings, $AllowArbitraryChainAbbrev) && ChainAbbrev::IsChainAbbrevOkay($Sn2Abbrev, $AllowSubstituents, $AllowRings, $AllowArbitraryChainAbbrev))) {
      warn "Warning: Ignoring GP compound abbreviation $HeadGroupAbbrev($ChainsAbbrev)\n";
      next ABBREV;
    }
    if (!IsGPChainsAbbrevSupported($ChainsAbbrev, 1)) {
      warn "Warning: Ignoring GP compound abbreviation $HeadGroupAbbrev($ChainsAbbrev)\n";
      next ABBREV;
    }
    my($WildCardInHeadGroup) = ($HeadGroupAbbrev =~ /\*/ ) ? 1 : 0;

    if (!($WildCardInHeadGroup || ChainAbbrev::IsWildCardInChainAbbrev($Sn1Abbrev) || ChainAbbrev::IsWildCardInChainAbbrev($Sn2Abbrev))) {
      my($TemplateHeadGroup) = GetGPTemplateID($HeadGroupAbbrev, $ChainsAbbrev);
      if (exists $GPTemplatesDataMap{$TemplateHeadGroup}) {
	$NewAbbrev = "$HeadGroupAbbrev($ChainsAbbrev)$AbbrevModifier";
	$AbbrevCount++;
	if ($WriteSDFile) {
	  _GenerateAndWriteCmdDataString($SDFileRef, $NewAbbrev);
	}
      }
      else {
	warn "Warning: Ignored GP compound abbreviation $Abbrev : Abbreviation doesn't match any template\n";
      }
      next ABBREV;
    }

    # Arbitrary acyl chain abbreviation is not supported with wild cards...
    if ($AllowArbitraryChainAbbrev) {
      warn "Warning: Ignoring GP compound abbreviation $Abbrev : Allow arbitrary chain abbreviation option is not supported with wild cards in any part of the abbreviation\n";
      next ABBREV;
    }

    # Substituents are not supported with wild cards...
    if (ChainAbbrev::IsSubstituentInChainAbbrev($Sn1Abbrev) || ChainAbbrev::IsSubstituentInChainAbbrev($Sn2Abbrev)) {
      warn "Warning: Ignoring GP compound abbreviation $Abbrev : Substituent specifications are not supported with wild cards in any part of the abbreviation\n";
      next ABBREV;
    }

    # Get expanded abbreviation for each position...
    my($Sn1ExpandedAbbrevRef) = ChainAbbrev::ExpandChainAbbrev($Sn1Abbrev);
    my($Sn2ExpandedAbbrevRef) = ChainAbbrev::ExpandChainAbbrev($Sn2Abbrev);

    if ($Sn2Abbrev =~ /^(\*|\*:\*)$/i) {
      # Add 0:0 to Sn2Abbrev list containing wild cards...
      unshift(@{$Sn2ExpandedAbbrevRef}, "0:0")
    }

    # Enumerate various possibilities...
    my($ExpandedAbbrev, $ExpandedSn1Abbrev, $ExpandedSn2Abbrev);

    if ($WildCardInHeadGroup) {
      my($SupportedHeadGroupAbbrev);
      for $SupportedHeadGroupAbbrev (sort { $a cmp $b } keys %GPSupportedHeadGroupMap ) {
	SN1ABBREV: for $ExpandedSn1Abbrev (@$Sn1ExpandedAbbrevRef) {
	    if ($SupportedHeadGroupAbbrev =~ /^PPA$/i && $ExpandedSn1Abbrev =~ /^(O-|P-)/i ) {
	      next SN1ABBREV;
	    }
	    SN2ABBREV: for $ExpandedSn2Abbrev (@$Sn2ExpandedAbbrevRef) {
	      if (ChainAbbrev::IsAlkenylChainAbbrev($ExpandedSn2Abbrev)) {
		next SN2ABBREV;
	      }
	      if (ChainAbbrev::IsAlkylChainAbbrev($ExpandedSn2Abbrev) && !ChainAbbrev::IsAlkylChainAbbrev($ExpandedSn1Abbrev)) {
		next SN2ABBREV;
	      }
	      $ExpandedAbbrev = $SupportedHeadGroupAbbrev . '(' . $ExpandedSn1Abbrev . '/' . $ExpandedSn2Abbrev . ')' . $AbbrevModifier;
	      $AbbrevCount++;
	      if ($WriteSDFile) {
		_GenerateAndWriteCmdDataString($SDFileRef, $ExpandedAbbrev);
	      }
	    }
	  }
	}
    }
    else {
      SN1ABBREV: for $ExpandedSn1Abbrev (@$Sn1ExpandedAbbrevRef) {
	if ($HeadGroupAbbrev =~ /^PPA$/i && $ExpandedSn1Abbrev =~ /^(O-|P-)/i ) {
	  next SN1ABBREV;
	}
	SN2ABBREV: for $ExpandedSn2Abbrev (@$Sn2ExpandedAbbrevRef) {
	  if (ChainAbbrev::IsAlkenylChainAbbrev($ExpandedSn2Abbrev)) {
	    next  SN2ABBREV;
	  }
	  if (ChainAbbrev::IsAlkylChainAbbrev($ExpandedSn2Abbrev) && !ChainAbbrev::IsAlkylChainAbbrev($ExpandedSn1Abbrev)) {
	    next SN2ABBREV;
	  }
	  $ExpandedAbbrev = $HeadGroupAbbrev . '(' . $ExpandedSn1Abbrev . '/' . $ExpandedSn2Abbrev . ')' . $AbbrevModifier;
	  $AbbrevCount++;
	  if ($WriteSDFile) {
	    _GenerateAndWriteCmdDataString($SDFileRef, $ExpandedAbbrev);
	  }
	}
      }
    }
  }

  if ($AbbrevCount) {
    print "Valid abbreviations count: $AbbrevCount\n";
  }
  else {
    print "No valid abbreviations found...\n";
  }
}

# Generate a SD file for valid abbreviations...
sub GenerateSDFile {
  my($SDFileName, $CmdAbbrevsRef) = @_;
  my($Abbrev);

  open SDFILE, ">$SDFileName" or die "Error: Couldn't open $SDFileName: $! \n";

  for $Abbrev (@$CmdAbbrevsRef) {
    _GenerateAndWriteCmdDataString(\*SDFILE, $Abbrev);
  }

  close SDFILE;
}

# Generate appropriate compound data string containing structure data and
# write it out to SD file...
sub _GenerateAndWriteCmdDataString {
  my($SDFileRef, $Abbrev) = @_;
  my($CmpdDataString);

    $CmpdDataString = _GenerateCmdDataString($Abbrev);
    print $SDFileRef "$CmpdDataString\n";
}

# Generate appropriate compound data string containing structure data...
sub _GenerateCmdDataString {
  my($Abbrev) = @_;
  my($CmpdDataString, $CmpdAbbrevTemplateDataMapRef);

  $CmpdDataString = '';

  # Setup template data for a specific compound abbreviation...
  $CmpdAbbrevTemplateDataMapRef = SetupGPCmpdAbbrevTemplateDataMap($Abbrev);

  # Generate structure data for chains...
  my($Sn1AtomBlockLinesRef, $Sn1BondBlockLinesRef) = GPStr::GenerateGPChainStrData('Sn1', $CmpdAbbrevTemplateDataMapRef);
  my($Sn2AtomBlockLinesRef, $Sn2BondBlockLinesRef) = GPStr::GenerateGPChainStrData('Sn2', $CmpdAbbrevTemplateDataMapRef);
  my($Sn3AtomBlockLinesRef, $Sn3BondBlockLinesRef) = GPStr::GenerateGPChainStrData('Sn3', $CmpdAbbrevTemplateDataMapRef);

  # Generate data block lines including various desriptors...
  my($OntologyDataLinesRef) = GenerateCmpdOntologySDDataLines($CmpdAbbrevTemplateDataMapRef);

  # Setup first 4 SD file lines string: Name, MiscInfo, Comments, Count
  $CmpdDataString .= "$Abbrev\n";
  $CmpdDataString .= LMAPSStr::GenerateCmpdMiscInfoLine(). "\n";
  $CmpdDataString .= LMAPSStr::GenerateCmpdCommentsLine() . "\n";
  $CmpdDataString .= ChainStr::GenerateCmpdCountsLine($CmpdAbbrevTemplateDataMapRef) . "\n";

  # Atom lines for template and chains...
  $CmpdDataString .= ChainStr::GenerateAtomBlockLines($CmpdAbbrevTemplateDataMapRef, $Sn1AtomBlockLinesRef, $Sn2AtomBlockLinesRef, $Sn3AtomBlockLinesRef) . "\n";

  # Bond lines for template and chains...
  $CmpdDataString .= ChainStr::GenerateBondBlockLines($CmpdAbbrevTemplateDataMapRef, $Sn1BondBlockLinesRef, $Sn2BondBlockLinesRef, $Sn3BondBlockLinesRef) . "\n";

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

# Generate atom and bond block lines for a chain
sub GenerateGPChainStrData {
  my($ChainType, $CmpdAbbrevTemplateDataMapRef) = @_;

  return ChainStr::GenerateChainStrData($ChainType, $CmpdAbbrevTemplateDataMapRef);
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

  _SetupAbbrevAndSystematicName($CmpdAbbrevTemplateDataMapRef, \%OntologyDataMap);
  _SetupLMClassificationInfo($CmpdAbbrevTemplateDataMapRef, \%OntologyDataMap);
  _SetupChainLengthAndMultipleBondsCountInfo($CmpdAbbrevTemplateDataMapRef, \%OntologyDataMap);

  return \%OntologyDataMap;
}

# Abbreviation and systematic name...
sub _SetupAbbrevAndSystematicName {
  my($CmpdAbbrevTemplateDataMapRef, $OntologyDataMapRef) = @_;
  my($ChainsSystematicName, $HeadGroupAbbrev, $HeadGroupName, $AbbrevModifier, $ChainsAbbrev, @ChainIndices, @ChainPositions);

  $HeadGroupAbbrev = $CmpdAbbrevTemplateDataMapRef->{AbbrevID};
  $HeadGroupName = $CmpdAbbrevTemplateDataMapRef->{HeadGroupName};

  $AbbrevModifier = $CmpdAbbrevTemplateDataMapRef->{AbbrevModifier};

  @ChainIndices = (0, 1);
  @ChainPositions = (1, 2);

  # For S stereochemistry, HeadGroup position is 1 instead of 3...
  if ($AbbrevModifier =~ /^S$/i) {
    $HeadGroupName =~ s/^glycero-3-/glycero-1-/;
    @ChainPositions = (2, 3);
  }

  $ChainsSystematicName = ChainAbbrev::SetupChainsSystematicName($CmpdAbbrevTemplateDataMapRef, \@ChainIndices, \@ChainPositions);

  $OntologyDataMapRef->{'Systematic Name'} = "${ChainsSystematicName}-sn-${HeadGroupName}";

  if ($AbbrevModifier) {
    $AbbrevModifier = "[$AbbrevModifier]";
  }
  $ChainsAbbrev = ChainAbbrev::SetupChainsAbbrev($CmpdAbbrevTemplateDataMapRef, \@ChainIndices);

  $OntologyDataMapRef->{Abbrev} = "${HeadGroupAbbrev}(${ChainsAbbrev})${AbbrevModifier}";
}

# LM classification info...
sub _SetupLMClassificationInfo {
  my($CmpdAbbrevTemplateDataMapRef, $OntologyDataMapRef) = @_;

  if (IsEmpty($CmpdAbbrevTemplateDataMapRef->{LMCategory})) {
    return;
  }

  $OntologyDataMapRef->{'LM Category'} = $CmpdAbbrevTemplateDataMapRef->{LMCategory};
  $OntologyDataMapRef->{'LM Main Class'} = $CmpdAbbrevTemplateDataMapRef->{LMCategory} . $CmpdAbbrevTemplateDataMapRef->{LMMainClass};
  $OntologyDataMapRef->{'LM Sub Class'} = $CmpdAbbrevTemplateDataMapRef->{LMCategory} . $CmpdAbbrevTemplateDataMapRef->{LMMainClass} . $CmpdAbbrevTemplateDataMapRef->{LMSubClass};
}

# Multiple bond counts...
sub _SetupChainLengthAndMultipleBondsCountInfo {
  my($CmpdAbbrevTemplateDataMapRef, $OntologyDataMapRef) = @_;
  my($ChainIndex, $ChainAbbrev, $ChainCount, $ChainLengthAbbrev, $DoubleBondCountAbbrev, $ChainLength, $DoubleBonds, @ChainIndices);

  @ChainIndices = (0, 1);

  for $ChainIndex (@ChainIndices) {
    $ChainCount = $ChainIndex + 1;

    $ChainAbbrev = $CmpdAbbrevTemplateDataMapRef->{SnChainAdd}[$ChainIndex] ? $CmpdAbbrevTemplateDataMapRef->{SnAbbrev}[$ChainIndex] : '0:0';

    ($ChainLengthAbbrev, $DoubleBondCountAbbrev) = ChainAbbrev::ParseChainAbbrev($ChainAbbrev);

    $ChainLength = $ChainLengthAbbrev ? $CmpdAbbrevTemplateDataMapRef->{SnChainLength}[$ChainIndex] : '0';
    $DoubleBonds = $ChainLengthAbbrev ? $CmpdAbbrevTemplateDataMapRef->{SnDblBondCount}[$ChainIndex] : '0';

    $OntologyDataMapRef->{"Sn${ChainCount} Chain Length"} = $ChainLength;
    $OntologyDataMapRef->{"Sn${ChainCount} Double Bonds"} = $DoubleBonds;
  }
}

# Return reference to GPTemplatesDataMap...
sub GetGPTemplatesData {
  return \%GPTemplatesDataMap;
}

# Return reference to GPSupportedHeadGroupMap...
sub GetGPSupportedHeadGroupMap {
  return \%GPSupportedHeadGroupMap;
}

# Get GP template ID...
sub GetGPTemplateID {
  my($HeadGroupAbbrev, $ChainsAbbrev) = @_;
  my($HeadGroupID);

  $HeadGroupID = "";
  if ($HeadGroupAbbrev) {
    my(@AbbrevWords) = split /\//, $ChainsAbbrev;
    my ($Sn1Abbrev, $Sn2Abbrev) = @AbbrevWords;

    if ($Sn2Abbrev eq "0:0") {
      $HeadGroupID = (ChainAbbrev::IsAlkylChainAbbrev($Sn1Abbrev)) ? "Sn1Alkyl" : ((ChainAbbrev::IsAlkenylChainAbbrev($Sn1Abbrev)) ? "Sn1Alkenyl" : "Sn1" ) ;
      $HeadGroupID = "$HeadGroupAbbrev" . "$HeadGroupID";
    }
    else {
      my($Sn1ChainType, $Sn2ChainType);

      $Sn1ChainType = ChainAbbrev::IsAlkylChainAbbrev($Sn1Abbrev) ? "Sn1Alkyl" : (ChainAbbrev::IsAlkenylChainAbbrev($Sn1Abbrev) ? "Sn1Alkenyl" : "Sn1");
      $Sn2ChainType = ChainAbbrev::IsAlkylChainAbbrev($Sn2Abbrev) ? "Sn2Alkyl" : (ChainAbbrev::IsAlkenylChainAbbrev($Sn2Abbrev) ? "Sn2Alkenyl" : "Sn2");

      $HeadGroupID = "${HeadGroupAbbrev}${Sn1ChainType}${Sn2ChainType}";
    }
  }
  return $HeadGroupID;
}

# Does template exist to handle this abbreviation?
#
# Based on proposed LM subclasses, here is what's not allowed for GP.
#
# For GP: alkenyl at sn2; alkyl at sn2 without alkyl at sn1
#
#
sub IsGPChainsAbbrevSupported {
  my($Abbrev, $PrintWarning) = @_;
  my($Sn1Abbrev, $Sn2Abbrev, @AbbrevList);

  $Sn1Abbrev = "0:0"; $Sn2Abbrev = "0:0";
  @AbbrevList = split /\//, $Abbrev;
  if (@AbbrevList != 2) {
    if ($PrintWarning) {
      warn "Warning: Ignoring GP compound abbreviation $Abbrev : Must contain chain abbreviation for sn1 and sn2 only\n";
    }
    return 0;
  }
  if (ChainAbbrev::IsAlkenylChainAbbrev($Sn2Abbrev)) {
    if ($PrintWarning) {
      warn "Warning: Ignoring GP compound abbreviation $Abbrev : alkenyl chain at sn2 position.\n";
    }
    return 0;
  }
  if (ChainAbbrev::IsAlkylChainAbbrev($Sn2Abbrev) && !ChainAbbrev::IsAlkylChainAbbrev($Sn1Abbrev)) {
    if ($PrintWarning) {
      warn "Warning: Ignoring GP compound abbreviation $Abbrev : alkyl chain at sn2 position without alkyl chain at sn1.\n";
    }
    return 0;
  }

  return 1;
}

# Parse abbreviation...
sub ParseGPAbbrev {
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
# Check out the validity of GP abbreviation...
sub ValidateGPAbbrev {
  my($Abbrev) = @_;

  if (!($Abbrev =~ /\(/ && $Abbrev =~ /\)/)) {
    warn "Warning: Ignored GP compound abbreviation $Abbrev due to incorrect format. Valid format: HeadGroupAbbrev(Chain1Abbrev/Chain2Abbrev)\n";
    return 0;
  }

  # Make sure head group is these...
  if ($Abbrev =~ /^\(/) {
    warn "Warning: Ignored GP compound abbreviation $Abbrev due to incorrect format. Valid format: HeadGroupAbbrev(Chain1Abbrev/Chain2Abbrev)\n";
    return 0;
  }

  # Number of '(' and ')' parenthesis must match...
  my(@LeftParenthesis) = $Abbrev =~ /\(/g;
  my(@RightParenthesis) = $Abbrev =~ /\)/g;
  if ($#LeftParenthesis != $#RightParenthesis) {
    warn "Warning: Ignored GP compound abbreviation $Abbrev due to incorrect format: Number of ( and ) doesn't match. Valid format: HeadGroupAbbrev(Chain1Abbrev/Chain2Abbrev).\n";
    return 0;
  }

  @LeftParenthesis = $Abbrev =~ /\[/g;
  @RightParenthesis = $Abbrev =~ /\]/g;
  if ($#LeftParenthesis != $#RightParenthesis) {
    warn "Warning: Ignored GP compound abbreviation $Abbrev due to incorrect format: Number of [ and ] doesn't match. Valid format: HeadGroupAbbrev(Chain1Abbrev/Chain2Abbrev).\n";
    return 0;
  }

  my($HeadGroup, $ChainsAbbrev, $AbbrevModifier) = ParseGPAbbrev($Abbrev);

  if ($HeadGroup !~ /\*/ ) {
    if (!(exists $GPSupportedHeadGroupMap{$HeadGroup})) {
      warn "Warning: Ignored GP compound abbreviation $Abbrev) : Unknown head group $HeadGroup.\n";
      return 0;
    }
  }

  my(@ChainsAbbrevList) = ();
  @ChainsAbbrevList = split /\//, $ChainsAbbrev;
  if (@ChainsAbbrevList != 2) {
    warn "Warning: Ignored GP compound abbreviation $Abbrev due to incorrect format: Didn't specify abbreviations for  both sn1 and sn2 chains.\n";
    return 0;
  }

  my($Sn1Abbrev, $Sn2Abbrev) = @ChainsAbbrevList;

  if ($Sn1Abbrev eq "0:0") {
    warn "Warning: Ignoring GP compound abbreviation $Abbrev : sn1 abbreviation value of 0:0 is not allowed.\n";
    return 0;
  }

  if ($AbbrevModifier) {
    if (!($AbbrevModifier =~ /\[/ && $AbbrevModifier =~ /\]/)) {
      warn "Warning: Ignored GP compound abbreviation $Abbrev due to incorrect stereochemistry abbreviation: It must be enclosed by []\n";
      return 0;
    }
    $AbbrevModifier =~ s/(\[|\])//g;
    if ($AbbrevModifier !~ /^(rac|R|S|U)$/) {
      warn "Warning: Ignored GP compound abbreviation $Abbrev due to incorrect stereochemistry abbreviation: Valid values: rac, R, S, U\n";
      return 0;
    }
  }
  return 1;
}

# Set up template data for a specific cmpd abbreviation and return a
# reference to the hash...
sub SetupGPCmpdAbbrevTemplateDataMap {
  my($AbbrevHeadGroup, $Abbrev, $ChainsAbbrev, $AbbrevModifier, $Sn1Abbrev, $Sn2Abbrev, $Sn3Abbrev, $TemplateID, %AbbrevTemplateDataMap);

  ($Abbrev) = @_;
  %AbbrevTemplateDataMap = ();

  ($AbbrevHeadGroup, $ChainsAbbrev, $AbbrevModifier) = ParseGPAbbrev($Abbrev);

  #$ChainsAbbrev =~ s/\)//g;
  ($Sn1Abbrev, $Sn2Abbrev) = split /\//, $ChainsAbbrev;
  $Sn3Abbrev = '';

  my(@SnChainAdd) = ();
  $TemplateID = GetGPTemplateID($AbbrevHeadGroup, $ChainsAbbrev);

  if ($Sn2Abbrev eq "0:0") {
    push @SnChainAdd, (1,0,0);
  }
  else {
    push @SnChainAdd, (1,1,0);
  }
  @{$AbbrevTemplateDataMap{SnChainAdd}} = ();
  push @{$AbbrevTemplateDataMap{SnChainAdd}}, @SnChainAdd;

  @{$AbbrevTemplateDataMap{SnAbbrev}} = ();
  push @{$AbbrevTemplateDataMap{SnAbbrev}}, ($Sn1Abbrev, $Sn2Abbrev, $Sn3Abbrev);

  $AbbrevModifier =~ s/(\[|\])//g;
  $AbbrevTemplateDataMap{AbbrevModifier} = $AbbrevModifier;

  ChainStr::SetupTemplateDataMap('GP', \%AbbrevTemplateDataMap, $GPTemplatesDataMap{$TemplateID});

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
  %GPTemplatesDataMap = ();

  # GPCho templates...

  my($GPChoSn1Sn2TemplateString)=<<ENDGPCHOSN1SN2TEMPLATE;
GPCho sn1 acyl and sn2 acyl template structure
  LipdMAPS02060609152D

 23 22  0  0  0  0  0  0  0  0999 V2000
   -1.5990    0.3628    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.3133    0.7740    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0278    0.3628    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.7420    0.7740    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.7420    1.5997    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1861   -0.3515    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0119   -0.3515    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.4563    0.3628    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8845    0.7753    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1701    0.3628    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.5985    0.3456    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.3129   -0.0669    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0274    0.3456    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.7418   -0.0669    0.0000 N   0  3  0  0  0  0  0  0  0  0  0  0
    4.4563    0.3456    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.7418   -0.8919    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.4563   -0.4794    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.8464    0.6518    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    0.4825    0.0210    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.8464    1.4040    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.7602   -0.7737    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.7602   -1.5997    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.4745   -0.3610    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  2  0      
  4  8  1  0      
  1  7  1  6      
  1  6  1  1      
  9  1  1  0      
 10  9  1  0      
 12 11  1  0      
 13 12  1  0      
 14 13  1  0      
 15 14  1  0      
 14 16  1  0      
 14 17  1  0      
 18 11  1  0      
 18 19  1  0      
 18 20  2  0      
 21 22  2  0      
 21 23  1  0      
 21  7  1  0      
 18 10  1  0      
M  CHG  2  14   1  19  -1
M  END

ENDGPCHOSN1SN2TEMPLATE

  my($GPChoSn1AlkylSn2TemplateString)=<<ENDGPCHOSN1ALKYLSN2TEMPLATE;
GPCho sn1 alkyl and sn2 acyl template structure
  LipdMAPS02060609152D

 20 19  0  0  0  0  0  0  0  0999 V2000
   -2.0899    0.4606    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.8042    0.8718    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.5187    0.4606    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6770   -0.2537    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5028   -0.2537    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3754    0.8731    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6610    0.4606    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.1076    0.4435    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.8220    0.0309    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.5365    0.4435    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.2509    0.0309    0.0000 N   0  3  0  0  0  0  0  0  0  0  0  0
    3.9654    0.4435    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.2509   -0.7941    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.9654   -0.3815    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3555    0.7496    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0084    0.1188    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.3555    1.5019    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.2511   -0.6759    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.2511   -1.5019    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.9654   -0.2632    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  1  5  1  6      
  1  4  1  1      
  6  1  1  0      
  7  6  1  0      
  9  8  1  0      
 10  9  1  0      
 11 10  1  0      
 12 11  1  0      
 11 13  1  0      
 11 14  1  0      
 15  8  1  0      
 15 16  1  0      
 15 17  2  0      
 18 19  2  0      
 18 20  1  0      
 18  5  1  0      
 15  7  1  0      
M  CHG  2  11   1  16  -1
M  END

ENDGPCHOSN1ALKYLSN2TEMPLATE

  my($GPChoSn1TemplateString)=<<ENDGPCHOSN1TEMPLATE;
GPCho sn1 acyl template structure
  LipdMAPS02060609152D

 20 19  0  0  0  0  0  0  0  0999 V2000
   -1.5990    0.0089    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.3133    0.4201    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0278    0.0089    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.7420    0.4201    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.7420    1.2458    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1861   -0.7054    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0119   -0.7054    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.4563    0.0089    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8845    0.4214    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1701    0.0089    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.5985   -0.0083    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.3129   -0.4208    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0274   -0.0083    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.7418   -0.4208    0.0000 N   0  3  0  0  0  0  0  0  0  0  0  0
    4.4563   -0.0083    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.7418   -1.2458    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.4563   -0.8333    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.8464    0.2979    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    0.4825   -0.3329    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.8464    1.0501    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  2  0      
  4  8  1  0      
  1  7  1  6      
  1  6  1  1      
  9  1  1  0      
 10  9  1  0      
 12 11  1  0      
 13 12  1  0      
 14 13  1  0      
 15 14  1  0      
 14 16  1  0      
 14 17  1  0      
 18 11  1  0      
 18 19  1  0      
 18 20  2  0      
 18 10  1  0      
M  CHG  2  14   1  19  -1
M  END

ENDGPCHOSN1TEMPLATE

  my($GPChoSn1AlkylSn2AlkylTemplateString)=<<ENDGPCHOSN1ALKYLSN2ALKYLTEMPLATE;
GPCho sn1 alkyl and sn2 alkyl template structure
  LipdMAPS02060609152D

 18 17  0  0  0  0  0  0  0  0999 V2000
   -2.3133    0.1067    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0276    0.5179    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.7421    0.1067    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9004   -0.6076    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.7262   -0.6076    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5988    0.5192    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8844    0.1067    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.8843    0.0896    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.5987   -0.3230    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.3132    0.0896    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0276   -0.3230    0.0000 N   0  3  0  0  0  0  0  0  0  0  0  0
    3.7421    0.0896    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0276   -1.1480    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.7421   -0.7354    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.1322    0.3957    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2318   -0.2351    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.1322    1.1480    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.4745   -1.0298    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  1  5  1  6      
  1  4  1  1      
  6  1  1  0      
  7  6  1  0      
  9  8  1  0      
 10  9  1  0      
 11 10  1  0      
 12 11  1  0      
 11 13  1  0      
 11 14  1  0      
 15  8  1  0      
 15 16  1  0      
 15 17  2  0      
 18  5  1  0      
 15  7  1  0      
M  CHG  2  11   1  16  -1
M  END

ENDGPCHOSN1ALKYLSN2ALKYLTEMPLATE

  my($GPChoSn1AlkylTemplateString)=<<ENDGPCHOSN1ALKYLTEMPLATE;
GPCho sn1 alkyl template structure
  LipdMAPS02060609152D

 17 16  0  0  0  0  0  0  0  0999 V2000
   -2.3133    0.1067    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0276    0.5179    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.7421    0.1067    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9004   -0.6076    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.7262   -0.6076    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5988    0.5192    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8844    0.1067    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.8843    0.0896    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.5987   -0.3230    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.3132    0.0896    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0276   -0.3230    0.0000 N   0  3  0  0  0  0  0  0  0  0  0  0
    3.7421    0.0896    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0276   -1.1480    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.7421   -0.7354    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.1322    0.3957    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2318   -0.2351    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.1322    1.1480    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  1  5  1  6      
  1  4  1  1      
  6  1  1  0      
  7  6  1  0      
  9  8  1  0      
 10  9  1  0      
 11 10  1  0      
 12 11  1  0      
 11 13  1  0      
 11 14  1  0      
 15  8  1  0      
 15 16  1  0      
 15 17  2  0      
 15  7  1  0      
M  CHG  2  11   1  16  -1
M  END

ENDGPCHOSN1ALKYLTEMPLATE

  # GPEtn templates...

  my($GPEtnSn1Sn2TemplateString)=<<ENDGPETNSN1SN2TEMPLATE;
GPEtn sn1 acyl and sn2 acyl template structure
  LipdMAPS02060609152D

 20 19  0  0  0  0  0  0  0  0999 V2000
   -1.2418    0.3628    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9561    0.7740    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6706    0.3628    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.3848    0.7740    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.3848    1.5997    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8289   -0.3515    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6547   -0.3515    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.0990    0.3628    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5273    0.7753    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.1871    0.3628    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.9557    0.3456    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.6701   -0.0669    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.3846    0.3456    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.0990   -0.0669    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    1.2036    0.6518    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    0.8397    0.0210    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.2036    1.4040    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4030   -0.7737    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4030   -1.5997    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1173   -0.3610    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  2  0      
  4  8  1  0      
  1  7  1  6      
  1  6  1  1      
  9  1  1  0      
 10  9  1  0      
 12 11  1  0      
 13 12  1  0      
 14 13  1  0      
 15 11  1  0      
 15 16  1  0      
 15 17  2  0      
 18 19  2  0      
 18 20  1  0      
 18  7  1  0      
 15 10  1  0      
M  END

ENDGPETNSN1SN2TEMPLATE

  my($GPEtnSn1AlkylSn2TemplateString)=<<ENDGPETNSN1ALKYLSN2TEMPLATE;
GPEtn sn1 alkyl and sn2 acyl template structure
  LipdMAPS02060609152D

 17 16  0  0  0  0  0  0  0  0999 V2000
   -1.7328    0.4607    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4471    0.8719    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1617    0.4607    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3198   -0.2537    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1457   -0.2537    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0182    0.8732    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3038    0.4607    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.4650    0.4435    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.1794    0.0310    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.8940    0.4435    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.6084    0.0310    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.7128    0.7497    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    0.3489    0.1189    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.7128    1.5020    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.8941   -0.6759    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.8941   -1.5020    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.6084   -0.2632    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  1  5  1  6      
  1  4  1  1      
  6  1  1  0      
  7  6  1  0      
  9  8  1  0      
 10  9  1  0      
 11 10  1  0      
 12  8  1  0      
 12 13  1  0      
 12 14  2  0      
 15 16  2  0      
 15 17  1  0      
 15  5  1  0      
 12  7  1  0      
M  END

ENDGPETNSN1ALKYLSN2TEMPLATE

  my($GPEtnSn1TemplateString)=<<ENDGPETNSN1TEMPLATE;
GPEtn sn1 acyl template structure
  LipdMAPS02060609152D

 17 16  0  0  0  0  0  0  0  0999 V2000
   -1.2418   -0.2613    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9561    0.1499    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6706   -0.2613    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.3848    0.1499    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.3848    0.9756    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8289   -0.9756    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6547   -0.9756    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.0990   -0.2613    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5273    0.1512    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.1871   -0.2613    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.9557   -0.2785    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.6701   -0.6910    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.3846   -0.2785    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.0990   -0.6910    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    1.2036    0.0277    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    0.8397   -0.6031    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.2036    0.7799    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  2  0      
  4  8  1  0      
  1  7  1  6      
  1  6  1  1      
  9  1  1  0      
 10  9  1  0      
 12 11  1  0      
 13 12  1  0      
 14 13  1  0      
 15 11  1  0      
 15 16  1  0      
 15 17  2  0      
 15 10  1  0      
M  END

ENDGPETNSN1TEMPLATE

  my($GPEtnSn1AlkylSn2AlkylTemplateString)=<<ENDGPETNSN1ALKYLSN2ALKYLTEMPLATE;
GPEtn sn1 alkyl and sn2 alkyl template structure
  LipdMAPS02060609152D

 15 14  0  0  0  0  0  0  0  0999 V2000
   -1.9563    0.0477    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6706    0.4589    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.3852    0.0477    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5432   -0.6668    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.3692   -0.6668    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2416    0.4602    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5272    0.0477    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.2417    0.0305    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.9562   -0.3821    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.6708    0.0305    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.3852   -0.3821    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.4895    0.3367    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    0.1256   -0.2942    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.4895    1.0890    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1176   -1.0890    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  1  5  1  6      
  1  4  1  1      
  6  1  1  0      
  7  6  1  0      
  9  8  1  0      
 10  9  1  0      
 11 10  1  0      
 12  8  1  0      
 12 13  1  0      
 12 14  2  0      
 15  5  1  0      
 12  7  1  0      
M  END

ENDGPETNSN1ALKYLSN2ALKYLTEMPLATE

  my($GPEtnSn1AlkylTemplateString)=<<ENDGPETNSN1ALKYLTEMPLATE;
GPEtn sn1 alkyl template structure
  LipdMAPS02060609152D

 14 13  0  0  0  0  0  0  0  0999 V2000
   -1.9563   -0.1635    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6706    0.2478    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.3852   -0.1635    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5432   -0.8779    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.3692   -0.8779    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2416    0.2491    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5272   -0.1635    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.2417   -0.1807    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.9562   -0.5932    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.6708   -0.1807    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.3852   -0.5932    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.4895    0.1256    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    0.1256   -0.5053    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.4895    0.8779    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  1  5  1  6      
  1  4  1  1      
  6  1  1  0      
  7  6  1  0      
  9  8  1  0      
 10  9  1  0      
 11 10  1  0      
 12  8  1  0      
 12 13  1  0      
 12 14  2  0      
 12  7  1  0      
M  END

ENDGPETNSN1ALKYLTEMPLATE

  # GPSer templates...

  my($GPSerSn1Sn2TemplateString)=<<ENDGPSERSN1SN2TEMPLATE;
GPSer sn1 acyl and sn2 acyl template structure
  LipdMAPS02060609152D

 24 23  0  0  0  0  0  0  0  0999 V2000
    3.7162    0.9186    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4815    0.3124    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1959    0.7236    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.9104    0.3124    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.6247    0.7236    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.6247    1.5493    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0686   -0.4020    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8945   -0.4020    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.3390    0.3124    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7670    0.7249    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0525    0.3124    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.7162    0.2952    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.4307   -0.1174    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.1452    0.2952    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.8597   -0.1174    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.9641    0.6014    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    0.6001   -0.0295    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.9641    1.3536    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6428   -0.8242    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6428   -1.6503    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.3572   -0.4115    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.8069    0.8812    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.7162    1.6503    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.3390    0.5590    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  2  3  1  0      
  3  4  1  0      
  4  5  1  0      
  5  6  2  0      
  5  9  1  0      
  2  8  1  6      
  2  7  1  1      
 10  2  1  0      
 11 10  1  0      
 13 12  1  0      
 14 13  1  0      
 15 14  1  0      
 16 12  1  0      
 16 17  1  0      
 16 18  2  0      
 19 20  2  0      
 19 21  1  0      
 19  8  1  0      
 16 11  1  0      
 14 22  1  6      
 14  1  1  1      
  1 23  2  0      
  1 24  1  0      
M  END

ENDGPSERSN1SN2TEMPLATE

  my($GPSerSn1AlkylSn2TemplateString)=<<ENDGPSERSN1ALKYLSN2TEMPLATE;
GPSer sn1 alkyl and sn2 acyl template structure
  LipdMAPS02060609152D

 21 20  0  0  0  0  0  0  0  0999 V2000
    3.2252    0.9186    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9723    0.3124    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6867    0.7236    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.4012    0.3124    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5595   -0.4020    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.3853   -0.4020    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2579    0.7249    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5434    0.3124    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.2253    0.2952    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.9397   -0.1174    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.6542    0.2952    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.3687   -0.1174    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.4732    0.6014    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    0.1092   -0.0295    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.4732    1.3536    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1336   -0.8242    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1336   -1.6502    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.8480   -0.4115    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.3159    0.8812    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.2252    1.6502    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.8480    0.5590    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  2  3  1  0      
  3  4  1  0      
  2  6  1  6      
  2  5  1  1      
  7  2  1  0      
  8  7  1  0      
 10  9  1  0      
 11 10  1  0      
 12 11  1  0      
 13  9  1  0      
 13 14  1  0      
 13 15  2  0      
 16 17  2  0      
 16 18  1  0      
 16  6  1  0      
 13  8  1  0      
 11 19  1  6      
 11  1  1  1      
  1 20  2  0      
  1 21  1  0      
M  END

ENDGPSERSN1ALKYLSN2TEMPLATE

  my($GPSerSn1TemplateString)=<<ENDGPSERSN1TEMPLATE;
GPSer sn1 acyl template structure
  LipdMAPS02060609152D

 21 20  0  0  0  0  0  0  0  0999 V2000
    3.7161    0.2944    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4815   -0.3117    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1958    0.0994    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.9103   -0.3117    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.6246    0.0994    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.6246    0.9251    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0686   -1.0261    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8944   -1.0261    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.3389   -0.3117    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7670    0.1007    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0525   -0.3117    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.7161   -0.3289    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.4306   -0.7415    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.1451   -0.3289    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.8596   -0.7415    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.9641   -0.0227    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    0.6001   -0.6536    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.9641    0.7294    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.8068    0.2570    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.7161    1.0261    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.3389   -0.0651    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  2  3  1  0      
  3  4  1  0      
  4  5  1  0      
  5  6  2  0      
  5  9  1  0      
  2  8  1  6      
  2  7  1  1      
 10  2  1  0      
 11 10  1  0      
 13 12  1  0      
 14 13  1  0      
 15 14  1  0      
 16 12  1  0      
 16 17  1  0      
 16 18  2  0      
 16 11  1  0      
 14 19  1  6      
 14  1  1  1      
  1 20  2  0      
  1 21  1  0      
M  END

ENDGPSERSN1TEMPLATE

  my($GPSerSn1AlkylSn2AlkylTemplateString)=<<ENDGPSERSN1ALKYLSN2ALKYLTEMPLATE;
GPSer sn1 alkyl and sn2 alkyl template structure
  LipdMAPS02060609152D

 19 18  0  0  0  0  0  0  0  0999 V2000
    3.0017    0.5056    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1956   -0.1006    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.9100    0.3106    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.6245   -0.1006    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7828   -0.8150    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6086   -0.8150    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4813    0.3119    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7668   -0.1006    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.0019   -0.1178    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.7162   -0.5304    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.4307   -0.1178    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.1452   -0.5304    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.2498    0.1884    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1142   -0.4425    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.2498    0.9406    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.3569   -1.2372    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0924    0.4682    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.0017    1.2372    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.6245    0.1460    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  2  3  1  0      
  3  4  1  0      
  2  6  1  6      
  2  5  1  1      
  7  2  1  0      
  8  7  1  0      
 10  9  1  0      
 11 10  1  0      
 12 11  1  0      
 13  9  1  0      
 13 14  1  0      
 13 15  2  0      
 16  6  1  0      
 13  8  1  0      
 11 17  1  6      
 11  1  1  1      
  1 18  2  0      
  1 19  1  0      
M  END

ENDGPSERSN1ALKYLSN2ALKYLTEMPLATE

  my($GPSerSn1AlkylTemplateString)=<<ENDGPSERSN1ALKYLTEMPLATE;
GPSer sn1 alkyl template structure
  LipdMAPS02060609152D

 18 17  0  0  0  0  0  0  0  0999 V2000
    3.0017    0.2945    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1956   -0.3117    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.9100    0.0995    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.6245   -0.3117    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7828   -1.0261    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6086   -1.0261    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4813    0.1008    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7668   -0.3117    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.0019   -0.3289    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.7162   -0.7415    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.4307   -0.3289    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.1452   -0.7415    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.2498   -0.0227    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1142   -0.6536    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.2498    0.7295    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.0924    0.2571    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.0017    1.0261    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.6245   -0.0651    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  2  3  1  0      
  3  4  1  0      
  2  6  1  6      
  2  5  1  1      
  7  2  1  0      
  8  7  1  0      
 10  9  1  0      
 11 10  1  0      
 12 11  1  0      
 13  9  1  0      
 13 14  1  0      
 13 15  2  0      
 13  8  1  0      
 11 16  1  6      
 11  1  1  1      
  1 17  2  0      
  1 18  1  0      
M  END

ENDGPSERSN1ALKYLTEMPLATE

  # GPGro templates...

  my($GPGroSn1Sn2TemplateString)=<<ENDGPGROSN1SN2TEMPLATE;
GPGro sn1 acyl and sn2 acyl template structure
  LipdMAPS02060609152D

 23 22  0  0  0  0  0  0  0  0999 V2000
    3.5985    0.9691    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5991    0.3629    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.3134    0.7741    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0279    0.3629    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.7422    0.7741    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.7422    1.5998    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1862   -0.3515    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0121   -0.3515    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.4565    0.3629    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8846    0.7754    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1701    0.3629    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.5985    0.3457    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.3130   -0.0669    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0275    0.3457    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.7420   -0.0669    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.8465    0.6519    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    0.4825    0.0210    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.8465    1.4041    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.7603   -0.7737    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.7603   -1.5998    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.4747   -0.3610    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.6892    0.9317    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.4565    0.3455    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  2  3  1  0      
  3  4  1  0      
  4  5  1  0      
  5  6  2  0      
  5  9  1  0      
  2  8  1  6      
  2  7  1  1      
 10  2  1  0      
 11 10  1  0      
 13 12  1  0      
 14 13  1  0      
 15 14  1  0      
 16 12  1  0      
 16 17  1  0      
 16 18  2  0      
 19 20  2  0      
 19 21  1  0      
 19  8  1  0      
 16 11  1  0      
 14 22  1  1      
 14  1  1  6      
 15 23  1  0      
M  END

ENDGPGROSN1SN2TEMPLATE

  my($GPGroSn1AlkylSn2TemplateString)=<<ENDGPGROSN1ALKYLSN2TEMPLATE;
GPGro sn1 alkyl and sn2 acyl template structure
  LipdMAPS02060609152D

 20 19  0  0  0  0  0  0  0  0999 V2000
   -2.0451    0.4676    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.7594    0.8788    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.4738    0.4676    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6322   -0.2469    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4580   -0.2469    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3305    0.8801    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6160    0.4676    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.0830    0.4367    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.7975    0.0242    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.5120    0.4367    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.2265    0.0242    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3310    0.7429    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0330    0.1120    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.3310    1.4950    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.2063   -0.6692    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.2063   -1.4950    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.9206   -0.2564    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.9206    0.4275    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.1737    1.0341    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.0229    1.0098    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  6  1  1  0      
  7  6  1  0      
  9  8  1  0      
 10  9  1  0      
 11 10  1  0      
 12  8  1  0      
 12 13  1  0      
 12 14  2  0      
 15 16  2  0      
 15 17  1  0      
 15  5  1  0      
 11 18  1  0      
 10 19  1  1      
 10 20  1  6      
  1  5  1  6      
  1  4  1  1      
 12  7  1  0      
M  END

ENDGPGROSN1ALKYLSN2TEMPLATE

  my($GPGroSn1TemplateString)=<<ENDGPGROSN1TEMPLATE;
GPGro sn1 acyl template structure
  LipdMAPS02060609152D

 20 19  0  0  0  0  0  0  0  0999 V2000
    3.5986    0.3449    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5991   -0.2612    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.3134    0.1499    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0280   -0.2612    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.7423    0.1499    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.7423    0.9756    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1862   -0.9756    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0120   -0.9756    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.4566   -0.2612    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8846    0.1512    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1701   -0.2612    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.5986   -0.2784    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.3131   -0.6910    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0276   -0.2784    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.7421   -0.6910    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.8466    0.0278    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    0.4825   -0.6031    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.8466    0.7799    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.6893    0.3075    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.4566   -0.2785    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  2  3  1  0      
  3  4  1  0      
  4  5  1  0      
  5  6  2  0      
  5  9  1  0      
  2  8  1  6      
  2  7  1  1      
 10  2  1  0      
 11 10  1  0      
 13 12  1  0      
 14 13  1  0      
 15 14  1  0      
 16 12  1  0      
 16 17  1  0      
 16 18  2  0      
 16 11  1  0      
 14 19  1  1      
 14  1  1  6      
 15 20  1  0      
M  END

ENDGPGROSN1TEMPLATE

  my($GPGroSn1AlkylSn2AlkylTemplateString)=<<ENDGPGROSN1ALKYLSN2ALKYLTEMPLATE;
GPGro sn1 alkyl and sn2 alkyl template structure
  LipdMAPS02060609152D

 18 17  0  0  0  0  0  0  0  0999 V2000
   -2.2684    0.0547    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.9827    0.4659    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.6971    0.0547    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8555   -0.6598    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6813   -0.6598    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5539    0.4672    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8394    0.0547    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.8596    0.0238    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.5741   -0.3887    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.2885    0.0238    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0030   -0.3887    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.1076    0.3300    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2564   -0.3009    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.1076    1.0821    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.4296   -1.0821    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.6971    0.0146    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.9502    0.6212    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.7994    0.5969    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  6  1  1  0      
  7  6  1  0      
  9  8  1  0      
 10  9  1  0      
 11 10  1  0      
 12  8  1  0      
 12 13  1  0      
 12 14  2  0      
 15  5  1  0      
 11 16  1  0      
 10 17  1  1      
 10 18  1  6      
  1  5  1  6      
  1  4  1  1      
 12  7  1  0      
M  END

ENDGPGROSN1ALKYLSN2ALKYLTEMPLATE

  my($GPGroSn1AlkylTemplateString)=<<ENDGPGROSN1ALKYLTEMPLATE;
GPGro sn1 alkyl template structure
  LipdMAPS02060609152D

 17 16  0  0  0  0  0  0  0  0999 V2000
   -2.2684   -0.1564    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.9827    0.2547    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.6971   -0.1564    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8555   -0.8709    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6813   -0.8709    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5539    0.2560    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8394   -0.1564    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.8596   -0.1873    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.5741   -0.5998    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.2885   -0.1873    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0030   -0.5998    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.1076    0.1188    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2564   -0.5120    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.1076    0.8709    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.6971   -0.1965    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.9502    0.4100    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.7994    0.3857    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  6  1  1  0      
  7  6  1  0      
  9  8  1  0      
 10  9  1  0      
 11 10  1  0      
 12  8  1  0      
 12 13  1  0      
 12 14  2  0      
 11 15  1  0      
 10 16  1  1      
 10 17  1  6      
  1  5  1  6      
  1  4  1  1      
 12  7  1  0      
M  END

ENDGPGROSN1ALKYLTEMPLATE

  # GPGroP templates...

  my($GPGroPSn1Sn2TemplateString)=<<ENDGPGROPSN1SN2TEMPLATE;
GPGroP sn1 acyl and sn2 acyl template structure
  LipdMAPS02060609152D

 27 26  0  0  0  0  0  0  0  0999 V2000
    2.8014    0.9691    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.3959    0.3629    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1102    0.7741    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.8246    0.3629    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.5389    0.7741    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.5389    1.5997    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9830   -0.3515    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.8089   -0.3515    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.2532    0.3629    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6814    0.7754    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9670    0.3629    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.8015    0.3457    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.5160   -0.0669    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.2305    0.3457    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.9449   -0.0669    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0496    0.6519    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3144    0.0210    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.0496    1.4040    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.5571   -0.7737    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.5571   -1.5997    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.2714   -0.3610    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.8922    0.9317    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.6594    0.3455    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.4563    0.5590    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    4.4563    1.3840    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.2532    0.3455    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.0438   -0.1555    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  2  3  1  0      
  3  4  1  0      
  4  5  1  0      
  5  6  2  0      
  5  9  1  0      
  2  8  1  6      
  2  7  1  1      
 10  2  1  0      
 11 10  1  0      
 13 12  1  0      
 14 13  1  0      
 15 14  1  0      
 16 12  1  0      
 16 17  1  0      
 16 18  2  0      
 19 20  2  0      
 19 21  1  0      
 19  8  1  0      
 16 11  1  0      
 14 22  1  1      
 14  1  1  6      
 15 23  1  0      
 23 24  1  0      
 24 25  2  0      
 24 26  1  0      
 24 27  1  0      
M  END

ENDGPGROPSN1SN2TEMPLATE

  my($GPGroPSn1AlkylSn2TemplateString)=<<ENDGPGROPSN1ALKYLSN2TEMPLATE;
GPGroP sn1 alkyl and sn2 acyl template structure
  LipdMAPS02060609152D

 24 23  0  0  0  0  0  0  0  0999 V2000
   -2.8419    0.4676    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.5562    0.8788    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.2706    0.4676    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4290   -0.2469    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -3.2548   -0.2469    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1273    0.8801    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4129    0.4676    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.2861    0.4367    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.0006    0.0242    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.7150    0.4367    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.4295    0.0242    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4659    0.7429    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8299    0.1120    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4659    1.4950    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.0031   -0.6692    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.0031   -1.4950    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.7174   -0.2564    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.1236    0.4275    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.3767    1.0341    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.2259    1.0098    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.9205    0.6410    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    3.9205    1.4660    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.5080   -0.0735    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.7174    0.4275    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  6  1  1  0      
  7  6  1  0      
  9  8  1  0      
 10  9  1  0      
 11 10  1  0      
 12  8  1  0      
 12 13  1  0      
 12 14  2  0      
 15 16  2  0      
 15 17  1  0      
 15  5  1  0      
 11 18  1  0      
 10 19  1  1      
 10 20  1  6      
  1  5  1  6      
  1  4  1  1      
 12  7  1  0      
 18 21  1  0      
 21 22  2  0      
 21 23  1  0      
 21 24  1  0      
M  END

ENDGPGROPSN1ALKYLSN2TEMPLATE

  my($GPGroPSn1TemplateString)=<<ENDGPGROPSN1TEMPLATE;
GPGroP sn1 acyl template structure
  LipdMAPS02060609152D

 24 23  0  0  0  0  0  0  0  0999 V2000
    2.8018    0.3449    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.3960   -0.2612    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1104    0.1499    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.8250   -0.2612    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.5393    0.1499    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.5393    0.9756    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9831   -0.9756    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.8090   -0.9756    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.2536   -0.2612    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6815    0.1512    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9670   -0.2612    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.8018   -0.2784    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.5163   -0.6910    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.2308   -0.2784    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.9453   -0.6910    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0497    0.0278    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3144   -0.6031    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.0497    0.7799    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.8925    0.3075    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.6598   -0.2785    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.4567   -0.0650    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    4.4567    0.7600    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.0442   -0.7795    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.2536   -0.2785    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  2  3  1  0      
  3  4  1  0      
  4  5  1  0      
  5  6  2  0      
  5  9  1  0      
  2  8  1  6      
  2  7  1  1      
 10  2  1  0      
 11 10  1  0      
 13 12  1  0      
 14 13  1  0      
 15 14  1  0      
 16 12  1  0      
 16 17  1  0      
 16 18  2  0      
 16 11  1  0      
 14 19  1  1      
 14  1  1  6      
 15 20  1  0      
 20 21  1  0      
 21 22  2  0      
 21 23  1  0      
 21 24  1  0      
M  END

ENDGPGROPSN1TEMPLATE

  my($GPGroPSn1AlkylSn2AlkylTemplateString)=<<ENDGPGROPSN1ALKYLSN2ALKYLTEMPLATE;
GPGroP sn1 alkyl and sn2 alkyl template structure
  LipdMAPS02060609152D

 22 21  0  0  0  0  0  0  0  0999 V2000
   -3.0655    0.0547    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.7798    0.4659    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.4943    0.0547    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6525   -0.6598    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -3.4784   -0.6598    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.3509    0.4672    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6364    0.0547    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.0628    0.0238    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.7773   -0.3887    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4918    0.0238    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.2063   -0.3887    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6893    0.3300    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0533   -0.3009    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6893    1.0822    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.2267   -1.0822    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.9005    0.0146    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.1535    0.6212    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.0027    0.5969    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.6974    0.2281    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    3.6974    1.0531    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.2849   -0.4863    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.4943    0.0146    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  6  1  1  0      
  7  6  1  0      
  9  8  1  0      
 10  9  1  0      
 11 10  1  0      
 12  8  1  0      
 12 13  1  0      
 12 14  2  0      
 15  5  1  0      
 11 16  1  0      
 10 17  1  1      
 10 18  1  6      
  1  5  1  6      
  1  4  1  1      
 12  7  1  0      
 16 19  1  0      
 19 20  2  0      
 19 21  1  0      
 19 22  1  0      
M  END

ENDGPGROPSN1ALKYLSN2ALKYLTEMPLATE

  my($GPGroPSn1AlkylTemplateString)=<<ENDGPGROPSN1ALKYLTEMPLATE;
GPGroP sn1 alkyl template structure
  LipdMAPS02060609152D

 21 20  0  0  0  0  0  0  0  0999 V2000
   -3.0655   -0.1564    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.7798    0.2547    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.4943   -0.1564    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6525   -0.8710    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -3.4784   -0.8710    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.3509    0.2560    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6364   -0.1564    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.0628   -0.1873    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.7773   -0.5998    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4918   -0.1873    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.2063   -0.5998    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6893    0.1188    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0533   -0.5120    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6893    0.8710    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.9005   -0.1965    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.1535    0.4100    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.0027    0.3857    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.6974    0.0170    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    3.6974    0.8420    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.2849   -0.6975    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.4943   -0.1965    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  6  1  1  0      
  7  6  1  0      
  9  8  1  0      
 10  9  1  0      
 11 10  1  0      
 12  8  1  0      
 12 13  1  0      
 12 14  2  0      
 11 15  1  0      
 10 16  1  1      
 10 17  1  6      
  1  5  1  6      
  1  4  1  1      
 12  7  1  0      
 15 18  1  0      
 18 19  2  0      
 18 20  1  0      
 18 21  1  0      
M  END

ENDGPGROPSN1ALKYLTEMPLATE

  # GPIns templates...

  my($GPInsSn1Sn2TemplateString)=<<ENDGPINSSN1SN2TEMPLATE;
GPIns sn1 acyl and sn2 template structure
  LipdMAPS02060609152D

 28 28  0  0  0  0  0  0  0  0999 V2000
   -2.6495    0.3628    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.3639    0.7740    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.0783    0.3628    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.7927    0.7740    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.7927    1.5998    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2367   -0.3516    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0624   -0.3516    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.5069    0.3628    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9350    0.7753    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2206    0.3628    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.7284    0.3414    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0237    0.6476    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3876    0.0168    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0237    1.3998    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.8108   -0.7738    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.8108   -1.5998    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.5250   -0.3611    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.5053    0.6654    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.1942    1.0156    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5132   -0.1695    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.8259    0.1841    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.1405   -0.1695    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.8182    1.0156    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.3973    0.8021    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.8219    1.2299    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.9865    1.2623    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.0540    0.1300    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.5069    0.8417    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  2  0      
  4  8  1  0      
  9  1  1  0      
 10  9  1  0      
 12 11  1  0      
 12 13  1  0      
 12 14  2  0      
 15 16  2  0      
 15 17  1  0      
 15  7  1  0      
 18 19  1  0      
 18 23  1  0      
 19 20  1  0      
 22 23  1  0      
 22 21  1  1      
 20 21  1  1      
 19 24  1  0      
 11 20  1  0      
 21 25  1  0      
 18 26  1  0      
 22 27  1  0      
 23 28  1  0      
 10 12  1  0      
  1  7  1  6      
  1  6  1  1      
M  END

ENDGPINSSN1SN2TEMPLATE

  my($GPInsSn1AlkylSn2TemplateString)=<<ENDGPINSSN1ALKYLSN2TEMPLATE;
GPIns sn1 alkyl template structure
  LipdMAPS02060609152D

 25 25  0  0  0  0  0  0  0  0999 V2000
   -3.1397    0.4627    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.8540    0.8738    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.5682    0.4627    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.7270   -0.2516    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -3.5526   -0.2516    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4253    0.8751    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7112    0.4627    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.2374    0.4413    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5146    0.7475    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8783    0.1168    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5146    1.4994    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.3007   -0.6737    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.3007   -1.4994    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.0148   -0.2611    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0136    0.7653    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.7029    1.1153    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0220   -0.0695    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.3345    0.2841    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.6487   -0.0695    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.3262    1.1153    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.9061    0.9018    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.3305    1.3296    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.4948    1.3620    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.5620    0.2300    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.0148    0.9414    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  6  1  1  0      
  7  6  1  0      
  9  8  1  0      
  9 10  1  0      
  9 11  2  0      
 12 13  2  0      
 12 14  1  0      
 12  5  1  0      
 15 16  1  0      
 15 20  1  0      
 16 17  1  0      
 19 20  1  0      
 19 18  1  1      
 17 18  1  1      
 16 21  1  0      
  8 17  1  0      
 18 22  1  0      
 15 23  1  0      
 19 24  1  0      
 20 25  1  0      
  7  9  1  0      
  1  5  1  6      
  1  4  1  1      
M  END

ENDGPINSSN1ALKYLSN2TEMPLATE

  my($GPInsSn1TemplateString)=<<ENDGPINSSN1TEMPLATE;
GPIns sn1 acyl template structure
  LipdMAPS02060609152D

 25 25  0  0  0  0  0  0  0  0999 V2000
   -2.6495   -0.2612    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.3639    0.1500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.0783   -0.2612    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.7927    0.1500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.7927    0.9757    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2367   -0.9757    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0624   -0.9757    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.5069   -0.2612    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9350    0.1513    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2206   -0.2612    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.7284   -0.2827    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0237    0.0236    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3876   -0.6073    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0237    0.7757    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.5053    0.0413    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.1942    0.3916    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5132   -0.7935    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.8259   -0.4399    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.1405   -0.7935    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.8182    0.3916    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.3973    0.1781    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.8219    0.6059    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.9865    0.6382    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.0540   -0.4940    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.5069    0.2177    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  2  0      
  4  8  1  0      
  9  1  1  0      
 10  9  1  0      
 12 11  1  0      
 12 13  1  0      
 12 14  2  0      
 15 16  1  0      
 15 20  1  0      
 16 17  1  0      
 19 20  1  0      
 19 18  1  1      
 17 18  1  1      
 16 21  1  0      
 11 17  1  0      
 18 22  1  0      
 15 23  1  0      
 19 24  1  0      
 20 25  1  0      
 10 12  1  0      
  1  7  1  6      
  1  6  1  1      
M  END

ENDGPINSSN1TEMPLATE

  my($GPInsSn1AlkylSn2AlkylTemplateString)=<<ENDGPINSSN1ALKYLSN2ALKYLTEMPLATE;
GPIns sn1 alkyl and sn2 alkyl template structure
  LipdMAPS02060609152D

 23 23  0  0  0  0  0  0  0  0999 V2000
   -3.3628    0.0498    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.0770    0.4609    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.7912    0.0498    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.9501   -0.6644    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -3.7757   -0.6644    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6484    0.4622    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9344    0.0498    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.0141    0.0284    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7379    0.3346    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1015   -0.2960    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7379    1.0865    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.5237   -1.0865    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.7901    0.3524    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4795    0.7024    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7986   -0.4823    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.1111   -0.1287    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.4252   -0.4823    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.1026    0.7024    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6828    0.4889    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.1071    0.9167    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.2713    0.9491    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.3384   -0.1828    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.7912    0.5285    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  6  1  1  0      
  7  6  1  0      
  9  8  1  0      
  9 10  1  0      
  9 11  2  0      
 12  5  1  0      
 13 14  1  0      
 13 18  1  0      
 14 15  1  0      
 17 18  1  0      
 17 16  1  1      
 15 16  1  1      
 14 19  1  0      
  8 15  1  0      
 16 20  1  0      
 13 21  1  0      
 17 22  1  0      
 18 23  1  0      
  7  9  1  0      
  1  5  1  6      
  1  4  1  1      
M  END

ENDGPINSSN1ALKYLSN2ALKYLTEMPLATE

  my($GPInsSn1AlkylTemplateString)=<<ENDGPINSSN1ALKYLTEMPLATE;
GPIns sn1 alkyl template structure
  LipdMAPS02060609152D

 22 22  0  0  0  0  0  0  0  0999 V2000
   -3.3628   -0.1612    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.0770    0.2499    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.7912   -0.1612    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.9501   -0.8754    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -3.7757   -0.8754    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6484    0.2512    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9344   -0.1612    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.0141   -0.1826    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7379    0.1236    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1015   -0.5071    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7379    0.8754    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.7901    0.1414    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4795    0.4914    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7986   -0.6934    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.1111   -0.3398    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.4252   -0.6934    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.1026    0.4914    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6828    0.2779    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.1071    0.7057    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.2713    0.7381    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.3384   -0.3939    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.7912    0.3175    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  6  1  1  0      
  7  6  1  0      
  9  8  1  0      
  9 10  1  0      
  9 11  2  0      
 12 13  1  0      
 12 17  1  0      
 13 14  1  0      
 16 17  1  0      
 16 15  1  1      
 14 15  1  1      
 13 18  1  0      
  8 14  1  0      
 15 19  1  0      
 12 20  1  0      
 16 21  1  0      
 17 22  1  0      
  7  9  1  0      
  1  5  1  6      
  1  4  1  1      
M  END

ENDGPINSSN1ALKYLTEMPLATE

  # GPInsP templates...

  my($GPInsPSn1Sn2TemplateString)=<<ENDGPINSPSN1SN2TEMPLATE;
GPInsP sn1 acyl and sn2 template structure
  LipdMAPS02060609152D

 32 32  0  0  0  0  0  0  0  0999 V2000
   -3.3858    0.3627    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.1001    0.7739    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.8144    0.3627    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.5287    0.7739    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.5287    1.5995    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.9731   -0.3515    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -3.7987   -0.3515    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -6.2427    0.3627    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6714    0.7752    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9572    0.3627    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0085    0.3413    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7605    0.6475    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1243    0.0168    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7605    1.3996    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.5469   -0.7737    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.5469   -1.5995    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.2610   -0.3610    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.7679    0.6653    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4571    1.0154    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7762   -0.1695    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0887    0.1841    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.4030   -0.1695    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.0806    1.0154    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6603    0.8020    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.0847    1.2297    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.2491    1.2621    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.3164    0.1300    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.7692    0.8416    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.4459    0.4326    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    5.4459    1.2576    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.0334   -0.2818    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    6.2427    0.2191    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  2  0      
  4  8  1  0      
  9  1  1  0      
 10  9  1  0      
 12 11  1  0      
 12 13  1  0      
 12 14  2  0      
 15 16  2  0      
 15 17  1  0      
 15  7  1  0      
 18 19  1  0      
 18 23  1  0      
 19 20  1  0      
 22 23  1  0      
 22 21  1  1      
 20 21  1  1      
 19 24  1  0      
 11 20  1  0      
 21 25  1  0      
 18 26  1  0      
 22 27  1  0      
 23 28  1  0      
 10 12  1  0      
  1  7  1  6      
  1  6  1  1      
 27 29  1  0      
 29 30  2  0      
 29 31  1  0      
 29 32  1  0      
M  END

ENDGPINSPSN1SN2TEMPLATE

  my($GPInsPSn1AlkylSn2TemplateString)=<<ENDGPINSPSN1ALKYLSN2TEMPLATE;
GPInsP sn1 alkyl template structure
  LipdMAPS02060609152D

 29 29  0  0  0  0  0  0  0  0999 V2000
   -3.9172    0.4627    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.6314    0.8737    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.3456    0.4627    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.5045   -0.2516    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -4.3301   -0.2516    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.2028    0.8750    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4888    0.4627    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5403    0.4413    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2922    0.7475    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6559    0.1168    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2922    1.4993    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.0781   -0.6737    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.0781   -1.4993    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.7922   -0.2611    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.2357    0.7653    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.9251    1.1152    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.2443   -0.0695    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5567    0.2841    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.8708   -0.0695    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.5483    1.1152    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.1284    0.9017    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.5527    1.3295    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.7169    1.3619    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.7840    0.2300    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.2368    0.9413    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.9953    0.5545    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    4.9953    1.3795    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.5828   -0.1599    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.7922    0.3410    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  6  1  1  0      
  7  6  1  0      
  9  8  1  0      
  9 10  1  0      
  9 11  2  0      
 12 13  2  0      
 12 14  1  0      
 12  5  1  0      
 15 16  1  0      
 15 20  1  0      
 16 17  1  0      
 19 20  1  0      
 19 18  1  1      
 17 18  1  1      
 16 21  1  0      
  8 17  1  0      
 18 22  1  0      
 15 23  1  0      
 19 24  1  0      
 20 25  1  0      
  7  9  1  0      
  1  5  1  6      
  1  4  1  1      
 24 26  1  0      
 26 27  2  0      
 26 28  1  0      
 26 29  1  0      
M  END

ENDGPINSPSN1ALKYLSN2TEMPLATE

  my($GPInsPSn1TemplateString)=<<ENDGPINSPSN1TEMPLATE;
GPInsP sn1 acyl template structure
  LipdMAPS02060609152D

 29 29  0  0  0  0  0  0  0  0999 V2000
   -3.3997   -0.2611    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.1140    0.1500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.8282   -0.2611    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.5425    0.1500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.5425    0.9755    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.9870   -0.9755    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -3.8125   -0.9755    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -6.2565   -0.2611    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6854    0.1513    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9711   -0.2611    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0225   -0.2826    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7745    0.0236    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1383   -0.6072    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7745    0.7755    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.7538    0.0413    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4429    0.3915    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7621   -0.7933    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0745   -0.4398    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.3888   -0.7933    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.0664    0.3915    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6462    0.1781    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.0705    0.6058    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.2349    0.6381    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.3021   -0.4939    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.7549    0.2177    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.4596   -0.1837    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    5.4596    0.6413    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.0471   -0.8982    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    6.2565   -0.3973    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  2  0      
  4  8  1  0      
  9  1  1  0      
 10  9  1  0      
 12 11  1  0      
 12 13  1  0      
 12 14  2  0      
 15 16  1  0      
 15 20  1  0      
 16 17  1  0      
 19 20  1  0      
 19 18  1  1      
 17 18  1  1      
 16 21  1  0      
 11 17  1  0      
 18 22  1  0      
 15 23  1  0      
 19 24  1  0      
 20 25  1  0      
 10 12  1  0      
  1  7  1  6      
  1  6  1  1      
 24 26  1  0      
 26 27  2  0      
 26 28  1  0      
 26 29  1  0      
M  END

ENDGPINSPSN1TEMPLATE

  my($GPInsPSn1AlkylSn2AlkylTemplateString)=<<ENDGPINSPSN1ALKYLSN2ALKYLTEMPLATE;
GPInsP sn1 alkyl and sn2 alkyl template structure
  LipdMAPS02060609152D

 27 27  0  0  0  0  0  0  0  0999 V2000
   -4.1547    0.0498    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.8689    0.4609    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.5831    0.0498    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.7420   -0.6644    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -4.5676   -0.6644    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.4402    0.4622    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.7262    0.0498    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7776    0.0284    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5297    0.3346    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8933   -0.2960    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5297    1.0865    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.3156   -1.0865    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.9985    0.3524    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6878    0.7024    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0069   -0.4823    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.3195   -0.1287    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.6336   -0.4823    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.3110    0.7024    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1089    0.4889    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.3155    0.9167    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.4797    0.9491    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.5468   -0.1828    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.9997    0.5285    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.7863    0.1493    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    4.7863    0.9743    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.3738   -0.5652    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.5831   -0.0642    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  6  1  1  0      
  7  6  1  0      
  9  8  1  0      
  9 10  1  0      
  9 11  2  0      
 12  5  1  0      
 13 14  1  0      
 13 18  1  0      
 14 15  1  0      
 17 18  1  0      
 17 16  1  1      
 15 16  1  1      
 14 19  1  0      
  8 15  1  0      
 16 20  1  0      
 13 21  1  0      
 17 22  1  0      
 18 23  1  0      
  7  9  1  0      
  1  5  1  6      
  1  4  1  1      
 22 24  1  0      
 24 25  2  0      
 24 26  1  0      
 24 27  1  0      
M  END

ENDGPINSPSN1ALKYLSN2ALKYLTEMPLATE

  my($GPInsPSn1AlkylTemplateString)=<<ENDGPINSPSN1ALKYLTEMPLATE;
GPInsP sn1 alkyl template structure
  LipdMAPS02060609152D

 26 26  0  0  0  0  0  0  0  0999 V2000
   -4.1018   -0.1612    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.8160    0.2499    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.5302   -0.1612    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.6891   -0.8754    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -4.5147   -0.8754    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.3873    0.2512    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6733   -0.1612    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7247   -0.1826    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4768    0.1236    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8404   -0.5071    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4768    0.8754    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.0514    0.1414    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7407    0.4914    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0598   -0.6934    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.3724   -0.3398    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.6865   -0.6934    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.3639    0.4914    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0560    0.2779    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.3684    0.7057    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.5326    0.7381    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.5997   -0.3939    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.0526    0.3175    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.7334   -0.0902    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    4.7334    0.7348    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.3209   -0.8046    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.5302   -0.3037    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  6  1  1  0      
  7  6  1  0      
  9  8  1  0      
  9 10  1  0      
  9 11  2  0      
 12 13  1  0      
 12 17  1  0      
 13 14  1  0      
 16 17  1  0      
 16 15  1  1      
 14 15  1  1      
 13 18  1  0      
  8 14  1  0      
 15 19  1  0      
 12 20  1  0      
 16 21  1  0      
 17 22  1  0      
  7  9  1  0      
  1  5  1  6      
  1  4  1  1      
 21 23  1  0      
 23 24  2  0      
 23 25  1  0      
 23 26  1  0      
M  END

ENDGPINSPSN1ALKYLTEMPLATE

  my($GPInsP3Sn1Sn2TemplateString)=<<ENDGPINSP3SN1SN2TEMPLATE;
GPInsP3 sn1 acyl and sn2 template structure
  LipdMAPS02060609152D

 32 32  0  0  0  0  0  0  0  0999 V2000
   -3.1377    0.3628    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.8521    0.7739    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.5666    0.3628    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.2810    0.7739    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.2810    1.5998    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.7249   -0.3517    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -3.5507   -0.3517    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.9951    0.3628    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4232    0.7753    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7088    0.3628    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.2401    0.3414    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5120    0.6475    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8757    0.0168    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5120    1.3998    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.2989   -0.7738    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.2989   -1.5998    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.0132   -0.3610    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0171    0.6654    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.7059    1.0156    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0249   -0.1695    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.3376    0.1841    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.6522   -0.1695    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.3298    1.0156    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.9090    0.8022    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.3336    1.2299    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.4982    1.2622    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.5657    0.1300    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.0187    0.8417    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.2990   -0.4354    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    5.9951   -0.0334    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.2990   -1.0857    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.4935    0.2906    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  2  0      
  4  8  1  0      
  9  1  1  0      
 10  9  1  0      
 12 11  1  0      
 12 13  1  0      
 12 14  2  0      
 15 16  2  0      
 15 17  1  0      
 15  7  1  0      
 18 19  1  0      
 18 23  1  0      
 19 20  1  0      
 22 23  1  0      
 22 21  1  1      
 20 21  1  1      
 19 24  1  0      
 11 20  1  0      
 21 25  1  0      
 18 26  1  0      
 22 27  1  0      
 23 28  1  0      
 10 12  1  0      
  1  7  1  6      
  1  6  1  1      
 29 30  1  0      
 29 31  2  0      
 29 32  1  0      
 27 29  1  0      
M  END

ENDGPINSP3SN1SN2TEMPLATE

  my($GPInsP3Sn1AlkylSn2TemplateString)=<<ENDGPINSP3SN1ALKYLSN2TEMPLATE;
GPInsP3 sn1 alkyl and sn2 template structure
  LipdMAPS02060609152D

 29 29  0  0  0  0  0  0  0  0999 V2000
   -3.6279    0.4627    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.3421    0.8737    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.0565    0.4627    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.2152   -0.2516    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -4.0408   -0.2516    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.9135    0.8751    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1993    0.4627    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2508    0.4413    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0027    0.7473    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3664    0.1168    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0027    1.4995    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.7888   -0.6737    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.7888   -1.4995    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.5030   -0.2609    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.5256    0.7652    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2147    1.1154    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5338   -0.0695    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.8463    0.2840    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.1606   -0.0695    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.8380    1.1154    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4180    0.9020    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.8423    1.3296    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.0066    1.3619    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.0739    0.2300    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.5268    0.9415    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.8070   -0.3353    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    5.5030    0.0666    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.8070   -0.9855    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.0015    0.3905    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  6  1  1  0      
  7  6  1  0      
  9  8  1  0      
  9 10  1  0      
  9 11  2  0      
 12 13  2  0      
 12 14  1  0      
 12  5  1  0      
 15 16  1  0      
 15 20  1  0      
 16 17  1  0      
 19 20  1  0      
 19 18  1  1      
 17 18  1  1      
 16 21  1  0      
  8 17  1  0      
 18 22  1  0      
 15 23  1  0      
 19 24  1  0      
 20 25  1  0      
  7  9  1  0      
  1  5  1  6      
  1  4  1  1      
 26 27  1  0      
 26 28  2  0      
 26 29  1  0      
 24 26  1  0      
M  END

ENDGPINSP3SN1ALKYLSN2TEMPLATE

  my($GPInsP3Sn1AlkylSn2AlkylTemplateString)=<<ENDGPINSP3SN1ALKYLSN2ALKYLTEMPLATE;
GPInsP3 sn1 alkyl and sn2 template structure
  LipdMAPS02060609152D

 27 27  0  0  0  0  0  0  0  0999 V2000
   -3.8511    0.2057    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.5654    0.6167    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.2797    0.2057    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.4384   -0.5086    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -4.2640   -0.5086    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1368    0.6181    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4225    0.2057    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4740    0.1843    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2260    0.4903    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5896   -0.1402    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2260    1.2425    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.0121   -0.9307    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.3024    0.5082    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.9914    0.8584    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3106   -0.3265    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.6230    0.0270    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.9373   -0.3265    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.6148    0.8584    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.1947    0.6450    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.6190    1.0726    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.7834    1.1049    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.8506   -0.0270    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.3035    0.6845    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.5838   -0.5923    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    5.2797   -0.1904    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.5838   -1.2425    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.7782    0.1335    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  6  1  1  0      
  7  6  1  0      
  9  8  1  0      
  9 10  1  0      
  9 11  2  0      
 12  5  1  0      
 13 14  1  0      
 13 18  1  0      
 14 15  1  0      
 17 18  1  0      
 17 16  1  1      
 15 16  1  1      
 14 19  1  0      
  8 15  1  0      
 16 20  1  0      
 13 21  1  0      
 17 22  1  0      
 18 23  1  0      
  7  9  1  0      
  1  5  1  6      
  1  4  1  1      
 24 25  1  0      
 24 26  2  0      
 24 27  1  0      
 22 24  1  0      
M  END

ENDGPINSP3SN1ALKYLSN2ALKYLTEMPLATE

my($GPInsP3Sn1TemplateString)=<<ENDGPINSP3SN1TEMPLATE;
GPInsP3 sn1 acyl template structure
  LipdMAPS02060609152D

 29 29  0  0  0  0  0  0  0  0999 V2000
   -3.1370    0.1057    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.8513    0.5167    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.5656    0.1057    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.2799    0.5167    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.2799    1.3425    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.7243   -0.6086    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -3.5499   -0.6086    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.9938    0.1057    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4227    0.5181    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7084    0.1057    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.2400    0.0843    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5119    0.3904    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8755   -0.2402    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5119    1.1425    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.0165    0.4083    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.7055    0.7584    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0247   -0.4265    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.3371   -0.0729    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.6514   -0.4265    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.3289    0.7584    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.9088    0.5450    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.3331    0.9726    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.4975    1.0049    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.5647   -0.1270    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.0176    0.5845    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.2979   -0.6923    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    5.9938   -0.2904    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.2979   -1.3425    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.4923    0.0335    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  2  0      
  4  8  1  0      
  9  1  1  0      
 10  9  1  0      
 12 11  1  0      
 12 13  1  0      
 12 14  2  0      
 15 16  1  0      
 15 20  1  0      
 16 17  1  0      
 19 20  1  0      
 19 18  1  1      
 17 18  1  1      
 16 21  1  0      
 11 17  1  0      
 18 22  1  0      
 15 23  1  0      
 19 24  1  0      
 20 25  1  0      
 10 12  1  0      
  1  7  1  6      
  1  6  1  1      
 26 27  1  0      
 26 28  2  0      
 26 29  1  0      
 24 26  1  0      
M  END

ENDGPINSP3SN1TEMPLATE

my($GPInsP3Sn1AlkylTemplateString)=<<ENDGPINSP3SN1ALKYLTEMPLATE;
GPInsP3 sn1 alkyl template structure
  LipdMAPS02060609152D

 26 26  0  0  0  0  0  0  0  0999 V2000
   -3.8520    0.2057    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.5665    0.6168    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.2810    0.2057    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.4392   -0.5087    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -4.2650   -0.5087    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1376    0.6182    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4231    0.2057    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4742    0.1843    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2263    0.4905    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5900   -0.1402    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2263    1.2428    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.3030    0.5084    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.9916    0.8586    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3107   -0.3266    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.6234    0.0271    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.9380   -0.3266    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.6157    0.8586    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.1947    0.6452    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.6194    1.0729    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.7841    1.1052    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.8515   -0.0270    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.3045    0.6847    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.5849   -0.5924    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    5.2810   -0.1904    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.5849   -1.2428    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.7793    0.1335    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  6  1  1  0      
  7  6  1  0      
  9  8  1  0      
  9 10  1  0      
  9 11  2  0      
 12 13  1  0      
 12 17  1  0      
 13 14  1  0      
 16 17  1  0      
 16 15  1  1      
 14 15  1  1      
 13 18  1  0      
  8 14  1  0      
 15 19  1  0      
 12 20  1  0      
 16 21  1  0      
 17 22  1  0      
  7  9  1  0      
  1  5  1  6      
  1  4  1  1      
 23 24  1  0      
 23 25  2  0      
 23 26  1  0      
 21 23  1  0      
M  END

ENDGPINSP3SN1ALKYLTEMPLATE

  my($GPInsP4Sn1Sn2TemplateString)=<<ENDGPINSP4SN1SN2TEMPLATE;
GPInsP4 sn1 acyl and sn2 template structure
  LipdMAPS02060609152D

 32 32  0  0  0  0  0  0  0  0999 V2000
   -3.3775    0.3148    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.0918    0.7259    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.8063    0.3148    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.5207    0.7259    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.5207    1.5518    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.9646   -0.3997    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -3.7904   -0.3997    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -6.2348    0.3148    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6630    0.7272    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9486    0.3148    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.0004    0.2934    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7517    0.5995    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1155   -0.0312    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7517    1.3518    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.5387   -0.8219    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.5387   -1.6478    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.2530   -0.4091    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.7773    0.6174    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4662    0.9676    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7852   -0.2175    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0979    0.1360    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.4124   -0.2175    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.0901    0.9676    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6693    0.7542    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.0939    1.1818    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.2585    1.2142    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.3259    0.0820    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.7789    0.7936    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.6094    1.0881    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    6.2348    0.7270    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.7637    0.5121    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.6094    1.6478    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  2  0      
  4  8  1  0      
  9  1  1  0      
 10  9  1  0      
 12 11  1  0      
 12 13  1  0      
 12 14  2  0      
 15 16  2  0      
 15 17  1  0      
 15  7  1  0      
 18 19  1  0      
 18 23  1  0      
 19 20  1  0      
 22 23  1  0      
 22 21  1  1      
 20 21  1  1      
 19 24  1  0      
 11 20  1  0      
 21 25  1  0      
 18 26  1  0      
 22 27  1  0      
 23 28  1  0      
 10 12  1  0      
  1  7  1  6      
  1  6  1  1      
 29 30  1  0      
 29 31  1  0      
 29 32  2  0      
 28 29  1  0      
M  END

ENDGPINSP4SN1SN2TEMPLATE

  my($GPInsP4Sn1AlkylSn2TemplateString)=<<ENDGPINSP4SN1ALKYLSN2TEMPLATE;
GPInsP4 sn1 alkyl and sn2 template structure
  LipdMAPS02060609152D

 29 29  0  0  0  0  0  0  0  0999 V2000
   -3.8673    0.3147    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.5815    0.7257    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.2958    0.3147    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.4546   -0.3996    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -4.2801   -0.3996    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1530    0.7270    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4388    0.3147    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4904    0.2933    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2423    0.5993    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6060   -0.0312    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2423    1.3514    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.0282   -0.8217    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.0282   -1.6474    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.7423   -0.4090    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.2858    0.6172    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.9750    0.9673    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.2942   -0.2174    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.6066    0.1360    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.9207   -0.2174    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.5982    0.9673    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.1784    0.7540    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.6026    1.1815    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.7668    1.2139    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.8340    0.0820    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.2868    0.7934    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.1171    1.0878    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    5.7423    0.7268    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.2714    0.5120    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.1171    1.6474    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  6  1  1  0      
  7  6  1  0      
  9  8  1  0      
  9 10  1  0      
  9 11  2  0      
 12 13  2  0      
 12 14  1  0      
 12  5  1  0      
 15 16  1  0      
 15 20  1  0      
 16 17  1  0      
 19 20  1  0      
 19 18  1  1      
 17 18  1  1      
 16 21  1  0      
  8 17  1  0      
 18 22  1  0      
 15 23  1  0      
 19 24  1  0      
 20 25  1  0      
  7  9  1  0      
  1  5  1  6      
  1  4  1  1      
 26 27  1  0      
 26 28  1  0      
 26 29  2  0      
 25 26  1  0      
M  END

ENDGPINSP4SN1ALKYLSN2TEMPLATE

  my($GPInsP4Sn1AlkylSn2AlkylTemplateString)=<<ENDGPINSP4SN1ALKYLSN2ALKYLTEMPLATE;
GPInsP4 sn1 alkyl and sn2 alkyl template structure
  LipdMAPS02060609152D

 27 27  0  0  0  0  0  0  0  0999 V2000
   -4.0908   -0.0982    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.8050    0.3129    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.5194   -0.0982    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.6781   -0.8125    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -4.5036   -0.8125    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.3764    0.3142    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6622   -0.0982    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7137   -0.1196    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4656    0.1865    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8294   -0.4441    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4656    0.9386    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.2518   -1.2346    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0627    0.2044    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7518    0.5545    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0710   -0.6303    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.3834   -0.2769    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.6976   -0.6303    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.3751    0.5545    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0449    0.3412    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.3794    0.7687    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.5437    0.8011    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.6110   -0.3309    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.0638    0.3806    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.8941    0.6750    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    5.5194    0.3140    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.0484    0.0992    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.8941    1.2346    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  6  1  1  0      
  7  6  1  0      
  9  8  1  0      
  9 10  1  0      
  9 11  2  0      
 12  5  1  0      
 13 14  1  0      
 13 18  1  0      
 14 15  1  0      
 17 18  1  0      
 17 16  1  1      
 15 16  1  1      
 14 19  1  0      
  8 15  1  0      
 16 20  1  0      
 13 21  1  0      
 17 22  1  0      
 18 23  1  0      
  7  9  1  0      
  1  5  1  6      
  1  4  1  1      
 24 25  1  0      
 24 26  1  0      
 24 27  2  0      
 23 24  1  0      
M  END

ENDGPINSP4SN1ALKYLSN2ALKYLTEMPLATE

  my($GPInsP4Sn1TemplateString)=<<ENDGPINSP4SN1TEMPLATE;
GPInsP4 sn1 acyl template structure
  LipdMAPS02060609152D

 29 29  0  0  0  0  0  0  0  0999 V2000
   -3.3766   -0.3092    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.0907    0.1018    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.8050   -0.3092    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.5192    0.1018    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.5192    0.9275    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.9638   -1.0235    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -3.7894   -1.0235    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -6.2331   -0.3092    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6623    0.1031    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9481   -0.3092    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.0004   -0.3306    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7515   -0.0245    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1152   -0.6551    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7515    0.7276    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.7765   -0.0066    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4658    0.3435    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7850   -0.8413    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0973   -0.4879    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.4115   -0.8413    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.0890    0.3435    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6691    0.1301    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.0933    0.5576    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.2576    0.5900    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.3247   -0.5419    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.7776    0.1695    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.6079    0.4639    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    6.2331    0.1029    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.7621   -0.1119    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.6079    1.0235    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  2  0      
  4  8  1  0      
  9  1  1  0      
 10  9  1  0      
 12 11  1  0      
 12 13  1  0      
 12 14  2  0      
 15 16  1  0      
 15 20  1  0      
 16 17  1  0      
 19 20  1  0      
 19 18  1  1      
 17 18  1  1      
 16 21  1  0      
 11 17  1  0      
 18 22  1  0      
 15 23  1  0      
 19 24  1  0      
 20 25  1  0      
 10 12  1  0      
  1  7  1  6      
  1  6  1  1      
 26 27  1  0      
 26 28  1  0      
 26 29  2  0      
 25 26  1  0      
M  END

ENDGPINSP4SN1TEMPLATE

my($GPInsP4Sn1AlkylTemplateString)=<<ENDGPINSP4SN1ALKYLTEMPLATE;
GPInsP4 sn1 alkyl template structure
  LipdMAPS02060609152D

 26 26  0  0  0  0  0  0  0  0999 V2000
   -4.0917   -0.3093    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.8060    0.1018    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.5205   -0.3093    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.6788   -1.0238    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -4.5046   -1.0238    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.3772    0.1031    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6628   -0.3093    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7138   -0.3307    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4659   -0.0245    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8297   -0.6553    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4659    0.7278    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.0630   -0.0066    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7519    0.3436    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0710   -0.8415    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.3836   -0.4880    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.6982   -0.8415    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.3758    0.3436    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0450    0.1301    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.3796    0.5577    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.5442    0.5902    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.6116   -0.5420    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.0646    0.1695    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.8951    0.4640    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    5.5205    0.1029    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.0494   -0.1119    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.8951    1.0238    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  6  1  1  0      
  7  6  1  0      
  9  8  1  0      
  9 10  1  0      
  9 11  2  0      
 12 13  1  0      
 12 17  1  0      
 13 14  1  0      
 16 17  1  0      
 16 15  1  1      
 14 15  1  1      
 13 18  1  0      
  8 14  1  0      
 15 19  1  0      
 12 20  1  0      
 16 21  1  0      
 17 22  1  0      
  7  9  1  0      
  1  5  1  6      
  1  4  1  1      
 23 24  1  0      
 23 25  1  0      
 23 26  2  0      
 22 23  1  0      
M  END

ENDGPINSP4SN1ALKYLTEMPLATE

  my($GPInsP5Sn1Sn2TemplateString)=<<ENDGPINSP5SN1SN2TEMPLATE;
GPInsP5 sn1 acyl and sn2 template structure
  LipdMAPS02060609152D

 32 32  0  0  0  0  0  0  0  0999 V2000
   -2.6496   -0.0479    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.3639    0.3632    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.0784   -0.0479    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.7928    0.3632    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.7928    1.1890    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2367   -0.7624    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0625   -0.7624    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.5069   -0.0479    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9350    0.3645    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2206   -0.0479    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.7283   -0.0693    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0238    0.2368    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3876   -0.3939    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0238    0.9891    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.8107   -1.1846    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.8107   -2.0106    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.5251   -0.7718    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.5052    0.2546    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.1941    0.6049    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5131   -0.5802    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.8258   -0.2267    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.1404   -0.5802    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.8180    0.6049    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.3972    0.3914    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.8218    0.8191    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.9864    0.8515    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.0538   -0.2808    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.5069    0.4309    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.7619    1.4508    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    5.3873    1.0897    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.9163    0.8747    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.7619    2.0106    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  2  0      
  4  8  1  0      
  9  1  1  0      
 10  9  1  0      
 12 11  1  0      
 12 13  1  0      
 12 14  2  0      
 15 16  2  0      
 15 17  1  0      
 15  7  1  0      
 18 19  1  0      
 18 23  1  0      
 19 20  1  0      
 22 23  1  0      
 22 21  1  1      
 20 21  1  1      
 19 24  1  0      
 11 20  1  0      
 21 25  1  0      
 18 26  1  0      
 22 27  1  0      
 23 28  1  0      
 10 12  1  0      
  1  7  1  6      
  1  6  1  1      
 29 30  1  0      
 29 31  1  0      
 29 32  2  0      
 26 29  1  0      
M  END

ENDGPINSP5SN1SN2TEMPLATE

  my($GPInsP5Sn1AlkylSn2TemplateString)=<<ENDGPINSP5SN1ALKYLSN2TEMPLATE;
GPInsP5 sn1 alkyl and sn2 template structure
  LipdMAPS02060609152D

 29 29  0  0  0  0  0  0  0  0999 V2000
   -3.1396   -0.0479    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.8537    0.3631    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.5681   -0.0479    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.7269   -0.7622    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -3.5524   -0.7622    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4252    0.3644    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7110   -0.0479    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.2373   -0.0693    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5146    0.2367    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8783   -0.3938    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5146    0.9888    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.3004   -1.1843    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.3004   -2.0101    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.0146   -0.7716    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0135    0.2545    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.7027    0.6047    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0219   -0.5800    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.3343   -0.2266    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.6485   -0.5800    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.3259    0.6047    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.9061    0.3913    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.3303    0.8189    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.4945    0.8513    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.5617   -0.2807    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.0146    0.4308    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.2698    1.4504    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    4.8951    1.0894    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.4242    0.8745    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.2698    2.0101    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  6  1  1  0      
  7  6  1  0      
  9  8  1  0      
  9 10  1  0      
  9 11  2  0      
 12 13  2  0      
 12 14  1  0      
 12  5  1  0      
 15 16  1  0      
 15 20  1  0      
 16 17  1  0      
 19 20  1  0      
 19 18  1  1      
 17 18  1  1      
 16 21  1  0      
  8 17  1  0      
 18 22  1  0      
 15 23  1  0      
 19 24  1  0      
 20 25  1  0      
  7  9  1  0      
  1  5  1  6      
  1  4  1  1      
 26 27  1  0      
 26 28  1  0      
 26 29  2  0      
 23 26  1  0      
M  END

ENDGPINSP5SN1ALKYLSN2TEMPLATE


  my($GPInsP5Sn1AlkylSn2AlkylTemplateString)=<<ENDGPINSP5SN1ALKYLSN2ALKYLTEMPLATE;
GPInsP5 sn1 alkyl and sn2 alkyl template structure
  LipdMAPS02060609152D

 27 27  0  0  0  0  0  0  0  0999 V2000
   -3.3630   -0.4608    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.0772   -0.0498    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.7916   -0.4608    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.9503   -1.1752    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -3.7759   -1.1752    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6486   -0.0485    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9344   -0.4608    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.0141   -0.4822    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7379   -0.1762    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1016   -0.8067    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7379    0.5759    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.5239   -1.5973    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.7904   -0.1584    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4795    0.1918    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7987   -0.9930    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.1112   -0.6395    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.4255   -0.9930    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.1029    0.1918    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6829   -0.0216    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.1072    0.4060    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.2714    0.4384    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.3387   -0.6936    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.7916    0.0179    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.0468    1.0376    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    4.6721    0.6765    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.2012    0.4616    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.0468    1.5973    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  6  1  1  0      
  7  6  1  0      
  9  8  1  0      
  9 10  1  0      
  9 11  2  0      
 12  5  1  0      
 13 14  1  0      
 13 18  1  0      
 14 15  1  0      
 17 18  1  0      
 17 16  1  1      
 15 16  1  1      
 14 19  1  0      
  8 15  1  0      
 16 20  1  0      
 13 21  1  0      
 17 22  1  0      
 18 23  1  0      
  7  9  1  0      
  1  5  1  6      
  1  4  1  1      
 24 25  1  0      
 24 26  1  0      
 24 27  2  0      
 21 24  1  0      
M  END

ENDGPINSP5SN1ALKYLSN2ALKYLTEMPLATE

  my($GPInsP5Sn1TemplateString)=<<ENDGPINSP5SN1TEMPLATE;
GPInsP5 sn1 acyl template structure
  LipdMAPS02060609152D

 29 29  0  0  0  0  0  0  0  0999 V2000
   -2.6489   -0.6718    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.3630   -0.2608    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.0773   -0.6718    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.7915   -0.2608    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.7915    0.5647    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2361   -1.3861    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0617   -1.3861    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.5054   -0.6718    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9345   -0.2595    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2203   -0.6718    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.7281   -0.6932    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0238   -0.3872    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3875   -1.0177    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0238    0.3649    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.5042   -0.3694    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.1935   -0.0192    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5127   -1.2040    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.8250   -0.8506    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.1393   -1.2040    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.8167   -0.0192    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.3968   -0.2326    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.8210    0.1949    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.9853    0.2273    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.0524   -0.9047    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.5054   -0.1931    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.7606    0.8265    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    5.3858    0.4655    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.9150    0.2505    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.7606    1.3861    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  2  0      
  4  8  1  0      
  9  1  1  0      
 10  9  1  0      
 12 11  1  0      
 12 13  1  0      
 12 14  2  0      
 15 16  1  0      
 15 20  1  0      
 16 17  1  0      
 19 20  1  0      
 19 18  1  1      
 17 18  1  1      
 16 21  1  0      
 11 17  1  0      
 18 22  1  0      
 15 23  1  0      
 19 24  1  0      
 20 25  1  0      
 10 12  1  0      
  1  7  1  6      
  1  6  1  1      
 26 27  1  0      
 26 28  1  0      
 26 29  2  0      
 23 26  1  0      
M  END

ENDGPINSP5SN1TEMPLATE


  my($GPInsP5Sn1AlkylTemplateString)=<<ENDGPINSP5SN1ALKYLTEMPLATE;
GPInsP5 sn1 alkyl template structure
  LipdMAPS02060609152D

 26 26  0  0  0  0  0  0  0  0999 V2000
   -3.3636   -0.6719    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.0778   -0.2609    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.7923   -0.6719    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.9507   -1.3864    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -3.7765   -1.3864    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6491   -0.2596    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9347   -0.6719    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.0141   -0.6933    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7380   -0.3873    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1018   -1.0179    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7380    0.3650    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.7907   -0.3695    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4797   -0.0192    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7988   -1.2042    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.1114   -0.8508    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.4259   -1.2042    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.1035   -0.0192    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6829   -0.2326    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.1074    0.1949    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.2719    0.2273    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.3392   -0.9049    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.7923   -0.1931    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.0473    0.8267    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    4.6727    0.4656    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.2018    0.2505    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.0473    1.3864    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  6  1  1  0      
  7  6  1  0      
  9  8  1  0      
  9 10  1  0      
  9 11  2  0      
 12 13  1  0      
 12 17  1  0      
 13 14  1  0      
 16 17  1  0      
 16 15  1  1      
 14 15  1  1      
 13 18  1  0      
  8 14  1  0      
 15 19  1  0      
 12 20  1  0      
 16 21  1  0      
 17 22  1  0      
  7  9  1  0      
  1  5  1  6      
  1  4  1  1      
 23 24  1  0      
 23 25  1  0      
 23 26  2  0      
 20 23  1  0      
M  END

ENDGPINSP5SN1ALKYLTEMPLATE

  # GPInsP34, GPInsP35 and GPIns45 templates...

  my($GPInsP34Sn1Sn2TemplateString)=<<ENDGPINSP34SN1SN2TEMPLATE;
GPInsP34 sn1 acyl and sn2 template structure
  LipdMAPS02060609152D

 36 36  0  0  0  0  0  0  0  0999 V2000
   -3.3775    0.3148    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.0918    0.7259    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.8063    0.3148    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.5207    0.7259    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.5207    1.5518    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.9646   -0.3997    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -3.7904   -0.3997    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -6.2348    0.3148    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6630    0.7272    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9486    0.3148    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.0004    0.2934    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7517    0.5995    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1155   -0.0312    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7517    1.3518    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.5387   -0.8219    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.5387   -1.6478    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.2530   -0.4091    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.7773    0.6174    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4662    0.9676    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7852   -0.2175    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0979    0.1360    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.4124   -0.2175    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.0901    0.9676    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6693    0.7542    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.0939    1.1818    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.2585    1.2142    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.3259    0.0820    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.7789    0.7936    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.6094    1.0881    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    6.2348    0.7270    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.7637    0.5121    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.6094    1.6478    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.0593   -0.4834    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    5.7553   -0.0814    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.0593   -1.1337    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.2538    0.2425    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  2  0      
  4  8  1  0      
  9  1  1  0      
 10  9  1  0      
 12 11  1  0      
 12 13  1  0      
 12 14  2  0      
 15 16  2  0      
 15 17  1  0      
 15  7  1  0      
 18 19  1  0      
 18 23  1  0      
 19 20  1  0      
 22 23  1  0      
 22 21  1  1      
 20 21  1  1      
 19 24  1  0      
 11 20  1  0      
 21 25  1  0      
 18 26  1  0      
 22 27  1  0      
 23 28  1  0      
 10 12  1  0      
  1  7  1  6      
  1  6  1  1      
 29 30  1  0      
 29 31  1  0      
 29 32  2  0      
 28 29  1  0      
 33 34  1  0      
 33 35  2  0      
 33 36  1  0      
 27 33  1  0      
M  END

ENDGPINSP34SN1SN2TEMPLATE

  my($GPInsP34Sn1AlkylSn2TemplateString)=<<ENDGPINSP34SN1ALKYLSN2TEMPLATE;
GPInsP34 sn1 alkyl and sn2 template structure
  LipdMAPS02060609152D

 33 33  0  0  0  0  0  0  0  0999 V2000
   -3.8684    0.3148    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.5827    0.7259    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.2972    0.3148    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.4555   -0.3997    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -4.2813   -0.3997    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1539    0.7272    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4395    0.3148    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4905    0.2934    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2426    0.5995    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6064   -0.0312    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2426    1.3518    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.0296   -0.8219    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.0296   -1.6478    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.7439   -0.4091    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.2864    0.6174    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.9753    0.9676    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.2943   -0.2175    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.6070    0.1360    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.9215   -0.2175    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.5992    0.9676    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.1784    0.7542    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.6030    1.1818    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.7676    1.2142    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.8350    0.0820    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.2880    0.7936    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.1185    1.0881    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    5.7439    0.7270    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.2728    0.5121    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.1185    1.6478    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.5684   -0.4834    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    5.2644   -0.0814    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.5684   -1.1337    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.7629    0.2425    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  6  1  1  0      
  7  6  1  0      
  9  8  1  0      
  9 10  1  0      
  9 11  2  0      
 12 13  2  0      
 12 14  1  0      
 12  5  1  0      
 15 16  1  0      
 15 20  1  0      
 16 17  1  0      
 19 20  1  0      
 19 18  1  1      
 17 18  1  1      
 16 21  1  0      
  8 17  1  0      
 18 22  1  0      
 15 23  1  0      
 19 24  1  0      
 20 25  1  0      
  7  9  1  0      
  1  5  1  6      
  1  4  1  1      
 26 27  1  0      
 26 28  1  0      
 26 29  2  0      
 25 26  1  0      
 30 31  1  0      
 30 32  2  0      
 30 33  1  0      
 24 30  1  0      
M  END

ENDGPINSP34SN1ALKYLSN2TEMPLATE

  my($GPInsP34Sn1AlkylSn2AlkylTemplateString)=<<ENDGPINSP34SN1ALKYLSN2ALKYLTEMPLATE;
GPInsP34 sn1 alkyl and sn2 alkyl template structure
  LipdMAPS02060609152D

 31 31  0  0  0  0  0  0  0  0999 V2000
   -4.0917    0.0577    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.8060    0.4688    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.5205    0.0577    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.6788   -0.6567    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -4.5046   -0.6567    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.3772    0.4701    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6628    0.0577    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7138    0.0364    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4659    0.3424    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8297   -0.2882    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4659    1.0947    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.2529   -1.0789    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0630    0.3604    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7519    0.7105    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0710   -0.4745    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.3836   -0.1210    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.6981   -0.4745    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.3758    0.7105    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0449    0.4971    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.3796    0.9247    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.5442    0.9571    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.6116   -0.1750    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.0646    0.5365    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.8951    0.8310    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    5.5205    0.4699    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.0494    0.2550    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.8951    1.3907    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.3450   -0.7404    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    5.0410   -0.3384    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.3450   -1.3907    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.5395   -0.0146    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  6  1  1  0      
  7  6  1  0      
  9  8  1  0      
  9 10  1  0      
  9 11  2  0      
 12  5  1  0      
 13 14  1  0      
 13 18  1  0      
 14 15  1  0      
 17 18  1  0      
 17 16  1  1      
 15 16  1  1      
 14 19  1  0      
  8 15  1  0      
 16 20  1  0      
 13 21  1  0      
 17 22  1  0      
 18 23  1  0      
  7  9  1  0      
  1  5  1  6      
  1  4  1  1      
 24 25  1  0      
 24 26  1  0      
 24 27  2  0      
 23 24  1  0      
 28 29  1  0      
 28 30  2  0      
 28 31  1  0      
 22 28  1  0      
M  END

ENDGPINSP34SN1ALKYLSN2ALKYLTEMPLATE

  my($GPInsP34Sn1TemplateString)=<<ENDGPINSP34SN1TEMPLATE;
GPInsP34 sn1 acyl template structure
  LipdMAPS02060609152D

 33 33  0  0  0  0  0  0  0  0999 V2000
   -3.3775    0.0577    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.0918    0.4688    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.8063    0.0577    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.5207    0.4688    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.5207    1.2947    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.9646   -0.6567    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -3.7904   -0.6567    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -6.2348    0.0577    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6630    0.4701    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9486    0.0577    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.0004    0.0364    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7517    0.3424    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1155   -0.2882    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7517    1.0947    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.7773    0.3603    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4662    0.7105    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7852   -0.4745    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0979   -0.1210    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.4124   -0.4745    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.0901    0.7105    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6693    0.4971    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.0939    0.9247    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.2585    0.9571    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.3259   -0.1750    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.7789    0.5365    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.6094    0.8310    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    6.2348    0.4699    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.7637    0.2550    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.6094    1.3907    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.0593   -0.7404    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    5.7553   -0.3384    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.0593   -1.3907    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.2538   -0.0145    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  2  0      
  4  8  1  0      
  9  1  1  0      
 10  9  1  0      
 12 11  1  0      
 12 13  1  0      
 12 14  2  0      
 15 16  1  0      
 15 20  1  0      
 16 17  1  0      
 19 20  1  0      
 19 18  1  1      
 17 18  1  1      
 16 21  1  0      
 11 17  1  0      
 18 22  1  0      
 15 23  1  0      
 19 24  1  0      
 20 25  1  0      
 10 12  1  0      
  1  7  1  6      
  1  6  1  1      
 26 27  1  0      
 26 28  1  0      
 26 29  2  0      
 25 26  1  0      
 30 31  1  0      
 30 32  2  0      
 30 33  1  0      
 24 30  1  0      
M  END

ENDGPINSP34SN1TEMPLATE

  my($GPInsP34Sn1AlkylTemplateString)=<<ENDGPINSP34SN1ALKYLTEMPLATE;
GPInsP34 sn1 alkyl template structure
  LipdMAPS02060609152D

 30 30  0  0  0  0  0  0  0  0999 V2000
   -4.0919    0.0577    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.8062    0.4688    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.5207    0.0577    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.6790   -0.6567    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -4.5048   -0.6567    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.3774    0.4701    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6629    0.0577    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7139    0.0364    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4660    0.3424    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8298   -0.2882    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4660    1.0947    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.0631    0.3603    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7520    0.7105    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0710   -0.4745    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.3837   -0.1210    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.6982   -0.4745    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.3760    0.7105    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0450    0.4971    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.3797    0.9247    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.5443    0.9571    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.6118   -0.1750    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.0648    0.5365    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.8953    0.8310    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    5.5207    0.4699    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.0496    0.2550    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.8953    1.3907    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.3452   -0.7404    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    5.0412   -0.3384    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.3452   -1.3907    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.5397   -0.0145    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  6  1  1  0      
  7  6  1  0      
  9  8  1  0      
  9 10  1  0      
  9 11  2  0      
 12 13  1  0      
 12 17  1  0      
 13 14  1  0      
 16 17  1  0      
 16 15  1  1      
 14 15  1  1      
 13 18  1  0      
  8 14  1  0      
 15 19  1  0      
 12 20  1  0      
 16 21  1  0      
 17 22  1  0      
  7  9  1  0      
  1  5  1  6      
  1  4  1  1      
 23 24  1  0      
 23 25  1  0      
 23 26  2  0      
 22 23  1  0      
 27 28  1  0      
 27 29  2  0      
 27 30  1  0      
 21 27  1  0      
M  END

ENDGPINSP34SN1ALKYLTEMPLATE

  my($GPInsP35Sn1Sn2TemplateString)=<<ENDGPINSP35SN1SN2TEMPLATE;
GPInsP35 sn1 acyl and sn2 template structure
  LipdMAPS02060609152D

 36 36  0  0  0  0  0  0  0  0999 V2000
   -3.1377   -0.0479    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.8521    0.3632    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.5666   -0.0479    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.2810    0.3632    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.2810    1.1890    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.7249   -0.7624    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -3.5507   -0.7624    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.9951   -0.0479    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4232    0.3645    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7088   -0.0479    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.2401   -0.0693    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5120    0.2368    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8757   -0.3939    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5120    0.9891    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.2989   -1.1846    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.2989   -2.0106    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.0132   -0.7718    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0171    0.2546    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.7059    0.6049    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0249   -0.5802    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.3376   -0.2267    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.6522   -0.5802    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.3298    0.6049    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.9090    0.3914    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.3336    0.8191    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.4982    0.8515    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.5657   -0.2808    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.0187    0.4309    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.2737    1.4508    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    4.8991    1.0897    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.4281    0.8747    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.2737    2.0106    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.2990   -0.8462    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    5.9951   -0.4442    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.2990   -1.4964    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.4935   -0.1202    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  2  0      
  4  8  1  0      
  9  1  1  0      
 10  9  1  0      
 12 11  1  0      
 12 13  1  0      
 12 14  2  0      
 15 16  2  0      
 15 17  1  0      
 15  7  1  0      
 18 19  1  0      
 18 23  1  0      
 19 20  1  0      
 22 23  1  0      
 22 21  1  1      
 20 21  1  1      
 19 24  1  0      
 11 20  1  0      
 21 25  1  0      
 18 26  1  0      
 22 27  1  0      
 23 28  1  0      
 10 12  1  0      
  1  7  1  6      
  1  6  1  1      
 29 30  1  0      
 29 31  1  0      
 29 32  2  0      
 26 29  1  0      
 33 34  1  0      
 33 35  2  0      
 33 36  1  0      
 27 33  1  0      
M  END

ENDGPINSP35SN1SN2TEMPLATE

  my($GPInsP35Sn1AlkylSn2TemplateString)=<<ENDGPINSP35SN1ALKYLSN2TEMPLATE;
GPInsP35 sn1 alkyl and sn2 template structure
  LipdMAPS02060609152D

 33 33  0  0  0  0  0  0  0  0999 V2000
   -3.6286   -0.0479    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.3430    0.3632    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.0575   -0.0479    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.2158   -0.7624    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -4.0416   -0.7624    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.9141    0.3645    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1997   -0.0479    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2508   -0.0693    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0029    0.2368    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3666   -0.3939    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0029    0.9891    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.7898   -1.1846    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.7898   -2.0106    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.5041   -0.7718    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.5261    0.2546    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2149    0.6049    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5339   -0.5802    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.8466   -0.2267    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.1612   -0.5802    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.8388    0.6049    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4180    0.3914    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.8426    0.8191    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.0072    0.8515    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.0747   -0.2808    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.5277    0.4309    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.7827    1.4508    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    4.4081    1.0897    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.9371    0.8747    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.7827    2.0106    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.8080   -0.8462    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    5.5041   -0.4442    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.8080   -1.4964    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.0025   -0.1202    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  6  1  1  0      
  7  6  1  0      
  9  8  1  0      
  9 10  1  0      
  9 11  2  0      
 12 13  2  0      
 12 14  1  0      
 12  5  1  0      
 15 16  1  0      
 15 20  1  0      
 16 17  1  0      
 19 20  1  0      
 19 18  1  1      
 17 18  1  1      
 16 21  1  0      
  8 17  1  0      
 18 22  1  0      
 15 23  1  0      
 19 24  1  0      
 20 25  1  0      
  7  9  1  0      
  1  5  1  6      
  1  4  1  1      
 26 27  1  0      
 26 28  1  0      
 26 29  2  0      
 23 26  1  0      
 30 31  1  0      
 30 32  2  0      
 30 33  1  0      
 24 30  1  0      
M  END

ENDGPINSP35SN1ALKYLSN2TEMPLATE

  my($GPInsP35Sn1AlkylSn2AlkylTemplateString)=<<ENDGPINSP35SN1ALKYLSN2ALKYLTEMPLATE;
GPInsP35 sn1 alkyl and sn2 alkyl template structure
  LipdMAPS02060609152D

 31 31  0  0  0  0  0  0  0  0999 V2000
   -3.8519   -0.3050    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.5663    0.1061    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.2808   -0.3050    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.4391   -1.0195    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -4.2649   -1.0195    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1374    0.1074    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4230   -0.3050    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4741   -0.3264    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2262   -0.0203    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5899   -0.6510    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2262    0.7320    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.0131   -1.4417    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.3028   -0.0025    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.9916    0.3478    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3106   -0.8373    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.6233   -0.4838    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.9379   -0.8373    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.6155    0.3478    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.1947    0.1343    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.6193    0.5620    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.7839    0.5944    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.8514   -0.5379    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.3044    0.1738    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.5594    1.1937    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    4.1848    0.8326    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.7138    0.6176    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.5594    1.7535    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.5847   -1.1033    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    5.2808   -0.7013    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.5847   -1.7535    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.7792   -0.3773    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  6  1  1  0      
  7  6  1  0      
  9  8  1  0      
  9 10  1  0      
  9 11  2  0      
 12  5  1  0      
 13 14  1  0      
 13 18  1  0      
 14 15  1  0      
 17 18  1  0      
 17 16  1  1      
 15 16  1  1      
 14 19  1  0      
  8 15  1  0      
 16 20  1  0      
 13 21  1  0      
 17 22  1  0      
 18 23  1  0      
  7  9  1  0      
  1  5  1  6      
  1  4  1  1      
 24 25  1  0      
 24 26  1  0      
 24 27  2  0      
 21 24  1  0      
 28 29  1  0      
 28 30  2  0      
 28 31  1  0      
 22 28  1  0      
M  END

ENDGPINSP35SN1ALKYLSN2ALKYLTEMPLATE

  my($GPInsP35Sn1TemplateString)=<<ENDGPINSP35SN1TEMPLATE;
GPInsP35 sn1 acyl template structure
  LipdMAPS02060609152D

 33 33  0  0  0  0  0  0  0  0999 V2000
   -3.1377   -0.3050    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.8521    0.1061    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.5666   -0.3050    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.2810    0.1061    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.2810    0.9319    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.7249   -1.0195    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -3.5507   -1.0195    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.9951   -0.3050    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4232    0.1074    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7088   -0.3050    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.2401   -0.3264    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5120   -0.0203    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8757   -0.6510    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5120    0.7320    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.0171   -0.0025    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.7059    0.3478    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0249   -0.8373    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.3376   -0.4838    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.6522   -0.8373    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.3298    0.3478    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.9090    0.1343    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.3336    0.5620    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.4982    0.5944    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.5657   -0.5379    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.0187    0.1738    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.2737    1.1937    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    4.8991    0.8326    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.4281    0.6176    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.2737    1.7535    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.2990   -1.1033    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    5.9951   -0.7013    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.2990   -1.7535    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.4935   -0.3773    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  2  0      
  4  8  1  0      
  9  1  1  0      
 10  9  1  0      
 12 11  1  0      
 12 13  1  0      
 12 14  2  0      
 15 16  1  0      
 15 20  1  0      
 16 17  1  0      
 19 20  1  0      
 19 18  1  1      
 17 18  1  1      
 16 21  1  0      
 11 17  1  0      
 18 22  1  0      
 15 23  1  0      
 19 24  1  0      
 20 25  1  0      
 10 12  1  0      
  1  7  1  6      
  1  6  1  1      
 26 27  1  0      
 26 28  1  0      
 26 29  2  0      
 23 26  1  0      
 30 31  1  0      
 30 32  2  0      
 30 33  1  0      
 24 30  1  0      
M  END

ENDGPINSP35SN1TEMPLATE

  my($GPInsP35Sn1AlkylTemplateString)=<<ENDGPINSP35SN1ALKYLTEMPLATE;
GPInsP35 sn1 alkyl template structure
  LipdMAPS02060609152D

 30 30  0  0  0  0  0  0  0  0999 V2000
   -3.8521   -0.3050    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.5665    0.1061    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.2810   -0.3050    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.4393   -1.0195    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -4.2651   -1.0195    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1375    0.1074    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4231   -0.3050    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4742   -0.3264    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2263   -0.0203    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5900   -0.6510    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2263    0.7320    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.3029   -0.0025    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.9917    0.3478    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3107   -0.8373    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.6234   -0.4838    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.9380   -0.8373    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.6157    0.3478    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.1948    0.1343    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.6194    0.5620    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.7840    0.5944    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.8516   -0.5379    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.3046    0.1738    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.5596    1.1937    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    4.1850    0.8326    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.7140    0.6176    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.5596    1.7536    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.5849   -1.1033    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    5.2810   -0.7013    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.5849   -1.7536    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.7794   -0.3773    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  6  1  1  0      
  7  6  1  0      
  9  8  1  0      
  9 10  1  0      
  9 11  2  0      
 12 13  1  0      
 12 17  1  0      
 13 14  1  0      
 16 17  1  0      
 16 15  1  1      
 14 15  1  1      
 13 18  1  0      
  8 14  1  0      
 15 19  1  0      
 12 20  1  0      
 16 21  1  0      
 17 22  1  0      
  7  9  1  0      
  1  5  1  6      
  1  4  1  1      
 23 24  1  0      
 23 25  1  0      
 23 26  2  0      
 20 23  1  0      
 27 28  1  0      
 27 29  2  0      
 27 30  1  0      
 21 27  1  0      
M  END

ENDGPINSP35SN1ALKYLTEMPLATE

  my($GPInsP45Sn1Sn2TemplateString)=<<ENDGPINSP45SN1SN2TEMPLATE;
GPInsP45 sn1 acyl and sn2 template structure
  LipdMAPS02060609152D

 36 36  0  0  0  0  0  0  0  0999 V2000
   -3.3775   -0.0479    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.0918    0.3632    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.8063   -0.0479    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.5207    0.3632    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.5207    1.1890    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.9646   -0.7624    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -3.7904   -0.7624    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -6.2348   -0.0479    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6630    0.3645    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9486   -0.0479    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.0004   -0.0693    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7517    0.2368    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1155   -0.3939    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7517    0.9891    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.5387   -1.1846    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.5387   -2.0106    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.2530   -0.7718    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.7773    0.2546    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4662    0.6049    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7852   -0.5802    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0979   -0.2267    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.4124   -0.5802    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.0901    0.6049    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6693    0.3914    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.0939    0.8191    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.2585    0.8515    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.3259   -0.2808    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.7789    0.4309    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.0340    1.4508    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    4.6594    1.0897    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.1884    0.8747    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.0340    2.0106    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.6094    0.7254    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    6.2348    0.3643    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.7637    0.1493    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.6094    1.2851    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  2  0      
  4  8  1  0      
  9  1  1  0      
 10  9  1  0      
 12 11  1  0      
 12 13  1  0      
 12 14  2  0      
 15 16  2  0      
 15 17  1  0      
 15  7  1  0      
 18 19  1  0      
 18 23  1  0      
 19 20  1  0      
 22 23  1  0      
 22 21  1  1      
 20 21  1  1      
 19 24  1  0      
 11 20  1  0      
 21 25  1  0      
 18 26  1  0      
 22 27  1  0      
 23 28  1  0      
 10 12  1  0      
  1  7  1  6      
  1  6  1  1      
 29 30  1  0      
 29 31  1  0      
 29 32  2  0      
 26 29  1  0      
 33 34  1  0      
 33 35  1  0      
 33 36  2  0      
 28 33  1  0      
M  END

ENDGPINSP45SN1SN2TEMPLATE

  my($GPInsP45Sn1AlkylSn2TemplateString)=<<ENDGPINSP45SN1ALKYLSN2TEMPLATE;
GPInsP45 sn1 alkyl and sn2 template structure
  LipdMAPS02060609152D

 33 33  0  0  0  0  0  0  0  0999 V2000
   -3.8684   -0.0479    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.5827    0.3632    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.2972   -0.0479    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.4555   -0.7624    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -4.2813   -0.7624    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1539    0.3645    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4395   -0.0479    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4905   -0.0693    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2426    0.2368    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6064   -0.3939    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2426    0.9891    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.0296   -1.1846    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.0296   -2.0106    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.7439   -0.7718    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.2864    0.2546    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.9753    0.6049    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.2943   -0.5802    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.6070   -0.2267    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.9215   -0.5802    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.5992    0.6049    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.1784    0.3914    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.6030    0.8191    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.7676    0.8515    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.8350   -0.2808    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.2880    0.4309    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.5431    1.4508    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    4.1685    1.0897    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.6975    0.8747    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.5431    2.0106    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.1185    0.7254    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    5.7439    0.3643    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.2728    0.1493    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.1185    1.2851    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  6  1  1  0      
  7  6  1  0      
  9  8  1  0      
  9 10  1  0      
  9 11  2  0      
 12 13  2  0      
 12 14  1  0      
 12  5  1  0      
 15 16  1  0      
 15 20  1  0      
 16 17  1  0      
 19 20  1  0      
 19 18  1  1      
 17 18  1  1      
 16 21  1  0      
  8 17  1  0      
 18 22  1  0      
 15 23  1  0      
 19 24  1  0      
 20 25  1  0      
  7  9  1  0      
  1  5  1  6      
  1  4  1  1      
 26 27  1  0      
 26 28  1  0      
 26 29  2  0      
 23 26  1  0      
 30 31  1  0      
 30 32  1  0      
 30 33  2  0      
 25 30  1  0      
M  END

ENDGPINSP45SN1ALKYLSN2TEMPLATE

  my($GPInsP45Sn1AlkylSn2AlkylTemplateString)=<<ENDGPINSP45SN1ALKYLSN2ALKYLTEMPLATE;
GPInsP45 sn1 alkyl and sn2 alkyl template structure
  LipdMAPS02060609152D

 31 31  0  0  0  0  0  0  0  0999 V2000
   -4.0917   -0.4609    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.8060   -0.0498    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.5205   -0.4609    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.6788   -1.1754    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -4.5046   -1.1754    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.3772   -0.0485    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6628   -0.4609    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7138   -0.4823    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4659   -0.1762    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8297   -0.8069    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4659    0.5761    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.2529   -1.5976    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0630   -0.1584    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7519    0.1919    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0710   -0.9932    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.3836   -0.6397    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.6981   -0.9932    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.3758    0.1919    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0450   -0.0216    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.3796    0.4061    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.5442    0.4385    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.6116   -0.6938    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.0646    0.0179    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.3197    1.0378    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    3.9451    0.6767    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.4741    0.4617    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.3197    1.5976    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.8951    0.3124    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    5.5205   -0.0487    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.0494   -0.2637    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.8951    0.8721    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  6  1  1  0      
  7  6  1  0      
  9  8  1  0      
  9 10  1  0      
  9 11  2  0      
 12  5  1  0      
 13 14  1  0      
 13 18  1  0      
 14 15  1  0      
 17 18  1  0      
 17 16  1  1      
 15 16  1  1      
 14 19  1  0      
  8 15  1  0      
 16 20  1  0      
 13 21  1  0      
 17 22  1  0      
 18 23  1  0      
  7  9  1  0      
  1  5  1  6      
  1  4  1  1      
 24 25  1  0      
 24 26  1  0      
 24 27  2  0      
 21 24  1  0      
 28 29  1  0      
 28 30  1  0      
 28 31  2  0      
 23 28  1  0      
M  END

ENDGPINSP45SN1ALKYLSN2ALKYLTEMPLATE

  my($GPInsP45Sn1TemplateString)=<<ENDGPINSP45SN1TEMPLATE;
GPInsP45 sn1 acyl template structure
  LipdMAPS02060609152D

 33 33  0  0  0  0  0  0  0  0999 V2000
   -3.3775   -0.6720    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.0918   -0.2609    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.8063   -0.6720    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.5207   -0.2609    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.5207    0.5649    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.9646   -1.3865    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -3.7904   -1.3865    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -6.2348   -0.6720    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6630   -0.2596    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9486   -0.6720    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.0004   -0.6934    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7517   -0.3873    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1155   -1.0180    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7517    0.3650    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.7773   -0.3695    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4662   -0.0192    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7852   -1.2043    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0979   -0.8508    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.4124   -1.2043    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.0901   -0.0192    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6693   -0.2327    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.0939    0.1950    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.2585    0.2274    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.3259   -0.9049    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.7789   -0.1932    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.0340    0.8267    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    4.6594    0.4656    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.1884    0.2506    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.0340    1.3865    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.6094    0.1013    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    6.2348   -0.2598    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.7637   -0.4748    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.6094    0.6610    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  2  0      
  4  8  1  0      
  9  1  1  0      
 10  9  1  0      
 12 11  1  0      
 12 13  1  0      
 12 14  2  0      
 15 16  1  0      
 15 20  1  0      
 16 17  1  0      
 19 20  1  0      
 19 18  1  1      
 17 18  1  1      
 16 21  1  0      
 11 17  1  0      
 18 22  1  0      
 15 23  1  0      
 19 24  1  0      
 20 25  1  0      
 10 12  1  0      
  1  7  1  6      
  1  6  1  1      
 26 27  1  0      
 26 28  1  0      
 26 29  2  0      
 23 26  1  0      
 30 31  1  0      
 30 32  1  0      
 30 33  2  0      
 25 30  1  0      
M  END

ENDGPINSP45SN1TEMPLATE

  my($GPInsP45Sn1AlkylTemplateString)=<<ENDGPINSP45SN1ALKYLTEMPLATE;
GPInsP45 sn1 alkyl template structure
  LipdMAPS02060609152D

 30 30  0  0  0  0  0  0  0  0999 V2000
   -4.0919   -0.6720    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.8062   -0.2609    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.5207   -0.6720    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.6790   -1.3865    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -4.5048   -1.3865    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.3774   -0.2596    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6629   -0.6720    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7139   -0.6934    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4660   -0.3873    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8298   -1.0180    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4660    0.3650    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.0631   -0.3695    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7520   -0.0192    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0710   -1.2043    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.3837   -0.8508    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.6982   -1.2043    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.3760   -0.0192    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0450   -0.2327    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.3797    0.1950    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.5443    0.2274    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.6118   -0.9049    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.0648   -0.1932    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.3199    0.8267    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    3.9453    0.4656    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.4743    0.2506    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.3199    1.3865    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.8953    0.1013    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    5.5207   -0.2598    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.0496   -0.4748    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.8953    0.6610    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  6  1  1  0      
  7  6  1  0      
  9  8  1  0      
  9 10  1  0      
  9 11  2  0      
 12 13  1  0      
 12 17  1  0      
 13 14  1  0      
 16 17  1  0      
 16 15  1  1      
 14 15  1  1      
 13 18  1  0      
  8 14  1  0      
 15 19  1  0      
 12 20  1  0      
 16 21  1  0      
 17 22  1  0      
  7  9  1  0      
  1  5  1  6      
  1  4  1  1      
 23 24  1  0      
 23 25  1  0      
 23 26  2  0      
 20 23  1  0      
 27 28  1  0      
 27 29  1  0      
 27 30  2  0      
 22 27  1  0      
M  END

ENDGPINSP45SN1ALKYLTEMPLATE

  # GPInsP345 templates...

  my($GPInsP345Sn1Sn2TemplateString)=<<ENDGPINSP345SN1SN2TEMPLATE;
GPInsP345 sn1 acyl and sn2 template structure
  LipdMAPS02060609152D

 40 40  0  0  0  0  0  0  0  0999 V2000
   -3.3609   -0.0477    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.0717    0.3614    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.7827   -0.0477    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.4936    0.3614    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.4936    1.1832    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.9501   -0.7587    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -3.7718   -0.7587    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -6.2042   -0.0477    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6499    0.3627    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9390   -0.0477    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.0004   -0.0690    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7480    0.2356    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1100   -0.3920    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7480    0.9842    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.5164   -1.1788    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.5164   -2.0007    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.2272   -0.7680    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.7637    0.2534    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4590    0.6019    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7813   -0.5774    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0876   -0.2256    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.3957   -0.5774    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.0700    0.6019    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6660    0.3895    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.0836    0.8151    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.2425    0.8473    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.3047   -0.2794    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.7555    0.4288    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.0142    1.4437    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    4.6365    1.0844    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.1678    0.8704    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.0142    2.0007    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.5819    0.7218    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    6.2042    0.3625    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.7354    0.1486    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.5819    1.2788    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.0345   -0.8420    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    5.7271   -0.4420    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.0345   -1.4891    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.2280   -0.1196    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  2  0      
  4  8  1  0      
  9  1  1  0      
 10  9  1  0      
 12 11  1  0      
 12 13  1  0      
 12 14  2  0      
 15 16  2  0      
 15 17  1  0      
 15  7  1  0      
 18 19  1  0      
 18 23  1  0      
 19 20  1  0      
 22 23  1  0      
 22 21  1  1      
 20 21  1  1      
 19 24  1  0      
 11 20  1  0      
 21 25  1  0      
 18 26  1  0      
 22 27  1  0      
 23 28  1  0      
 10 12  1  0      
  1  7  1  6      
  1  6  1  1      
 29 30  1  0      
 29 31  1  0      
 29 32  2  0      
 26 29  1  0      
 33 34  1  0      
 33 35  1  0      
 33 36  2  0      
 28 33  1  0      
 37 38  1  0      
 37 39  2  0      
 37 40  1  0      
 27 37  1  0      
M  END

ENDGPINSP345SN1SN2TEMPLATE

  my($GPInsP345Sn1AlkylSn2TemplateString)=<<ENDGPINSP345SN1ALKYLSN2TEMPLATE;
GPInsP345 sn1 alkyl and sn2 template structure
  LipdMAPS02060609152D

 37 37  0  0  0  0  0  0  0  0999 V2000
   -3.8684   -0.0479    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.5827    0.3632    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.2972   -0.0479    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.4556   -0.7624    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -4.2813   -0.7624    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1539    0.3645    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4395   -0.0479    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4905   -0.0693    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2426    0.2368    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6064   -0.3939    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2426    0.9891    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.0296   -1.1846    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.0296   -2.0106    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.7439   -0.7718    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.2864    0.2546    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.9753    0.6049    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.2942   -0.5802    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.6070   -0.2267    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.9215   -0.5802    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.5992    0.6049    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.1784    0.3914    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.6030    0.8191    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.7676    0.8515    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.8350   -0.2808    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.2880    0.4309    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.5431    1.4508    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    4.1685    1.0897    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.6974    0.8747    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.5431    2.0106    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.1185    0.7254    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    5.7439    0.3643    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.2728    0.1493    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.1185    1.2851    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.5684   -0.8462    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    5.2644   -0.4442    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.5684   -1.4964    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.7629   -0.1202    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  6  1  1  0      
  7  6  1  0      
  9  8  1  0      
  9 10  1  0      
  9 11  2  0      
 12 13  2  0      
 12 14  1  0      
 12  5  1  0      
 15 16  1  0      
 15 20  1  0      
 16 17  1  0      
 19 20  1  0      
 19 18  1  1      
 17 18  1  1      
 16 21  1  0      
  8 17  1  0      
 18 22  1  0      
 15 23  1  0      
 19 24  1  0      
 20 25  1  0      
  7  9  1  0      
  1  5  1  6      
  1  4  1  1      
 26 27  1  0      
 26 28  1  0      
 26 29  2  0      
 23 26  1  0      
 30 31  1  0      
 30 32  1  0      
 30 33  2  0      
 25 30  1  0      
 34 35  1  0      
 34 36  2  0      
 34 37  1  0      
 24 34  1  0      
M  END

ENDGPINSP345SN1ALKYLSN2TEMPLATE

  my($GPInsP345Sn1AlkylSn2AlkylTemplateString)=<<ENDGPINSP345SN1ALKYLSN2ALKYLTEMPLATE;
GPInsP345 sn1 alkyl and sn2 alkyl template structure
  LipdMAPS02060609152D

 35 35  0  0  0  0  0  0  0  0999 V2000
   -4.0917   -0.3050    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.8060    0.1061    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.5205   -0.3050    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.6789   -1.0195    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -4.5046   -1.0195    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.3772    0.1074    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6628   -0.3050    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7138   -0.3264    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4659   -0.0203    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8297   -0.6510    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4659    0.7320    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.2529   -1.4417    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0630   -0.0025    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7519    0.3478    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0708   -0.8373    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.3836   -0.4838    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.6981   -0.8373    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.3758    0.3478    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0450    0.1343    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.3796    0.5620    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.5442    0.5944    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.6116   -0.5379    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.0646    0.1738    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.3197    1.1937    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    3.9451    0.8326    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.4740    0.6176    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.3197    1.7535    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.8951    0.4683    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    5.5205    0.1072    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.0494   -0.1078    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.8951    1.0280    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.3450   -1.1033    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    5.0410   -0.7013    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.3450   -1.7535    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.5395   -0.3773    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  6  1  1  0      
  7  6  1  0      
  9  8  1  0      
  9 10  1  0      
  9 11  2  0      
 12  5  1  0      
 13 14  1  0      
 13 18  1  0      
 14 15  1  0      
 17 18  1  0      
 17 16  1  1      
 15 16  1  1      
 14 19  1  0      
  8 15  1  0      
 16 20  1  0      
 13 21  1  0      
 17 22  1  0      
 18 23  1  0      
  7  9  1  0      
  1  5  1  6      
  1  4  1  1      
 24 25  1  0      
 24 26  1  0      
 24 27  2  0      
 21 24  1  0      
 28 29  1  0      
 28 30  1  0      
 28 31  2  0      
 23 28  1  0      
 32 33  1  0      
 32 34  2  0      
 32 35  1  0      
 22 32  1  0      
M  END

ENDGPINSP345SN1ALKYLSN2ALKYLTEMPLATE

  my($GPInsP345Sn1TemplateString)=<<ENDGPINSP345SN1TEMPLATE;
GPInsP345 sn1 acyl  template structure
  LipdMAPS02060609152D

 37 37  0  0  0  0  0  0  0  0999 V2000
   -3.3775   -0.3050    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.0918    0.1061    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.8063   -0.3050    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.5207    0.1061    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.5207    0.9320    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.9646   -1.0195    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -3.7904   -1.0195    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -6.2348   -0.3050    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6630    0.1074    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9486   -0.3050    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.0004   -0.3264    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7517   -0.0203    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1155   -0.6510    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7517    0.7320    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.7773   -0.0024    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4662    0.3478    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7852   -0.8373    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0979   -0.4838    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.4124   -0.8373    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.0901    0.3478    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6693    0.1344    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.0939    0.5621    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.2585    0.5944    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.3259   -0.5378    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.7789    0.1739    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.0340    1.1938    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    4.6594    0.8327    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.1884    0.6176    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.0340    1.7535    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.6094    0.4683    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    6.2348    0.1072    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.7637   -0.1077    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.6094    1.0280    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.0593   -1.1032    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    5.7553   -0.7012    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.0593   -1.7535    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.2538   -0.3773    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  2  0      
  4  8  1  0      
  9  1  1  0      
 10  9  1  0      
 12 11  1  0      
 12 13  1  0      
 12 14  2  0      
 15 16  1  0      
 15 20  1  0      
 16 17  1  0      
 19 20  1  0      
 19 18  1  1      
 17 18  1  1      
 16 21  1  0      
 11 17  1  0      
 18 22  1  0      
 15 23  1  0      
 19 24  1  0      
 20 25  1  0      
 10 12  1  0      
  1  7  1  6      
  1  6  1  1      
 26 27  1  0      
 26 28  1  0      
 26 29  2  0      
 23 26  1  0      
 30 31  1  0      
 30 32  1  0      
 30 33  2  0      
 25 30  1  0      
 34 35  1  0      
 34 36  2  0      
 34 37  1  0      
 24 34  1  0      
M  END

ENDGPINSP345SN1TEMPLATE

 my($GPInsP345Sn1AlkylTemplateString)=<<ENDGPINSP345SN1ALKYLTEMPLATE;
GPInsP345 sn1 alkyl  template structure
  LipdMAPS02060609152D

 34 34  0  0  0  0  0  0  0  0999 V2000
   -4.0923   -0.3050    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.8067    0.1061    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.5213   -0.3050    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.6793   -1.0196    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -4.5053   -1.0196    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.3777    0.1074    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6632   -0.3050    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7139   -0.3264    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4661   -0.0203    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8300   -0.6511    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4661    0.7321    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.0633   -0.0024    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7521    0.3478    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0710   -0.8374    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.3838   -0.4839    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.6985   -0.8374    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.3763    0.3478    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0450    0.1344    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.3798    0.5622    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.5446    0.5945    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.6121   -0.5379    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.0652    0.1739    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.3202    1.1940    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    3.9457    0.8328    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.4746    0.6177    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.3202    1.7537    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.8958    0.4684    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    5.5213    0.1072    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.0501   -0.1077    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.8958    1.0281    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.3456   -1.1033    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    5.0417   -0.7013    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.3456   -1.7537    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.5402   -0.3774    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  6  1  1  0      
  7  6  1  0      
  9  8  1  0      
  9 10  1  0      
  9 11  2  0      
 12 13  1  0      
 12 17  1  0      
 13 14  1  0      
 16 17  1  0      
 16 15  1  1      
 14 15  1  1      
 13 18  1  0      
  8 14  1  0      
 15 19  1  0      
 12 20  1  0      
 16 21  1  0      
 17 22  1  0      
  7  9  1  0      
  1  5  1  6      
  1  4  1  1      
 23 24  1  0      
 23 25  1  0      
 23 26  2  0      
 20 23  1  0      
 27 28  1  0      
 27 29  1  0      
 27 30  2  0      
 22 27  1  0      
 31 32  1  0      
 31 33  2  0      
 31 34  1  0      
 21 31  1  0      
M  END

ENDGPINSP345SN1ALKYLTEMPLATE

  # GPA templates...

  my($GPASn1Sn2TemplateString)=<<ENDGPASN1SN2TEMPLATE;
GPA sn1 acyl and sn2 acyl template structure
  LipdMAPS02060609152D

 17 16  0  0  0  0  0  0  0  0999 V2000
   -0.1701    0.3628    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8844    0.7740    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5989    0.3628    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.3131    0.7740    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.3131    1.5997    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.2428   -0.3515    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5830   -0.3515    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0274    0.3628    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5444    0.7753    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2588    0.3628    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.0274    0.3456    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.2753    0.6518    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    1.9114    0.0210    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.2753    1.4040    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3313   -0.7737    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3313   -1.5997    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0456   -0.3610    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  2  0      
  4  8  1  0      
  1  7  1  6      
  1  6  1  1      
  9  1  1  0      
 10  9  1  0      
 12 11  1  0      
 12 13  1  0      
 12 14  2  0      
 15 16  2  0      
 15 17  1  0      
 15  7  1  0      
 12 10  1  0      
M  END

ENDGPASN1SN2TEMPLATE

  my($GPASn1AlkylSn2TemplateString)=<<ENDGPASN1ALKYLSN2TEMPLATE;
GPA sn1 alkyl and sn2 acyl template structure
  LipdMAPS02060609152D

 14 13  0  0  0  0  0  0  0  0999 V2000
   -0.6610    0.4607    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3754    0.8719    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0899    0.4607    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2481   -0.2537    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0740   -0.2537    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.0535    0.8732    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7679    0.4607    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.5366    0.4435    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.7845    0.7497    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    1.4206    0.1189    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.7845    1.5019    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8223   -0.6759    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8223   -1.5019    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5366   -0.2632    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  1  5  1  6      
  1  4  1  1      
  6  1  1  0      
  7  6  1  0      
  9  8  1  0      
  9 10  1  0      
  9 11  2  0      
 12 13  2  0      
 12 14  1  0      
 12  5  1  0      
  9  7  1  0      
M  END

ENDGPASN1ALKYLSN2TEMPLATE

  my($GPASn1TemplateString)=<<ENDGPASN1TEMPLATE;
GPA sn1 acyl template structure
  LipdMAPS02060609152D

 14 13  0  0  0  0  0  0  0  0999 V2000
   -0.1701   -0.2613    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8845    0.1499    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5990   -0.2613    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.3132    0.1499    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.3132    0.9757    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.2428   -0.9757    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5830   -0.9757    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0276   -0.2613    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5444    0.1512    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2589   -0.2613    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.0276   -0.2785    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.2754    0.0277    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    1.9115   -0.6031    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.2754    0.7799    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  2  0      
  4  8  1  0      
  1  7  1  6      
  1  6  1  1      
  9  1  1  0      
 10  9  1  0      
 12 11  1  0      
 12 13  1  0      
 12 14  2  0      
 12 10  1  0      
M  END

ENDGPASN1TEMPLATE

  my($GPASn1AlkylSn2AlkylTemplateString)=<<ENDGPASN1ALKYLSN2ALKYLTEMPLATE;
GPA sn1 alkyl and sn2 alkyl template structure
  LipdMAPS02060609152D

 12 11  0  0  0  0  0  0  0  0999 V2000
   -0.8844    0.0477    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5988    0.4589    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.3134    0.0477    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4715   -0.6667    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2974   -0.6667    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1699    0.4602    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5446    0.0477    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.3134    0.0305    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.5612    0.3367    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    1.1973   -0.2941    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.5612    1.0890    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0458   -1.0890    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  1  5  1  6      
  1  4  1  1      
  6  1  1  0      
  7  6  1  0      
  9  8  1  0      
  9 10  1  0      
  9 11  2  0      
 12  5  1  0      
  9  7  1  0      
M  END

ENDGPASN1ALKYLSN2ALKYLTEMPLATE

  my($GPASn1AlkylTemplateString)=<<ENDGPASN1ALKYLTEMPLATE;
GPA sn1 alkyl template structure
  LipdMAPS02060609152D

 11 10  0  0  0  0  0  0  0  0999 V2000
   -0.8844   -0.1634    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5988    0.2478    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.3134   -0.1634    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4715   -0.8779    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2974   -0.8779    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1699    0.2491    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5446   -0.1634    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.3134   -0.1806    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.5612    0.1256    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    1.1973   -0.5052    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.5612    0.8779    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  1  5  1  6      
  1  4  1  1      
  6  1  1  0      
  7  6  1  0      
  9  8  1  0      
  9 10  1  0      
  9 11  2  0      
  9  7  1  0      
M  END

ENDGPASN1ALKYLTEMPLATE

  # GPP templates...

  my($GPPSn1Sn2TemplateString)=<<ENDGPPSN1SN2TEMPLATE;
GPP sn1 acyl and sn2 acyl template structure
  LipdMAPS02060609152D

 21 20  0  0  0  0  0  0  0  0999 V2000
   -0.9670    0.3628    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6813    0.7740    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.3959    0.3628    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1101    0.7740    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1101    1.5998    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5541   -0.3515    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3799   -0.3515    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.8245    0.3628    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2525    0.7753    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4620    0.3628    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.2307    0.3456    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.4785    0.6518    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    1.1146    0.0210    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.4785    1.4041    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1283   -0.7737    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1283   -1.5998    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.8426   -0.3610    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0276    0.5591    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    3.0276    1.3841    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.8245    0.3456    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.6151   -0.1553    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  2  0      
  4  8  1  0      
  1  7  1  6      
  1  6  1  1      
  9  1  1  0      
 10  9  1  0      
 12 11  1  0      
 12 13  1  0      
 12 14  2  0      
 15 16  2  0      
 15 17  1  0      
 15  7  1  0      
 12 10  1  0      
 11 18  1  0      
 18 19  2  0      
 18 20  1  0      
 18 21  1  0      
M  END

ENDGPPSN1SN2TEMPLATE

  my($GPPSn1TemplateString)=<<ENDGPPSN1TEMPLATE;
GPP sn1 acyl template structure
  LipdMAPS02060609152D

 18 17  0  0  0  0  0  0  0  0999 V2000
   -0.9671   -0.2615    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6821    0.1500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.3971   -0.2615    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1118    0.1500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1118    0.9764    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5539   -0.9764    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3803   -0.9764    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.8268   -0.2615    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2521    0.1513    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4630   -0.2615    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.2330   -0.2787    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.4802    0.0277    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    1.1161   -0.6036    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.4802    0.7805    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.0299   -0.0652    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    3.0299    0.7598    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.8268   -0.2787    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.6174   -0.7797    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  2  0      
  4  8  1  0      
  1  7  1  6      
  1  6  1  1      
  9  1  1  0      
 10  9  1  0      
 12 11  1  0      
 12 13  1  0      
 12 14  2  0      
 12 10  1  0      
 11 15  1  0      
 15 16  2  0      
 15 17  1  0      
 15 18  1  0      
M  END

ENDGPPSN1TEMPLATE

  # CDP templates...

  my($CDPSn1Sn2TemplateString)=<<ENDCDGPSN1SN2TEMPLATE;
CDGP sn1 acyl and sn2 acyl template structure
  LipdMAPS02060609152D

 37 38  0  0  0  0  0  0  0  0999 V2000
   -3.1782   -0.3161    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.8926    0.0951    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.6071   -0.3161    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.3214    0.0951    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.3214    0.9208    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.7653   -1.0305    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -3.5912   -1.0305    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -6.0357   -0.3161    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4637    0.0964    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7492   -0.3161    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.0195   -0.3333    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7327   -0.0271    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0966   -0.6579    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7327    0.7251    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.3395   -1.4527    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.3395   -2.2787    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.0538   -1.0400    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.7631   -1.9334    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.7631   -2.4164    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.1464   -2.3574    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.5753   -1.5126    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.3358   -1.5126    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.1464   -1.9334    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.4573   -1.2098    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.5753   -0.0704    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    5.3219    0.3602    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.3219    1.2213    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    4.5753    1.6518    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.8289    1.2213    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.8289    0.3602    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.5753    2.4164    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    6.0357   -0.0523    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.8163   -0.1198    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    0.8163    0.7052    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.4038   -0.8343    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.6132   -0.3333    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.3277   -0.7458    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  2  0      
  4  8  1  0      
  1  7  1  6      
  1  6  1  1      
  9  1  1  0      
 10  9  1  0      
 12 11  1  0      
 12 13  1  0      
 12 14  2  0      
 15 16  2  0      
 15 17  1  0      
 15  7  1  0      
 12 10  1  0      
 18 23  1  1      
 22 18  1  1      
 21 23  1  1      
 18 19  1  0      
 20 23  1  0      
 21 24  1  0      
 22 24  1  0      
 21 25  1  0      
 25 26  1  0      
 26 27  1  0      
 27 28  2  0      
 28 29  1  0      
 29 30  2  0      
 25 30  1  0      
 28 31  1  0      
 26 32  2  0      
 11 33  1  0      
 33 34  2  0      
 33 35  1  0      
 33 36  1  0      
 36 37  1  0      
 37 22  1  0      
M  END

ENDCDGPSN1SN2TEMPLATE

  my($CDPSn1AlkylSn2TemplateString)=<<ENDCDPSN1ALKYLSN2TEMPLATE;
CDP sn1 alkyl and sn2 acyl template structure
  LipdMAPS02060609152D

 34 35  0  0  0  0  0  0  0  0999 V2000
   -3.6674   -0.3007    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.3819    0.1105    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.0964   -0.3007    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.2545   -1.0152    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -4.0805   -1.0152    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.9529    0.1118    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2385   -0.3007    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2218   -0.0117    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5857   -0.6425    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2218    0.7405    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.8288   -1.4374    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.8288   -2.2634    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.5431   -1.0247    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4738   -0.3333    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.2699   -1.9334    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.2699   -2.4164    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.6532   -2.3574    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.0821   -1.5126    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.8426   -1.5126    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.6532   -1.9334    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.9641   -1.2098    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.0821   -0.0704    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    4.8287    0.3602    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.8287    1.2213    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    4.0821    1.6518    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.3356    1.2213    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.3356    0.3602    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.0821    2.4164    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.3231   -0.1198    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    0.3231    0.7052    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0894   -0.8343    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.1200   -0.3333    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.8345   -0.7458    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.5431   -0.0523    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  1  5  1  6      
  1  4  1  1      
  6  1  1  0      
  7  6  1  0      
  8  9  1  0      
  8 10  2  0      
 11 12  2  0      
 11 13  1  0      
 11  5  1  0      
  8  7  1  0      
 15 20  1  1      
 19 15  1  1      
 18 20  1  1      
 15 16  1  0      
 17 20  1  0      
 18 21  1  0      
 19 21  1  0      
 18 22  1  0      
 22 23  1  0      
 23 24  1  0      
 24 25  2  0      
 25 26  1  0      
 26 27  2  0      
 22 27  1  0      
 25 28  1  0      
 14 29  1  0      
 29 30  2  0      
 29 31  1  0      
 29 32  1  0      
 32 33  1  0      
 33 19  1  0      
 23 34  2  0      
  8 14  1  0      
M  END

ENDCDPSN1ALKYLSN2TEMPLATE

  my($CDPSn1TemplateString)=<<ENDCDPSN1TEMPLATE;
CDP sn1 acyl template structure
  LipdMAPS02060609152D

 34 35  0  0  0  0  0  0  0  0999 V2000
   -3.1775   -0.2994    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.8925    0.1121    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.6075   -0.2994    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.3223    0.1121    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.3223    0.9385    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.7643   -1.0143    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -3.5908   -1.0143    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -6.0372   -0.2994    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4625    0.1134    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7475   -0.2994    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7302   -0.0102    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0944   -0.6415    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7302    0.7426    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.0203   -0.3333    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.7639   -1.9334    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.7639   -2.4164    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.1473   -2.3574    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.5762   -1.5126    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.3367   -1.5126    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.1473   -1.9334    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.4582   -1.2098    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.5762   -0.0704    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    5.3227    0.3602    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.3227    1.2213    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    4.5762    1.6518    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.8297    1.2213    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.8297    0.3602    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.5762    2.4164    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.8172   -0.1198    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    0.8172    0.7052    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.4047   -0.8343    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.6140   -0.3333    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.3285   -0.7458    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.0372   -0.0523    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
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
 15 20  1  1      
 19 15  1  1      
 18 20  1  1      
 15 16  1  0      
 17 20  1  0      
 18 21  1  0      
 19 21  1  0      
 18 22  1  0      
 22 23  1  0      
 23 24  1  0      
 24 25  2  0      
 25 26  1  0      
 26 27  2  0      
 22 27  1  0      
 25 28  1  0      
 14 29  1  0      
 29 30  2  0      
 29 31  1  0      
 29 32  1  0      
 32 33  1  0      
 33 19  1  0      
 23 34  2  0      
 11 14  1  0      
M  END

ENDCDPSN1TEMPLATE

  my($CDPSn1AlkylSn2AlkylTemplateString)=<<ENDCDPSN1ALKYLSN2ALKYLTEMPLATE;
CDP sn1 alkyl and sn2 alkyl template structure
  LipdMAPS02060609152D

 32 33  0  0  0  0  0  0  0  0999 V2000
   -3.8871   -0.3001    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.6015    0.1111    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.3160   -0.3001    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.4742   -1.0145    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -4.3001   -1.0145    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1726    0.1124    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4581   -0.3001    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4416   -0.0111    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8054   -0.6419    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4416    0.7412    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.0484   -1.4368    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7009   -0.3333    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.0427   -1.9334    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0427   -2.4164    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.4261   -2.3574    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.8550   -1.5126    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.6155   -1.5126    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.4261   -1.9334    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.7370   -1.2098    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.8550   -0.0704    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    4.6016    0.3602    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.6016    1.2213    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    3.8550    1.6518    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.1085    1.2213    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.1085    0.3602    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.8550    2.4164    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.0960   -0.1198    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    0.0960    0.7052    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3165   -0.8343    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.8929   -0.3333    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.6073   -0.7458    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.3160   -0.0523    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  1  5  1  6      
  1  4  1  1      
  6  1  1  0      
  7  6  1  0      
  8  9  1  0      
  8 10  2  0      
 11  5  1  0      
  8  7  1  0      
 13 18  1  1      
 17 13  1  1      
 16 18  1  1      
 13 14  1  0      
 15 18  1  0      
 16 19  1  0      
 17 19  1  0      
 16 20  1  0      
 20 21  1  0      
 21 22  1  0      
 22 23  2  0      
 23 24  1  0      
 24 25  2  0      
 20 25  1  0      
 23 26  1  0      
 12 27  1  0      
 27 28  2  0      
 27 29  1  0      
 27 30  1  0      
 30 31  1  0      
 31 17  1  0      
 21 32  2  0      
  8 12  1  0      
M  END

ENDCDPSN1ALKYLSN2ALKYLTEMPLATE

  my($CDPSn1AlkylTemplateString)=<<ENDCDPSN1ALKYLTEMPLATE;
CDP sn1 alkyl template structure
  LipdMAPS02060609152D

 31 32  0  0  0  0  0  0  0  0999 V2000
   -3.8768   -0.3044    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.5911    0.1068    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.3057   -0.3044    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.4639   -1.0189    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -4.2898   -1.0189    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1623    0.1081    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4478   -0.3044    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4312   -0.0154    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7951   -0.6462    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4312    0.7368    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7112   -0.3333    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.0324   -1.9334    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0324   -2.4164    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.4158   -2.3574    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.8447   -1.5126    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.6052   -1.5126    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.4158   -1.9334    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.7267   -1.2098    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.8447   -0.0704    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    4.5913    0.3602    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.5913    1.2213    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    3.8447    1.6518    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0982    1.2213    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0982    0.3602    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.8447    2.4164    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.0857   -0.1198    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    0.0857    0.7052    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3268   -0.8343    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.8826   -0.3333    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.5970   -0.7458    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.3057   -0.0523    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  1  5  1  6      
  1  4  1  1      
  6  1  1  0      
  7  6  1  0      
  8  9  1  0      
  8 10  2  0      
  8  7  1  0      
 12 17  1  1      
 16 12  1  1      
 15 17  1  1      
 12 13  1  0      
 14 17  1  0      
 15 18  1  0      
 16 18  1  0      
 15 19  1  0      
 19 20  1  0      
 20 21  1  0      
 21 22  2  0      
 22 23  1  0      
 23 24  2  0      
 19 24  1  0      
 22 25  1  0      
 11 26  1  0      
 26 27  2  0      
 26 28  1  0      
 26 29  1  0      
 29 30  1  0      
 30 16  1  0      
 20 31  2  0      
  8 11  1  0      
M  END

ENDCDPSN1ALKYLTEMPLATE

  # GPnCho templates...

  my($GPnChoSn1Sn2TemplateString)=<<ENDGPNCHOSN1SN2TEMPLATE;
GPnCho sn1 acyl and sn2 acyl template structure
  LipdMAPS03010609152D

 22 21  0  0  0  0  0  0  0  0999 V2000
   -1.2230    0.3628    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9373    0.7740    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6518    0.3628    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.3660    0.7740    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.3660    1.5997    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8101   -0.3515    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6359   -0.3515    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.0803    0.3628    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5085    0.7753    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.2059    0.3628    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.2224    0.6518    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    0.8585    0.0210    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.2224    1.4040    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.3842   -0.7737    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.3842   -1.5997    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0985   -0.3610    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.9369    0.2393    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.6513    0.6518    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.3658    0.2393    0.0000 N   0  3  0  0  0  0  0  0  0  0  0  0
    3.3658   -0.5857    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.0803    0.6518    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.0802   -0.1732    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
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
 17 11  1  0      
 18 17  1  0      
 19 18  1  0      
 19 20  1  0      
 19 21  1  0      
 19 22  1  0      
M  CHG  1  19   1
M  END

ENDGPNCHOSN1SN2TEMPLATE

  my($GPnChoSn1AlkylSn2TemplateString)=<<ENDGPNCHOSN1ALKYLSN2TEMPLATE;
GPnCho sn1 alkyl and sn2 acyl template structure
  LipdMAPS02060609152D

 19 18  0  0  0  0  0  0  0  0999 V2000
   -1.7139    0.4606    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4282    0.8718    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1427    0.4606    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3010   -0.2537    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1268   -0.2537    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9994    0.8731    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2850    0.4606    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.7315    0.7496    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    0.3676    0.1188    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.7315    1.5019    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.8751   -0.6759    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.8751   -1.5019    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.5894   -0.2632    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4460    0.3371    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.1605    0.7496    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.8750    0.3371    0.0000 N   0  3  0  0  0  0  0  0  0  0  0  0
    2.8750   -0.4879    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.5894    0.7496    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.5894   -0.0754    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  1  5  1  6      
  1  4  1  1      
  6  1  1  0      
  7  6  1  0      
  8  9  1  0      
  8 10  2  0      
 11 12  2  0      
 11 13  1  0      
 11  5  1  0      
  8  7  1  0      
 14  8  1  0      
 15 14  1  0      
 16 15  1  0      
 16 17  1  0      
 16 18  1  0      
 16 19  1  0      
M  CHG  1  16   1
M  END

ENDGPNCHOSN1ALKYLSN2TEMPLATE

  my($GPnChoSn1TemplateString)=<<ENDGPNCHOSN1TEMPLATE;
GPnCho sn1 acyl template structure
  LipdMAPS02060609152D

 19 18  0  0  0  0  0  0  0  0999 V2000
   -1.2229   -0.1443    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9373    0.2669    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6519   -0.1443    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.3661    0.2669    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.3661    1.0927    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8100   -0.8586    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6359   -0.8586    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.0805   -0.1443    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5084    0.2682    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.2061   -0.1443    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.2226    0.1447    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    0.8587   -0.4861    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.2226    0.8970    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.9371   -0.2678    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.6516    0.1447    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.3660   -0.2678    0.0000 N   0  3  0  0  0  0  0  0  0  0  0  0
    4.0805    0.1447    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.3660   -1.0927    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.0805   -0.6802    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
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
 14 11  1  0      
 15 14  1  0      
 16 15  1  0      
 16 17  1  0      
 16 18  1  0      
 16 19  1  0      
M  CHG  1  16   1
M  END

ENDGPNCHOSN1TEMPLATE

  my($GPnChoSn1AlkylSn2AlkylTemplateString)=<<ENDGPNCHOSN1ALKYLSN2ALKYLTEMPLATE;
GPnCho sn1 alkyl and sn2 alkyl template structure
  LipdMAPS02060609152D

 17 16  0  0  0  0  0  0  0  0999 V2000
   -2.0487    0.0476    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.7630    0.4588    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.4775    0.0476    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6358   -0.6667    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4616   -0.6667    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3342    0.4601    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6198    0.0476    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.3968    0.3366    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    0.0328   -0.2942    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.3968    1.0889    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.2099   -1.0889    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.3342   -0.0759    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0487    0.3366    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.7631   -0.0759    0.0000 N   0  3  0  0  0  0  0  0  0  0  0  0
    2.7631   -0.9008    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.4775    0.3366    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.4775   -0.4884    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  1  5  1  6      
  1  4  1  1      
  6  1  1  0      
  7  6  1  0      
  8  9  1  0      
  8 10  2  0      
 11  5  1  0      
  8  7  1  0      
 12  8  1  0      
 13 12  1  0      
 14 13  1  0      
 14 15  1  0      
 14 16  1  0      
 14 17  1  0      
M  CHG  1  14   1
M  END

ENDGPNCHOSN1ALKYLSN2ALKYLTEMPLATE

  my($GPnChoSn1AlkylTemplateString)=<<ENDGPNCHOSN1ALKYLTEMPLATE;
GPnCho sn1 alkyl template structure
  LipdMAPS02060609152D

 16 15  0  0  0  0  0  0  0  0999 V2000
   -2.0487   -0.0464    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.7630    0.3648    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.4775   -0.0464    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6358   -0.7608    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4616   -0.7608    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3342    0.3661    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6198   -0.0464    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.3968    0.2426    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    0.0328   -0.3882    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.3968    0.9949    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.3342   -0.1699    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0487    0.2426    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.7631   -0.1699    0.0000 N   0  3  0  0  0  0  0  0  0  0  0  0
    3.4775    0.2425    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.7631   -0.9949    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.4775   -0.5824    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  1  5  1  6      
  1  4  1  1      
  6  1  1  0      
  7  6  1  0      
  8  9  1  0      
  8 10  2  0      
  8  7  1  0      
 11  8  1  0      
 12 11  1  0      
 13 12  1  0      
 13 14  1  0      
 13 15  1  0      
 13 16  1  0      
M  CHG  1  13   1
M  END

ENDGPNCHOSN1ALKYLTEMPLATE

  # GPnEtn templates...

  my($GPnEtnSn1Sn2TemplateString)=<<ENDGPNETNSN1SN2TEMPLATE;
GPnEtn sn1 acyl and sn2 acyl template structure
  LipdMAPS02060609152D

 19 18  0  0  0  0  0  0  0  0999 V2000
   -0.8658    0.3628    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5801    0.7741    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2947    0.3628    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0090    0.7741    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0090    1.5998    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4529   -0.3515    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2787   -0.3515    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.7232    0.3628    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1512    0.7754    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5632    0.3628    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.5798    0.6518    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    1.2159    0.0210    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.5798    1.4041    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0271   -0.7738    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0271   -1.5998    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.7414   -0.3610    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.2943    0.2393    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0087    0.6518    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.7232    0.2393    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
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
 17 11  1  0      
 18 17  1  0      
 19 18  1  0      
M  END

ENDGPNETNSN1SN2TEMPLATE

  my($GPnEtnSn1AlkylSn2TemplateString)=<<ENDGPNETNSN1ALKYLSN2TEMPLATE;
GPnEtn sn1 alkyl and sn2 acyl template structure
  LipdMAPS02060609152D

 16 15  0  0  0  0  0  0  0  0999 V2000
   -1.3567    0.4607    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0711    0.8720    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.7857    0.4607    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9437   -0.2537    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7696   -0.2537    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6421    0.8733    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0724    0.4607    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.0890    0.7497    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    0.7251    0.1189    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.0890    1.5021    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5181   -0.6759    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5181   -1.5021    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.2324   -0.2632    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.8035    0.3372    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.5180    0.7497    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.2324    0.3372    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  1  5  1  6      
  1  4  1  1      
  6  1  1  0      
  7  6  1  0      
  8  9  1  0      
  8 10  2  0      
 11 12  2  0      
 11 13  1  0      
 11  5  1  0      
  8  7  1  0      
 14  8  1  0      
 15 14  1  0      
 16 15  1  0      
M  END

ENDGPNETNSN1ALKYLSN2TEMPLATE

  my($GPnEtnSn1TemplateString)=<<ENDGPNETNSN1TEMPLATE;
GPnEtn sn1 acyl template structure
  LipdMAPS02060609152D

 16 15  0  0  0  0  0  0  0  0999 V2000
   -0.8658   -0.2613    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5801    0.1499    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2947   -0.2613    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0090    0.1499    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0090    0.9757    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4529   -0.9757    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2787   -0.9757    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.7232   -0.2613    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1512    0.1512    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5632   -0.2613    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.5798    0.0277    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    1.2159   -0.6031    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.5798    0.7800    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.2943   -0.3848    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0087    0.0277    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.7232   -0.3848    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
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
 14 11  1  0      
 15 14  1  0      
 16 15  1  0      
M  END

ENDGPNETNSN1TEMPLATE

  my($GPnEtnSn1AlkylSn2AlkylTemplateString)=<<ENDGPNETNSN1ALKYLSN2ALKYLTEMPLATE;
GPnEtn sn1 alkyl and sn2 alkyl template structure
  LipdMAPS02060609152D

 15 14  0  0  0  0  0  0  0  0999 V2000
   -1.9563    0.0477    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6706    0.4589    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.3852    0.0477    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5432   -0.6668    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.3692   -0.6668    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2416    0.4602    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5272    0.0477    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.2417    0.0305    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.9562   -0.3821    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.6708    0.0305    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.3852   -0.3821    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.4895    0.3367    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    0.1256   -0.2942    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.4895    1.0890    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1176   -1.0890    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  1  5  1  6      
  1  4  1  1      
  6  1  1  0      
  7  6  1  0      
  9  8  1  0      
 10  9  1  0      
 11 10  1  0      
 12  8  1  0      
 12 13  1  0      
 12 14  2  0      
 15  5  1  0      
 12  7  1  0      
M  END

ENDGPNETNSN1ALKYLSN2ALKYLTEMPLATE

  my($GPnEtnSn1AlkylTemplateString)=<<ENDGPNETNSN1ALKYLTEMPLATE;
GPnEtn sn1 alkyl template structure
  LipdMAPS02060609152D

 17 15  0  0  0  0  0  0  0  0999 V2000
   -3.8949    0.3438    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.6092    0.7551    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.3238    0.3438    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.4818   -0.3706    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -4.3078   -0.3706    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1802    0.7564    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4658    0.3438    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4491    0.6329    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8130    0.0020    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4491    1.3852    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7346    0.2204    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0201    0.6329    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6943    0.2204    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    3.1803   -0.9727    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.8948   -1.3852    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.6094   -0.9727    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.3238   -1.3852    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  1  5  1  6      
  1  4  1  1      
  6  1  1  0      
  7  6  1  0      
  8  9  1  0      
  8 10  2  0      
  8  7  1  0      
 11  8  1  0      
 12 11  1  0      
 13 12  1  0      
 15 14  1  0      
 16 15  1  0      
 17 16  1  0      
M  END

ENDGPNETNSN1ALKYLTEMPLATE

  # GPnEtnNMe2 templates...

  my($GPEtnNMe2Sn1Sn2TemplateString)=<<ENDGPNETNNME2SN1SN2TEMPLATE;
GPEtnNMe2 sn1 acyl and sn2 acyl template structure
  LipdMAPS02060609152D

 22 21  0  0  0  0  0  0  0  0999 V2000
   -1.5991    0.3628    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.3135    0.7741    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0280    0.3628    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.7423    0.7741    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.7423    1.5998    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1862   -0.3515    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0120   -0.3515    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.4565    0.3628    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8846    0.7754    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1701    0.3628    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.5986    0.3456    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.3131   -0.0669    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0276    0.3456    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.7421   -0.0669    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.8465    0.6518    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    0.4825    0.0210    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.8465    1.4041    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.7604   -0.7738    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.7604   -1.5998    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.4747   -0.3610    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.4565    0.3455    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.7421   -0.8918    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  2  0      
  4  8  1  0      
  1  7  1  6      
  1  6  1  1      
  9  1  1  0      
 10  9  1  0      
 12 11  1  0      
 13 12  1  0      
 14 13  1  0      
 15 11  1  0      
 15 16  1  0      
 15 17  2  0      
 18 19  2  0      
 18 20  1  0      
 18  7  1  0      
 15 10  1  0      
 14 21  1  0      
 14 22  1  0      
M  END

ENDGPNETNNME2SN1SN2TEMPLATE

  my($GPEtnNMe2Sn1AlkylSn2TemplateString)=<<ENDGPNETNNME2SN1ALKYLSN2TEMPLATE;
GPEtnNMe2 sn1 alkyl and sn2 acyl template structure
  LipdMAPS02060609152D

 19 18  0  0  0  0  0  0  0  0999 V2000
   -2.0901    0.4607    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.8045    0.8720    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.5191    0.4607    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6771   -0.2537    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5030   -0.2537    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3755    0.8733    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6610    0.4607    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.1079    0.4435    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.8223    0.0310    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.5370    0.4435    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.2514    0.0310    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.3556    0.7497    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0083    0.1189    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.3556    1.5021    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.2515   -0.6759    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.2515   -1.5021    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.9658   -0.2632    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.9658    0.4434    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.2514   -0.7939    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  1  5  1  6      
  1  4  1  1      
  6  1  1  0      
  7  6  1  0      
  9  8  1  0      
 10  9  1  0      
 11 10  1  0      
 12  8  1  0      
 12 13  1  0      
 12 14  2  0      
 15 16  2  0      
 15 17  1  0      
 15  5  1  0      
 12  7  1  0      
 11 18  1  0      
 11 19  1  0      
M  END

ENDGPNETNNME2SN1ALKYLSN2TEMPLATE

  my($GPEtnNMe2Sn1TemplateString)=<<ENDGPNETNNME2SN1TEMPLATE;
GPEtnNMe2 sn1 acyl template structure
  LipdMAPS02060609152D

 19 18  0  0  0  0  0  0  0  0999 V2000
   -1.5991    0.0088    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.3135    0.4201    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0280    0.0088    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.7423    0.4201    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.7423    1.2458    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1862   -0.7055    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0120   -0.7055    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.4565    0.0088    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8846    0.4214    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1701    0.0088    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.5986   -0.0084    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.3131   -0.4209    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0276   -0.0084    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.7421   -0.4209    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.8465    0.2979    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    0.4825   -0.3330    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.8465    1.0501    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.4565   -0.0085    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.7421   -1.2458    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  2  0      
  4  8  1  0      
  1  7  1  6      
  1  6  1  1      
  9  1  1  0      
 10  9  1  0      
 12 11  1  0      
 13 12  1  0      
 14 13  1  0      
 15 11  1  0      
 15 16  1  0      
 15 17  2  0      
 15 10  1  0      
 14 18  1  0      
 14 19  1  0      
M  END

ENDGPNETNNME2SN1TEMPLATE

  my($GPEtnNMe2Sn1AlkylSn2AlkylTemplateString)=<<ENDGPNETNNME2SN1ALKYLSN2ALKYLTEMPLATE;
GPEtnNMe2 sn1 alkyl and sn2 alkyl template structure
  LipdMAPS02060609152D

 17 16  0  0  0  0  0  0  0  0999 V2000
   -2.3134    0.1067    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0276    0.5179    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.7421    0.1067    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9003   -0.6077    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.7262   -0.6077    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5987    0.5192    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8844    0.1067    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.8843    0.0895    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.5988   -0.3230    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.3133    0.0895    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0277   -0.3230    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.1322    0.3957    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2316   -0.2351    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.1322    1.1479    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.4746   -1.0299    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.7421    0.0893    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0275   -1.1479    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  1  5  1  6      
  1  4  1  1      
  6  1  1  0      
  7  6  1  0      
  9  8  1  0      
 10  9  1  0      
 11 10  1  0      
 12  8  1  0      
 12 13  1  0      
 12 14  2  0      
 15  5  1  0      
 12  7  1  0      
 11 16  1  0      
 11 17  1  0      
M  END

ENDGPNETNNME2SN1ALKYLSN2ALKYLTEMPLATE

  my($GPEtnNMe2Sn1AlkylTemplateString)=<<ENDGPNETNNME2SN1ALKYLTEMPLATE;
GPEtnNMe2 sn1 alkyl template structure
  LipdMAPS02060609152D

 16 15  0  0  0  0  0  0  0  0999 V2000
   -2.3135    0.1066    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0279    0.5179    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.7425    0.1066    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9004   -0.6078    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.7264   -0.6078    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5988    0.5192    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8844    0.1066    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.8845    0.0894    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.5990   -0.3231    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.3136    0.0894    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0280   -0.3231    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.1323    0.3957    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2316   -0.2352    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.1323    1.1480    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.7425    0.0893    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0280   -1.1480    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  1  5  1  6      
  1  4  1  1      
  6  1  1  0      
  7  6  1  0      
  9  8  1  0      
 10  9  1  0      
 11 10  1  0      
 12  8  1  0      
 12 13  1  0      
 12 14  2  0      
 12  7  1  0      
 11 15  1  0      
 11 16  1  0      
M  END

ENDGPNETNNME2SN1ALKYLTEMPLATE

  # GPnEtnNMe templates...

  my($GPEtnNMeSn1Sn2TemplateString)=<<ENDGPNETNNMESN1SN2TEMPLATE;
GPEtnNMe sn1 acyl and sn2 acyl template structure
  LipdMAPS02060609152D

 21 20  0  0  0  0  0  0  0  0999 V2000
   -1.5991    0.3628    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.3135    0.7741    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0280    0.3628    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.7423    0.7741    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.7423    1.5998    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1862   -0.3515    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0120   -0.3515    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.4565    0.3628    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8846    0.7754    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1701    0.3628    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.5986    0.3456    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.3131   -0.0669    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0276    0.3456    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.7421   -0.0669    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.8465    0.6518    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    0.4825    0.0210    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.8465    1.4041    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.7604   -0.7738    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.7604   -1.5998    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.4747   -0.3610    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.4565    0.3455    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  2  0      
  4  8  1  0      
  1  7  1  6      
  1  6  1  1      
  9  1  1  0      
 10  9  1  0      
 12 11  1  0      
 13 12  1  0      
 14 13  1  0      
 15 11  1  0      
 15 16  1  0      
 15 17  2  0      
 18 19  2  0      
 18 20  1  0      
 18  7  1  0      
 15 10  1  0      
 14 21  1  0      
M  END

ENDGPNETNNMESN1SN2TEMPLATE

  my($GPEtnNMeSn1AlkylSn2TemplateString)=<<ENDGPNETNNMESN1ALKYLSN2TEMPLATE;
GPEtnNMe sn1 alkyl and sn2 acyl template structure
  LipdMAPS02060609152D

 18 17  0  0  0  0  0  0  0  0999 V2000
   -2.0901    0.4607    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.8045    0.8720    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.5191    0.4607    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6771   -0.2537    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5030   -0.2537    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3755    0.8733    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6610    0.4607    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.1079    0.4435    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.8223    0.0310    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.5370    0.4435    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.2514    0.0310    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.3556    0.7497    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0083    0.1189    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.3556    1.5021    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.2515   -0.6759    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.2515   -1.5021    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.9658   -0.2632    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.9658    0.4434    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  1  5  1  6      
  1  4  1  1      
  6  1  1  0      
  7  6  1  0      
  9  8  1  0      
 10  9  1  0      
 11 10  1  0      
 12  8  1  0      
 12 13  1  0      
 12 14  2  0      
 15 16  2  0      
 15 17  1  0      
 15  5  1  0      
 12  7  1  0      
 11 18  1  0      
M  END

ENDGPNETNNMESN1ALKYLSN2TEMPLATE

  my($GPEtnNMeSn1TemplateString)=<<ENDGPNETNNMESN1TEMPLATE;
GPEtnNMe sn1 acyl template structure
  LipdMAPS02060609152D

 18 17  0  0  0  0  0  0  0  0999 V2000
   -1.5991   -0.2613    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.3135    0.1499    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0280   -0.2613    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.7423    0.1499    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.7423    0.9757    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1862   -0.9757    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0120   -0.9757    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.4565   -0.2613    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8846    0.1512    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1701   -0.2613    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.5986   -0.2785    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.3131   -0.6911    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0276   -0.2785    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.7421   -0.6911    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.8465    0.0277    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    0.4825   -0.6031    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.8465    0.7800    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.4565   -0.2786    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  2  0      
  4  8  1  0      
  1  7  1  6      
  1  6  1  1      
  9  1  1  0      
 10  9  1  0      
 12 11  1  0      
 13 12  1  0      
 14 13  1  0      
 15 11  1  0      
 15 16  1  0      
 15 17  2  0      
 15 10  1  0      
 14 18  1  0      
M  END

ENDGPNETNNMESN1TEMPLATE

  my($GPEtnNMeSn1AlkylSn2AlkylTemplateString)=<<ENDGPNETNNMESN1ALKYLSN2ALKYLTEMPLATE;
GPEtnNMe sn1 alkyl and sn2 alkyl template structure
  LipdMAPS02060609152D

 16 15  0  0  0  0  0  0  0  0999 V2000
   -2.3134    0.0477    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0276    0.4589    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.7421    0.0477    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9003   -0.6667    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.7262   -0.6667    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5987    0.4602    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8844    0.0477    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.8843    0.0305    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.5988   -0.3821    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.3133    0.0305    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0277   -0.3821    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.1322    0.3367    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2316   -0.2942    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.1322    1.0889    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.4746   -1.0889    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.7421    0.0303    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  1  5  1  6      
  1  4  1  1      
  6  1  1  0      
  7  6  1  0      
  9  8  1  0      
 10  9  1  0      
 11 10  1  0      
 12  8  1  0      
 12 13  1  0      
 12 14  2  0      
 15  5  1  0      
 12  7  1  0      
 11 16  1  0      
M  END

ENDGPNETNNMESN1ALKYLSN2ALKYLTEMPLATE

  my($GPEtnNMeSn1AlkylTemplateString)=<<ENDGPNETNNMESN1ALKYLTEMPLATE;
GPEtnNMe sn1 alkyl template structure
  LipdMAPS02060609152D

 15 14  0  0  0  0  0  0  0  0999 V2000
   -2.3135   -0.1635    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0279    0.2478    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.7425   -0.1635    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9004   -0.8779    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.7264   -0.8779    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5988    0.2491    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8844   -0.1635    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.8845   -0.1807    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.5990   -0.5932    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.3136   -0.1807    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0280   -0.5932    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.1323    0.1256    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2316   -0.5053    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.1323    0.8779    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.7425   -0.1808    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  1  5  1  6      
  1  4  1  1      
  6  1  1  0      
  7  6  1  0      
  9  8  1  0      
 10  9  1  0      
 11 10  1  0      
 12  8  1  0      
 12 13  1  0      
 12 14  2  0      
 12  7  1  0      
 11 15  1  0      
M  END

ENDGPNETNNMESN1ALKYLTEMPLATE


# Format: ID => AbbrevID|HeadGroupName|Y1Sn1|Y2Sn1|Y1Sn2|Y2Sn2|Y1Sn3|Y2Sn3|Sn1ChainAtomNum|Sn2ChainAtomNum|Sn3ChainAtomNum|Sn1ChainCarbons|Sn2ChainCarbons|Sn3ChainCarbons|Sn2CAtomNum|Sn2OAtomNum|Sn2HAtomNum|LMCategory|LMMainClass|LMSubClass|TemplateString
  # Note: Template structure for Both Alkyl and Alkenyl substituted GP is the same; Different entries are used to facilitate assigning LM identity
  %GPTemplatesDataMap = (
			       "PCSn1Sn2" => "PC|glycero-3-phosphocholine|0.3628|0.7753|-0.7737|-0.3611|0|0|8|23|0|2|2|0|1|7|6|GP|01|01|$GPChoSn1Sn2TemplateString",
			       "PCSn1AlkylSn2" => "PC|glycero-3-phosphocholine|0.4606|0.8718|-0.6759|-0.2632|0|0|3|20|0|0|2|0|1|5|4|GP|01|02|$GPChoSn1AlkylSn2TemplateString",
			       "PCSn1AlkenylSn2" => "PC|glycero-3-phosphocholine|0.4606|0.8718|-0.6759|-0.2632|0|0|3|20|0|0|2|0|1|5|4|GP|01|03|$GPChoSn1AlkylSn2TemplateString",
			       "PCSn1AlkylSn2Alkyl" => "PC|glycero-3-phosphocholine|0.1067|0.5179|-1.0298|-0.6076|0|0|3|18|0|0|1|0|1|5|4|GP|01|04|$GPChoSn1AlkylSn2AlkylTemplateString",
			       "PCSn1" => "PC|glycero-3-phosphocholine|0.0089|0.4214|0|0|0|0|8|0|0|2|0|0|1|7|6|GP|01|05|$GPChoSn1TemplateString",
			       "PCSn1Alkyl" => "PC|glycero-3-phosphocholine|0.1067|0.5179|0|0|0|0|3|0|0|0|0|0|1|5|4|GP|01|06|$GPChoSn1AlkylTemplateString",
			       "PCSn1Alkenyl" => "PC|glycero-3-phosphocholine|0.1067|0.5179|0|0|0|0|3|0|0|0|0|0|1|5|4|GP|01|07|$GPChoSn1AlkylTemplateString",


			       "PESn1Sn2" => "PE|glycero-3-phosphoethanolamine|0.3628|0.7753|-0.7737|-0.3611|0|0|8|20|0|2|2|0|1|7|6|GP|02|01|$GPEtnSn1Sn2TemplateString",
			       "PESn1AlkylSn2" => "PE|glycero-3-phosphoethanolamine|0.4607|0.8719|-0.6759|-0.2632|0|0|3|17|0|0|2|0|1|5|4|GP|02|02|$GPEtnSn1AlkylSn2TemplateString",
			       "PESn1AlkenylSn2" => "PE|glycero-3-phosphoethanolamine|0.4607|0.8719|-0.6759|-0.2632|0|0|3|17|0|0|2|0|1|5|4|GP|02|03|$GPEtnSn1AlkylSn2TemplateString",
			       "PESn1AlkylSn2Alkyl" => "PE|glycero-3-phosphoethanolamine|0.0477|0.4588|-1.0890|-0.6668|0|0|3|15|0|0|1|0|1|5|4|GP|02|04|$GPEtnSn1AlkylSn2AlkylTemplateString",
			       "PESn1" => "PE|glycero-3-phosphoethanolamine|-0.2613|0.1512|0|0|0|0|8|0|0|2|0|0|1|7|6|GP|02|05|$GPEtnSn1TemplateString",
			       "PESn1Alkyl" => "PE|glycero-3-phosphoethanolamine|-0.1635|0.2478|0|0|0|0|3|0|0|0|0|0|1|5|4|GP|02|06|$GPEtnSn1AlkylTemplateString",
			       "PESn1Alkenyl" => "PE|glycero-3-phosphoethanolamine|-0.1635|0.2478|0|0|0|0|3|0|0|0|0|0|1|5|4|GP|02|07|$GPEtnSn1AlkylTemplateString",

			       "PENMe2Sn1Sn2" => "PENMe2|glycero-3-phospho-N,N-dimethylethanolamine|0.3628|0.7753|-0.7737|-0.3611|0|0|8|20|0|2|2|0|1|7|6|GP|02|01|$GPEtnNMe2Sn1Sn2TemplateString",
			       "PENMe2Sn1AlkylSn2" => "PENMe2|glycero-3-phospho-N,N-dimethylethanolamine|0.4607|0.8719|-0.6759|-0.2632|0|0|3|17|0|0|2|0|1|5|4|GP|02|02|$GPEtnNMe2Sn1AlkylSn2TemplateString",
			       "PENMe2Sn1AlkenylSn2" => "PENMe2|glycero-3-phospho-N,N-dimethylethanolamine|0.4607|0.8719|-0.6759|-0.2632|0|0|3|17|0|0|2|0|1|5|4|GP|02|03|$GPEtnNMe2Sn1AlkylSn2TemplateString",
			       "PENMe2Sn1AlkylSn2Alkyl" => "PENMe2|glycero-3-phospho-N,N-dimethylethanolamine| 0.1067|0.5179|-1.0299|-0.6077|0|0|3|15|0|0|1|0|1|5|4|GP|02|04|$GPEtnNMe2Sn1AlkylSn2AlkylTemplateString",
			       "PENMe2Sn1" => "PENMe2|glycero-3-phospho-N,N-dimethylethanolamine|0.0088|0.4214|0|0|0|0|8|0|0|2|0|0|1|7|6|GP|02|05|$GPEtnNMe2Sn1TemplateString",
			       "PENMe2Sn1Alkyl" => "PENMe2|glycero-3-phospho-N,N-dimethylethanolamine|0.1066|0.5179|0|0|0|0|3|0|0|0|0|0|1|5|4|GP|02|06|$GPEtnNMe2Sn1AlkylTemplateString",
			       "PENMe2Sn1Alkenyl" => "PENMe2|glycero-3-phospho-N,N-dimethylethanolamine|0.1066|0.5179|0|0|0|0|3|0|0|0|0|0|1|5|4|GP|02|07|$GPEtnNMe2Sn1AlkylTemplateString",

			       "PENMeSn1Sn2" => "PENMe|glycero-3-phospho-N-methylethanolamine|0.3628|0.7753|-0.7737|-0.3611|0|0|8|20|0|2|2|0|1|7|6|GP|02|01|$GPEtnNMeSn1Sn2TemplateString",
			       "PENMeSn1AlkylSn2" => "PENMe|glycero-3-phospho-N-methylethanolamine|0.4607|0.8719|-0.6759|-0.2632|0|0|3|17|0|0|2|0|1|5|4|GP|02|02|$GPEtnNMeSn1AlkylSn2TemplateString",
			       "PENMeSn1AlkenylSn2" => "PENMe|glycero-3-phospho-N-methylethanolamine|0.4607|0.8719|-0.6759|-0.2632|0|0|3|17|0|0|2|0|1|5|4|GP|02|03|$GPEtnNMeSn1AlkylSn2TemplateString",
			       "PENMeSn1AlkylSn2Alkyl" => "PENMe|glycero-3-phospho-N-methylethanolamine| 0.1067|0.5179|-1.0299|-0.6077|0|0|3|15|0|0|1|0|1|5|4|GP|02|04|$GPEtnNMeSn1AlkylSn2AlkylTemplateString",
			       "PENMeSn1" => "PENMe|glycero-3-phospho-N-methylethanolamine|-0.2613|0.1499|0|0|0|0|8|0|0|2|0|0|1|7|6|GP|02|05|$GPEtnNMeSn1TemplateString",
			       "PENMeSn1Alkyl" => "PENMe|glycero-3-phospho-N-methylethanolamine|-0.1635|0.2478|0|0|0|0|3|0|0|0|0|0|1|5|4|GP|02|06|$GPEtnNMeSn1AlkylTemplateString",
			       "PENMeSn1Alkenyl" => "PENMe|glycero-3-phospho-N-methylethanolamine|-0.1635|0.2478|0|0|0|0|3|0|0|0|0|0|1|5|4|GP|02|07|$GPEtnNMeSn1AlkylTemplateString",

			       "PSSn1Sn2" => "PS|glycero-3-phosphoserine|0.3124|0.7249|-0.8242|-0.4115|0|0|9|21|0|2|2|0|2|8|7|GP|03|01|$GPSerSn1Sn2TemplateString",
			       "PSSn1AlkylSn2" => "PS|glycero-3-phosphoserine|0.3124|0.7236|-0.8242|-0.4115|0|0|4|18|0|0|2|0|2|6|5|GP|03|02|$GPSerSn1AlkylSn2TemplateString",
			       "PSSn1AlkenylSn2" => "PS|glycero-3-phosphoserine|0.3124|0.7236|-0.8242|-0.4115|0|0|4|18|0|0|2|0|2|6|5|GP|03|03|$GPSerSn1AlkylSn2TemplateString",
			       "PSSn1AlkylSn2Alkyl" => "PS|glycero-3-phosphoserine|-0.1006|0.3106|-1.2372|-0.8150|0|0|4|16|0|0|1|0|2|6|5|GP|03|04|$GPSerSn1AlkylSn2AlkylTemplateString",
			       "PSSn1" => "PS|glycero-3-phosphoserine|-0.3117|0.1007|0|0|0|0|9|0|0|2|0|0|2|8|7|GP|03|05|$GPSerSn1TemplateString",
			       "PSSn1Alkyl" => "PS|glycero-3-phosphoserine|-0.3117|0.0995|0|0|0|0|4|0|0|0|0|0|2|6|5|GP|03|06|$GPSerSn1AlkylTemplateString",
			       "PSSn1Alkenyl" => "PS|glycero-3-phosphoserine|-0.3117|0.0995|0|0|0|0|4|0|0|0|0|0|2|6|5|GP|03|07|$GPSerSn1AlkylTemplateString",


			       "PGSn1Sn2" => "PG|glycero-3-phospho-(1\'-sn-glycerol)|0.3124|0.7249|-0.8242|-0.4115|0|0|9|21|0|2|2|0|2|8|7|GP|04|01|$GPGroSn1Sn2TemplateString",
			       "PGSn1AlkylSn2" => "PG|glycero-3-phospho-(1\'-sn-glycerol)|0.4676|0.8788|-0.6692|-0.2564|0|0|3|17|0|0|2|0|1|5|4|GP|04|02|$GPGroSn1AlkylSn2TemplateString",
			       "PGSn1AlkenylSn2" => "PG|glycero-3-phospho-(1\'-sn-glycerol)|0.4676|0.8788|-0.6692|-0.2564|0|0|3|17|0|0|2|0|1|5|4|GP|04|03|$GPGroSn1AlkylSn2TemplateString",
			       "PGSn1AlkylSn2Alkyl" => "PG|glycero-3-phospho-(1\'-sn-glycerol)|0.0547|0.4659|-1.0821|-0.6598|0|0|3|15|0|0|1|0|1|5|4|GP|04|04|$GPGroSn1AlkylSn2AlkylTemplateString",
			       "PGSn1" => "PG|glycero-3-phospho-(1\'-sn-glycerol)|-0.3117|0.1007|0|0|0|0|9|0|0|2|0|0|2|8|7|GP|04|05|$GPGroSn1TemplateString",
			       "PGSn1Alkyl" => "PG|glycero-3-phospho-(1\'-sn-glycerol)|-0.1564|0.2547|0|0|0|0|3|0|0|0|0|0|1|5|4|GP|04|06|$GPGroSn1AlkylTemplateString",
			       "PGSn1Alkenyl" => "PG|glycero-3-phospho-(1\'-sn-glycerol)|-0.1564|0.2547|0|0|0|0|3|0|0|0|0|0|1|5|4|GP|04|07|$GPGroSn1AlkylTemplateString",

			       "PGPSn1Sn2" => "PGP|glycero-3-phospho-(1\'-sn-glycerol-3\'-phosphate)|0.3124|0.7249|-0.8242|-0.4115|0|0|9|21|0|2|2|0|2|8|7|GP|05|01|$GPGroPSn1Sn2TemplateString",
			       "PGPSn1AlkylSn2" => "PGP|glycero-3-phospho-(1\'-sn-glycerol-3\'-phosphate)|0.4676|0.8788|-0.6692|-0.2564|0|0|3|17|0|0|2|0|1|5|4|GP|05|02|$GPGroPSn1AlkylSn2TemplateString",
			       "PGPSn1AlkenylSn2" => "PGP|glycero-3-phospho-(1\'-sn-glycerol-3\'-phosphate)|0.4676|0.8788|-0.6692|-0.2564|0|0|3|17|0|0|2|0|1|5|4|GP|05|03|$GPGroPSn1AlkylSn2TemplateString",
			       "PGPSn1AlkylSn2Alkyl" => "PGP|glycero-3-phospho-(1\'-sn-glycerol-3\'-phosphate)|0.0547|0.4659|-1.0821|-0.6598|0|0|3|15|0|0|1|0|1|5|4|GP|05|04|$GPGroPSn1AlkylSn2AlkylTemplateString",
			       "PGPSn1" => "PGP|glycero-3-phospho-(1\'-sn-glycerol-3\'-phosphate)|-0.3117|0.1007|0|0|0|0|9|0|0|2|0|0|2|8|7|GP|05|05|$GPGroPSn1TemplateString",
			       "PGPSn1Alkyl" => "PGP|glycero-3-phospho-(1\'-sn-glycerol-3\'-phosphate)|-0.1564|0.2547|0|0|0|0|3|0|0|0|0|0|1|5|4|GP|05|06|$GPGroPSn1AlkylTemplateString",
			       "PGPSn1Alkenyl" => "PGP|glycero-3-phospho-(1\'-sn-glycerol-3\'-phosphate)|-0.1564|0.2547|0|0|0|0|3|0|0|0|0|0|1|5|4|GP|05|07|$GPGroPSn1AlkylTemplateString",


			       "PISn1Sn2" => "PI|glycero-3-phospho-(1\'-myo-inositol)|0.3628|0.7740|-0.7738|-0.3611|0|0|8|17|0|2|2|0|1|7|6|GP|06|01|$GPInsSn1Sn2TemplateString",
			       "PISn1AlkylSn2" => "PI|glycero-3-phospho-(1\'-myo-inositol)|0.4627|0.8739|-0.6737|-0.2611|0|0|3|14|0|0|2|0|1|5|4|GP|06|02|$GPInsSn1AlkylSn2TemplateString",
			       "PISn1AlkenylSn2" => "PI|glycero-3-phospho-(1\'-myo-inositol)|0.4627|0.8739|-0.6737|-0.2611|0|0|3|14|0|0|2|0|1|5|4|GP|06|03|$GPInsSn1AlkylSn2TemplateString",
			       "PISn1AlkylSn2Alkyl" => "PI|glycero-3-phospho-(1\'-myo-inositol)|0.0498|0.4609|-1.0865|-0.6644|0|0|3|12|0|0|1|0|1|5|4|GP|06|04|$GPInsSn1AlkylSn2AlkylTemplateString",
			       "PISn1" => "PI|glycero-3-phospho-(1\'-myo-inositol)|-0.2612|0.1540|0|0|0|0|8|0|0|2|0|0|1|7|6|GP|06|05|$GPInsSn1TemplateString",
			       "PISn1Alkyl" => "PI|glycero-3-phospho-(1\'-myo-inositol)|-0.1612|0.2499|0|0|0|0|3|0|0|0|0|0|1|5|4|GP|06|06|$GPInsSn1AlkylTemplateString",
			       "PISn1Alkenyl" => "PI|glycero-3-phospho-(1\'-myo-inositol)|-0.1612|0.2499|0|0|0|0|3|0|0|0|0|0|1|5|4|GP|06|07|$GPInsSn1AlkylTemplateString",

			       "PIPSn1Sn2" => "PIP|glycero-3-phospho-(1\'-myo-inositol-3\'-phosphate)|0.3628|0.7740|-0.7738|-0.3611|0|0|8|17|0|2|2|0|1|7|6|GP|07|01|$GPInsPSn1Sn2TemplateString",
			       "PIPSn1AlkylSn2" => "PIP|glycero-3-phospho-(1\'-myo-inositol-3\'-phosphate)|0.4627|0.8739|-0.6737|-0.2611|0|0|3|14|0|0|2|0|1|5|4|GP|07|02|$GPInsPSn1AlkylSn2TemplateString",
			       "PIPSn1AlkenylSn2" => "PIP|glycero-3-phospho-(1\'-myo-inositol-3\'-phosphate)|0.4627|0.8739|-0.6737|-0.2611|0|0|3|14|0|0|2|0|1|5|4|GP|07|03|$GPInsPSn1AlkylSn2TemplateString",
			       "PIPSn1AlkylSn2Alkyl" => "PIP|glycero-3-phospho-(1\'-myo-inositol-3\'-phosphate)|0.0498|0.4609|-1.0865|-0.6644|0|0|3|12|0|0|1|0|1|5|4|GP|07|04|$GPInsPSn1AlkylSn2AlkylTemplateString",
			       "PIPSn1" => "PIP|glycero-3-phospho-(1\'-myo-inositol-3\'-phosphate)|-0.2612|0.1540|0|0|0|0|8|0|0|2|0|0|1|7|6|GP|07|05|$GPInsPSn1TemplateString",
			       "PIPSn1Alkyl" => "PIP|glycero-3-phospho-(1\'-myo-inositol-3\'-phosphate)|-0.1612|0.2499|0|0|0|0|3|0|0|0|0|0|1|5|4|GP|07|06|$GPInsPSn1AlkylTemplateString",
			       "PIPSn1Alkenyl" => "PIP|glycero-3-phospho-(1\'-myo-inositol-3\'-phosphate)|-0.1612|0.2499|0|0|0|0|3|0|0|0|0|0|1|5|4|GP|07|07|$GPInsPSn1AlkylTemplateString",

			       "PIP[3\']Sn1Sn2" => "PIP[3\']|glycero-3-phospho-(1\'-myo-inositol-3\'-phosphate)|0.3628|0.7739|-0.7738|-0.3610|0|0|8|17|0|2|2|0|1|7|6|GP|07|01|$GPInsP3Sn1Sn2TemplateString",
			       "PIP[3\']Sn1AlkylSn2" => "PIP[3\']|glycero-3-phospho-(1\'-myo-inositol-3\'-phosphate)|0.4627|0.8737|-0.6737|-0.2609|0|0|3|14|0|0|2|0|1|5|4|GP|07|02|$GPInsP3Sn1AlkylSn2TemplateString",
			       "PIP[3\']Sn1AlkenylSn2" => "PIP[3\']|glycero-3-phospho-(1\'-myo-inositol-3\'-phosphate)|0.4627|0.8737|-0.6737|-0.2609|0|0|3|14|0|0|2|0|1|5|4|GP|07|03|$GPInsP3Sn1AlkylSn2TemplateString",
			       "PIP[3\']Sn1AlkylSn2Alkyl" => "PIP[3\']|glycero-3-phospho-(1\'-myo-inositol-3\'-phosphate)|0.2057|0.6167|-0.9307|-0.5086|0|0|3|12|0|0|1|0|1|5|4|GP|07|04|$GPInsP3Sn1AlkylSn2AlkylTemplateString",
			       "PIP[3\']Sn1" => "PIP[3\']|glycero-3-phospho-(1\'-myo-inositol-3\'-phosphate)|0.1057|0.5167|0|0|0|0|8|0|0|2|0|0|1|7|6|GP|07|05|$GPInsP3Sn1TemplateString",
			       "PIP[3\']Sn1Alkyl" => "PIP[3\']|glycero-3-phospho-(1\'-myo-inositol-3\'-phosphate)|0.2057|0.6168|0|0|0|0|3|0|0|0|0|0|1|5|4|GP|07|06|$GPInsP3Sn1AlkylTemplateString",
			       "PIP[3\']Sn1Alkenyl" => "PIP[3\']|glycero-3-phospho-(1\'-myo-inositol-3\'-phosphate)|0.2057|0.6168|0|0|0|0|3|0|0|0|0|0|1|5|4|GP|07|07|$GPInsP3Sn1AlkylTemplateString",

			       "PIP[4\']Sn1Sn2" => "PIP[4\']|glycero-3-phospho-(1\'-myo-inositol-4\'-phosphate)|0.3148|0.7259|-0.8219|-0.4091|0|0|8|17|0|2|2|0|1|7|6|GP|07|01|$GPInsP4Sn1Sn2TemplateString",
			       "PIP[4\']Sn1AlkylSn2" => "PIP[4\']|glycero-3-phospho-(1\'-myo-inositol-4\'-phosphate)|0.3147|0.7257|-0.8217|-0.4090|0|0|3|14|0|0|2|0|1|5|4|GP|07|02|$GPInsP4Sn1AlkylSn2TemplateString",
			       "PIP[4\']Sn1AlkenylSn2" => "PIP[4\']|glycero-3-phospho-(1\'-myo-inositol-4\'-phosphate)|0.3147|0.7257|-0.8217|-0.4090|0|0|3|14|0|0|2|0|1|5|4|GP|07|03|$GPInsP4Sn1AlkylSn2TemplateString",
			       "PIP[4\']Sn1AlkylSn2Alkyl" => "PIP[4\']|glycero-3-phospho-(1\'-myo-inositol-4\'-phosphate)|-0.0982|0.3129|-1.2346|-0.8125|0|0|3|12|0|0|1|0|1|5|4|GP|07|04|$GPInsP4Sn1AlkylSn2AlkylTemplateString",
			       "PIP[4\']Sn1" => "PIP[4\']|glycero-3-phospho-(1\'-myo-inositol-4\'-phosphate)|-0.3092|0.1018|0|0|0|0|8|0|0|2|0|0|1|7|6|GP|07|05|$GPInsP4Sn1TemplateString",
			       "PIP[4\']Sn1Alkyl" => "PIP[4\']|glycero-3-phospho-(1\'-myo-inositol-4\'-phosphate)|-0.3093|0.1018|0|0|0|0|3|0|0|0|0|0|1|5|4|GP|07|06|$GPInsP4Sn1AlkylTemplateString",
			       "PIP[4\']Sn1Alkenyl" => "PIP[4\']|glycero-3-phospho-(1\'-myo-inositol-4\'-phosphate)|-0.3093|0.1018|0|0|0|0|3|0|0|0|0|0|1|5|4|GP|07|07|$GPInsP4Sn1AlkylTemplateString",

			       "PIP[5\']Sn1Sn2" => "PIP[5\']|glycero-3-phospho-(1\'-myo-inositol-5\'-phosphate)|-0.0479|0.3632|-1.1846|-0.7718|0|0|8|17|0|2|2|0|1|7|6|GP|07|01|$GPInsP5Sn1Sn2TemplateString",
			       "PIP[5\']Sn1AlkylSn2" => "PIP[5\']|glycero-3-phospho-(1\'-myo-inositol-5\'-phosphate)|-0.0479|0.3631|-1.1843|-0.7716|0|0|3|14|0|0|2|0|1|5|4|GP|07|02|$GPInsP5Sn1AlkylSn2TemplateString",
			       "PIP[5\']Sn1AlkenylSn2" => "PIP[5\']|glycero-3-phospho-(1\'-myo-inositol-5\'-phosphate)|-0.0479|0.3631|-1.1843|-0.7716|0|0|3|14|0|0|2|0|1|5|4|GP|07|03|$GPInsP5Sn1AlkylSn2TemplateString",
			       "PIP[5\']Sn1AlkylSn2Alkyl" => "PIP[5\']|glycero-3-phospho-(1\'-myo-inositol-5\'-phosphate)|-0.4608|-0.0498|-1.5973|-1.1752|0|0|3|12|0|0|1|0|1|5|4|GP|07|04|$GPInsP5Sn1AlkylSn2AlkylTemplateString",
			       "PIP[5\']Sn1" => "PIP[5\']|glycero-3-phospho-(1\'-myo-inositol-5\'-phosphate)|-0.6718|-0.2608|0|0|0|0|8|0|0|2|0|0|1|7|6|GP|07|05|$GPInsP5Sn1TemplateString",
			       "PIP[5\']Sn1Alkyl" => "PIP[5\']|glycero-3-phospho-(1\'-myo-inositol-5\'-phosphate)|-0.6719|-0.2609|0|0|0|0|3|0|0|0|0|0|1|5|4|GP|07|06|$GPInsP5Sn1AlkylTemplateString",
			       "PIP[5\']Sn1Alkenyl" => "PIP[5\']|glycero-3-phospho-(1\'-myo-inositol-5\'-phosphate)|-0.6719|-0.2609|0|0|0|0|3|0|0|0|0|0|1|5|4|GP|07|07|$GPInsP5Sn1AlkylTemplateString",

			       "PIP2[3\',4\']Sn1Sn2" => "PIP2[3\',4\']|glycero-3-phospho-(1\'-myo-inositol-3\',4\'-biphosphate)|0.3148|0.7259|-0.8219|-0.4091|0|0|8|17|0|2|2|0|1|7|6|GP|08|01|$GPInsP34Sn1Sn2TemplateString",
			       "PIP2[3\',4\']Sn1AlkylSn2" => "PIP2[3\',4\']|glycero-3-phospho-(1\'-myo-inositol-3\',4\'-biphosphate)|0.3148|0.7259|-0.8219|-0.4091|0|0|3|14|0|0|2|0|1|5|4|GP|08|02|$GPInsP34Sn1AlkylSn2TemplateString",
			       "PIP2[3\',4\']Sn1AlkenylSn2" => "PIP2[3\',4\']|glycero-3-phospho-(1\'-myo-inositol-3\',4\'-biphosphate)|0.3148|0.7259|-0.8219|-0.4091|0|0|3|14|0|0|2|0|1|5|4|GP|08|03|$GPInsP34Sn1AlkylSn2TemplateString",
			       "PIP2[3\',4\']Sn1" => "PIP2[3\',4\']|glycero-3-phospho-(1\'-myo-inositol-3\',4\'-biphosphate)|0.0577|0.4688|0|0|0|0|8|0|0|2|0|0|1|7|6|GP|08|04|$GPInsP34Sn1TemplateString",
			       "PIP2[3\',4\']Sn1Alkyl" => "PIP2[3\',4\']|glycero-3-phospho-(1\'-myo-inositol-3\',4\'-biphosphate)|0.0577|0.4688|0|0|0|0|3|0|0|0|0|0|1|5|4|GP|08|05|$GPInsP34Sn1AlkylTemplateString",
			       "PIP2[3\',4\']Sn1Alkenyl" => "PIP2[3\',4\']|glycero-3-phospho-(1\'-myo-inositol-3\',4\'-biphosphate)|0.0577|0.4688|0|0|0|0|3|0|0|0|0|0|1|5|4|GP|08|06|$GPInsP34Sn1AlkylTemplateString",
			       "PIP2[3\',4\']Sn1AlkylSn2Alkyl" => "PIP2[3\',4\']|glycero-3-phospho-(1\'-myo-inositol-3\',4\'-biphosphate)|0.0577|0.4688|-1.0789|-0.6567|0|0|3|12|0|0|1|0|1|5|4|GP|08||$GPInsP34Sn1AlkylSn2AlkylTemplateString",

			       "PIP2[3\',5\']Sn1Sn2" => "PIP2[3\',5\']|glycero-3-phospho-(1\'-myo-inositol-3\',5\'-biphosphate)|-0.0479|0.3632|-1.1846|-0.7718|0|0|8|17|0|2|2|0|1|7|6|GP|08|01|$GPInsP35Sn1Sn2TemplateString",
			       "PIP2[3\',5\']Sn1AlkylSn2" => "PIP2[3\',5\']|glycero-3-phospho-(1\'-myo-inositol-3\',5\'-biphosphate)|-0.0479|0.3632|-1.1846|-0.7718|0|0|3|14|0|0|2|0|1|5|4|GP|08|02|$GPInsP35Sn1AlkylSn2TemplateString",
			       "PIP2[3\',5\']Sn1AlkenylSn2" => "PIP2[3\',5\']|glycero-3-phospho-(1\'-myo-inositol-3\',5\'-biphosphate)|-0.0479|0.3632|-1.1846|-0.7718|0|0|3|14|0|0|2|0|1|5|4|GP|08|03|$GPInsP35Sn1AlkylSn2TemplateString",
			       "PIP2[3\',5\']Sn1" => "PIP2[3\',5\']|glycero-3-phospho-(1\'-myo-inositol-3\',5\'-biphosphate)|-0.3050|0.1061|0|0|0|0|8|0|0|2|0|0|1|7|6|GP|08|04|$GPInsP35Sn1TemplateString",
			       "PIP2[3\',5\']Sn1Alkyl" => "PIP2[3\',5\']|glycero-3-phospho-(1\'-myo-inositol-3\',5\'-biphosphate)|-0.3050|0.1061|0|0|0|0|3|0|0|0|0|0|1|5|4|GP|08|05|$GPInsP35Sn1AlkylTemplateString",
			       "PIP2[3\',5\']Sn1Alkenyl" => "PIP2[3\',5\']|glycero-3-phospho-(1\'-myo-inositol-3\',5\'-biphosphate)|-0.3050|0.1061|0|0|0|0|3|0|0|0|0|0|1|5|4|GP|08|06|$GPInsP35Sn1AlkylTemplateString",
			       "PIP2[3\',5\']Sn1AlkylSn2Alkyl" => "PIP2[3\',5\']|glycero-3-phospho-(1\'-myo-inositol-3\',5\'-biphosphate)|-0.3050|0.1061|-1.4417|-1.0195|0|0|3|12|0|0|1|0|1|5|4|GP|08|0|$GPInsP35Sn1AlkylSn2AlkylTemplateString",

			       "PIP2[4\',5\']Sn1Sn2" => "PIP2[4\',5\']|glycero-3-phospho-(1\'-myo-inositol-4\',5\'-biphosphate)|-0.0479|0.3632|-1.1846|-0.7718|0|0|8|17|0|2|2|0|1|7|6|GP|08|01|$GPInsP45Sn1Sn2TemplateString",
			       "PIP2[4\',5\']Sn1AlkylSn2" => "PIP2[4\',5\']|glycero-3-phospho-(1\'-myo-inositol-4\',5\'-biphosphate)|-0.0479|0.3632|-1.1846|-0.7718|0|0|3|14|0|0|2|0|1|5|4|GP|08|02|$GPInsP45Sn1AlkylSn2TemplateString",
			       "PIP2[4\',5\']Sn1AlkenylSn2" => "PIP2[4\',5\']|glycero-3-phospho-(1\'-myo-inositol-4\',5\'-biphosphate)|-0.0479|0.3632|-1.1846|-0.7718|0|0|3|14|0|0|2|0|1|5|4|GP|08|03|$GPInsP45Sn1AlkylSn2TemplateString",
			       "PIP2[4\',5\']Sn1" => "PIP2[4\',5\']|glycero-3-phospho-(1\'-myo-inositol-4\',5\'-biphosphate)|-0.6720|-0.2609|0|0|0|0|8|0|0|2|0|0|1|7|6|GP|08|04|$GPInsP45Sn1TemplateString",
			       "PIP2[4\',5\']Sn1Alkyl" => "PIP2[4\',5\']|glycero-3-phospho-(1\'-myo-inositol-4\',5\'-biphosphate)|-0.6720|-0.2609|0|0|0|0|3|0|0|0|0|0|1|5|4|GP|08|05|$GPInsP45Sn1AlkylTemplateString",
			       "PIP2[4\',5\']Sn1Alkenyl" => "PIP2[4\',5\']|glycero-3-phospho-(1\'-myo-inositol-4\',5\'-biphosphate)|-0.6720|-0.2609|0|0|0|0|3|0|0|0|0|0|1|5|4|GP|08|06|$GPInsP45Sn1AlkylTemplateString",
			       "PIP2[4\',5\']Sn1AlkylSn2Alkyl" => "PIP2[4\',5\']|glycero-3-phospho-(1\'-myo-inositol-4\',5\'-biphosphate)|-0.4609|-0.0498|-1.5976|-1.1754|0|0|3|12|0|0|1|0|1|5|4|GP|08|0|$GPInsP45Sn1AlkylSn2AlkylTemplateString",

			       "PIP3[3\',4\',5\']Sn1Sn2" => "PIP3[3\',4\',5\']|glycero-3-phospho-(1\'-myo-inositol-3\',4\',5\'-triphosphate)|-0.0479|0.3632|-1.1846|-0.7718|0|0|8|17|0|2|2|0|1|7|6|GP|09|01|$GPInsP345Sn1Sn2TemplateString",
			       "PIP3[3\',4\',5\']Sn1AlkylSn2" => "PIP3[3\',4\',5\']|glycero-3-phospho-(1\'-myo-inositol-3\',4\',5\'-triphosphate)|-0.0479|0.3632|-1.1846|-0.7718|0|0|3|14|0|0|2|0|1|5|4|GP|09|02|$GPInsP345Sn1AlkylSn2TemplateString",
			       "PIP3[3\',4\',5\']Sn1AlkenylSn2" => "PIP3[3\',4\',5\']|glycero-3-phospho-(1\'-myo-inositol-3\',4\',5\'-triphosphate)|-0.0479|0.3632|-1.1846|-0.7718|0|0|3|14|0|0|2|0|1|5|4|GP|09|03|$GPInsP345Sn1AlkylSn2TemplateString",
			       "PIP3[3\',4\',5\']Sn1" => "PIP3[3\',4\',5\']|glycero-3-phospho-(1\'-myo-inositol-3\',4\',5\'-triphosphate)|-0.3050|0.1061|0|0|0|0|8|0|0|2|0|0|1|7|6|GP|09|04|$GPInsP345Sn1TemplateString",
			       "PIP3[3\',4\',5\']Sn1Alkyl" => "PIP3[3\',4\',5\']|glycero-3-phospho-(1\'-myo-inositol-3\',4\',5\'-triphosphate)|-0.3050|0.1061|0|0|0|0|3|0|0|0|0|0|1|5|4|GP|09|05|$GPInsP345Sn1AlkylTemplateString",
			       "PIP3[3\',4\',5\']Sn1Alkenyl" => "PIP3[3\',4\',5\']|glycero-3-phospho-(1\'-myo-inositol-3\',4\',5\'-triphosphate)|-0.3050|0.1061|0|0|0|0|3|0|0|0|0|0|1|5|4|GP|09|06|$GPInsP345Sn1AlkylTemplateString",
			       "PIP3[3\',4\',5\']Sn1AlkylSn2Alkyl" => "PIP3[3\',4\',5\']|glycero-3-phospho-(1\'-myo-inositol-3\',4\',5\'-triphosphate)|-0.3050|0.1061|-1.4417|-1.0195|0|0|3|12|0|0|1|0|1|5|4|GP|09|0|$GPInsP345Sn1AlkylSn2AlkylTemplateString",

			       "PASn1Sn2" => "PA|glycero-3-phosphate|0.3628|0.7753|-0.7737|-0.3611|0|0|8|17|0|2|2|0|1|7|6|GP|10|01|$GPASn1Sn2TemplateString",
			       "PASn1AlkylSn2" => "PA|glycero-3-phosphate|0.4607|0.8719|-0.6759|-0.2632|0|0|3|14|0|0|2|0|1|5|4|GP|10|02|$GPASn1AlkylSn2TemplateString",
			       "PASn1AlkenylSn2" => "PA|glycero-3-phosphate|0.4607|0.8719|-0.6759|-0.2632|0|0|3|14|0|0|2|0|1|5|4|GP|10|03|$GPASn1AlkylSn2TemplateString",
			       "PASn1AlkylSn2Alkyl" => "PA|glycero-3-phosphate|0.0477|0.4589|-1.0890|-0.6667|0|0|3|12|0|0|1|0|1|5|4|GP|10|04|$GPASn1AlkylSn2AlkylTemplateString",
			       "PASn1" => "PA|glycero-3-phosphate|-0.2613|0.15120|0|0|0|0|8|0|0|2|0|0|1|7|6|GP|10|05|$GPASn1TemplateString",
			       "PASn1Alkyl" => "PA|glycero-3-phosphate|-0.1634|0.2478|0|0|0|0|3|0|0|0|0|0|1|5|4|GP|10|06|$GPASn1AlkylTemplateString",
			       "PASn1Alkenyl" => "PA|glycero-3-phosphate|-0.1634|0.2478|0|0|0|0|3|0|0|0|0|0|1|5|4|GP|10|07|$GPASn1AlkylTemplateString",

			       "PPASn1Sn2" => "PPA|glycero-3-pyrophosphate|0.3628|0.7753|-0.7737|-0.3611|0|0|8|17|0|2|2|0|1|7|6|GP|11|01|$GPPSn1Sn2TemplateString",
			       "PPASn1" => "PPA|glycero-3-pyrophosphate|-0.2613|0.15120|0|0|0|0|8|0|0|2|0|0|1|7|6|GP|11|02|$GPPSn1TemplateString",

			       "PnCSn1Sn2" => "PnC|glycero-3-phosphonocholine|0.3628|0.7753|-0.7737|-0.3611|0|0|8|16|0|2|2|0|1|7|6|GP|16|01|$GPnChoSn1Sn2TemplateString",
			       "PnCSn1AlkylSn2" => "PnC|glycero-3-phosphonocholine|0.4606|0.8718|-0.6759|-0.2632|0|0|3|13|0|0|2|0|1|5|4|GP|16|02|$GPnChoSn1AlkylSn2TemplateString",
			       "PnCSn1AlkenylSn2" => "PnC|glycero-3-phosphonocholine|0.4606|0.8718|-0.6759|-0.2632|0|0|3|13|0|0|2|0|1|5|4|GP|16|03|$GPnChoSn1AlkylSn2TemplateString",
			       "PnCSn1AlkylSn2Alkyl" => "PnC|glycero-3-phosphonocholine|0.1067|0.5179|-1.0298|-0.6076|0|0|3|11|0|0|1|0|1|5|4|GP|16|04|$GPnChoSn1AlkylSn2AlkylTemplateString",
			       "PnCSn1" => "PnC|glycero-3-phosphonocholine|-0.1443|0.2682|0|0|0|0|8|0|0|2|0|0|1|7|6|GP|16|05|$GPnChoSn1TemplateString",
			       "PnCSn1Alkyl" => "PnC|glycero-3-phosphonocholine|-0.0464|0.3648|0|0|0|0|3|0|0|0|0|0|1|5|4|GP|16|06|$GPnChoSn1AlkylTemplateString",

			       "PnCSn1Alkenyl" => "PnC|glycero-3-phosphonocholine|-0.0464|0.3648|0|0|0|0|3|0|0|0|0|0|1|5|4|GP|16|07|$GPnChoSn1AlkylTemplateString",

			       "PnESn1Sn2" => "PnE|glycero-3-phosphonoethanolamine|0.3628|0.7753|-0.7737|-0.3611|0|0|8|16|0|2|2|0|1|7|6|GP|17|01|$GPnEtnSn1Sn2TemplateString",
			       "PnESn1AlkylSn2" => "PnE|glycero-3-phosphonoethanolamine|0.4607|0.8719|-0.6759|-0.2632|0|0|3|13|0|0|2|0|1|5|4|GP|17|02|$GPnEtnSn1AlkylSn2TemplateString",
			       "PnESn1AlkenylSn2" => "PnE|glycero-3-phosphonoethanolamine|0.4607|0.8719|-0.6759|-0.2632|0|0|3|13|0|0|2|0|1|5|4|GP|17|03|$GPnEtnSn1AlkylSn2TemplateString",
			       "PnESn1AlkylSn2Alkyl" => "PnE|glycero-3-phosphonoethanolamine|0.0477|0.4588|-1.0890|-0.6668|0|0|3|11|0|0|1|0|1|5|4|GP|17|04|$GPnEtnSn1AlkylSn2AlkylTemplateString",
			       "PnESn1" => "PnE|glycero-3-phosphonoethanolamine|-0.2613|0.1499|0|0|0|0|8|0|0|2|0|0|1|7|6|GP|17|05|$GPnEtnSn1TemplateString",
			       "PnESn1Alkyl" => "PnE|glycero-3-phosphonoethanolamine|0.3438|0.7551|0|0|0|0|3|0|0|0|0|0|1|5|4|GP|17|06|$GPnEtnSn1AlkylTemplateString",
			       "PnESn1Alkenyl" => "PnE|glycero-3-phosphonoethanolamine|0.3438|0.7551|0|0|0|0|3|0|0|0|0|0|1|5|4|GP|17|07|$GPnEtnSn1AlkylTemplateString",

			       );
}
# Template format: ID => AbbrevID|HeadGroupName|Y1Sn1|Y2Sn1|Y1Sn2|Y2Sn2|Y1Sn3|Y2Sn3|Sn1ChainAtomNum|Sn2ChainAtomNum|Sn3ChainAtomNum|Sn1ChainCarbons|Sn2ChainCarbons|Sn3ChainCarbons|Sn2CAtomNum|Sn2OAtomNum|Sn2HAtomNum|LMCategory|LMMainClass|LMSubClass|TemplateString

# Initialize supported head groups...
sub _InitializeSupportedHeadGroupsData {
  my($GPType, $GPHeadGroup);
  %GPSupportedHeadGroupMap = ();

  for $GPType (keys %GPTemplatesDataMap) {
    ($GPHeadGroup) = split /\|/, $GPTemplatesDataMap{$GPType};
    if (!(exists $GPSupportedHeadGroupMap{$GPHeadGroup})) {
      $GPSupportedHeadGroupMap{$GPHeadGroup} = $GPHeadGroup;
    }
  }
}


1;

__END__

=head1 NAME

GPStr - Glycerolipids (GP) structure generation methods

=head1 SYNOPSIS

use GPStr;

use GPStr qw(:all);

=head1 DESCRIPTION

GPStr module provides these methods:

    GenerateCmpdOntologyData - Generate ontology data
    GenerateCmpdOntologySDDataLines - Generate ontology data lines for SD file
    GenerateGPChainStrData - Generate chain structure data
    GenerateSDFile - Generate SD file
    GetGPTemplatesData - Get templates data
    GetGPSupportedHeadGroupMap - Get supported headgroups data
    GetGPTemplateID - Get templates ID
    IsGPChainsAbbrevSupported - Is it a supported GP abbreviation
    ParseGPAbbrev - Parse GP abbreviation
    ProcessGPCmpdAbbrevs - Process GP abbreviation
    SetupGPCmpdAbbrevTemplateDataMap - Setup template structure data map
    ValidateGPAbbrev - Validate GP abbreviation

=head1 METHODS

=over 4

=item B<GenerateCmpdOntologyData>

    $DataHashRef = GenerateCmpdOntologyData($CmpDataRef);

Return a reference to a hash containing ontology data with hash keys and values
corresponding to property names and values.

=item B<GenerateCmpdOntologySDDataLines>

    $DataLinesArrayRef = GenerateCmpdOntologySDDataLines($CmpdDataRef);

Return a reference to an array containing ontology data lines suitable for
generate SD file data block.

=item B<GenerateGPChainStrData>

    ($AtomLinesArrayRef, $BondLinesArrayRef) =
       GenerateGPChainStrData($ChainType, $CmpdDataRef);

Return array references containing atom and bond data lines for SD file. Appropriate atom
and bond data lines are generated using chain type and abbreviation template data.

=item B<GenerateSDFile>

    GenerateSDFile($SDFileName, $CmdAbbrevsRef);

Generate a SD file for compound abbreviations. Structure data for specified abbreviation
is generated sequentially and written to SD file.

=item B<GetGPTemplatesData>

    $TemplatesDataRef = GetGPTemplatesData();

Return a reference to a hash containing GP templates data

=item B<GetGPSupportedHeadGroupMap>

    $SupportedHeadGroupDataRef = GetGPSupportedHeadGroupMap();

Return a reference to a hash containing supported head groups data.

=item B<GetGPTemplateID>

    $HeadGroupID = GetGPTemplateID($HeadGroupAbbrev, $ChainsAbbrev);

Return a supported template ID for compound abbreviation.

=item B<IsGPChainsAbbrevSupported>

    $Status = IsGPChainsAbbrevSupported($Abbrev, $PrintWarning);

Return 1 or 0 based on whether GP abbreviated is supported. For unsupported GP abbreviations,
a warning is printed unless PrintWarning flag is set.

=item B<ParseGPAbbrev>

    ($HeadGroup, $ChainsAbbrev, $AbbrevModifier) =
       ParseGPAbbrev($Abbrev);

Parse GP abbreviation and return these values: HeadGroup, ChainsAbbrev,
AbbrevModifier.

=item B<ProcessGPCmpdAbbrevs>

    ProcessGPCmpdAbbrevs($CmpdAbbrevsRef, $AllowArbitraryChainAbbrev,
                         $WriteSDFile, $SDFileName);

Process specified GP abbreviations to generate structures and write them out either
a SD file or simply report number of valid abbreviations.

=item B<SetupGPCmpdAbbrevTemplateDataMap>

    $AbbrevTemplateDataMapRef =
       SetupGPCmpdAbbrevTemplateDataMap($Abbrev);

Return a reference to a hash containing template data for compound abbreviation. The
template data is used to generate SD file for compound abbreviation.

=item B<ValidateGPAbbrev>

    $Status = ValidateGPAbbrev($Abbrev);

Return 1 or 0 based on whether a GP abbreviation is valid.

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
