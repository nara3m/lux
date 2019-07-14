package SPStr;
#
# File: SPStr.pm
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
use SPChainAbbrev;

our ($VERSION, @ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);

$VERSION = '2.00';
@ISA = qw(Exporter);
@EXPORT = qw();
@EXPORT_OK = qw(GenerateCmpdOntologyData GenerateCmpdOntologySDDataLines GenerateSPChainStrData GenerateSDFile GetSPTemplatesData GetSPSupportedHeadGroupMap GetSPTemplateID IsSPChainsAbbrevSupported ParseSPAbbrev ProcessSPCmpdAbbrevs SetupSPCmpdAbbrevTemplateDataMap ValidateSPAbbrev);
%EXPORT_TAGS = (
		all  => [@EXPORT, @EXPORT_OK],
		);

# Initialize SP data...
my(%SPTemplatesDataMap, %SPSupportedHeadGroupMap);
_InitializeData();

# Process specified compound abbreviations containing any wild cards
# and count the number of valid abbreviations along with optional generation
# of a SD file.
#
sub ProcessSPCmpdAbbrevs {
  my($CmpdAbbrevsRef, $AllowArbitraryChainAbbrev, $WriteSDFile, $SDFileName) = @_;

  $AllowArbitraryChainAbbrev = defined $AllowArbitraryChainAbbrev ? $AllowArbitraryChainAbbrev : 0;

  $WriteSDFile = defined $WriteSDFile ? $WriteSDFile : 0;
  if ($WriteSDFile &&  IsEmpty($SDFileName)) {
    warn "Warning: SPStr::ProcessSPCmpdAbbrevs: No SD file name specified. Ingoring structure generation.\n";
    return;
  }

  if ($WriteSDFile) {
    print "Generating SD file $SDFileName...\n";

    open SDFILE, ">$SDFileName" or die "Error: Couldn't open $SDFileName: $! \n";

    _ProcessSPCmpdAbbrevs($CmpdAbbrevsRef, $AllowArbitraryChainAbbrev, $WriteSDFile, \*SDFILE);

    close SDFILE;
  }
  else {
    _ProcessSPCmpdAbbrevs($CmpdAbbrevsRef, $AllowArbitraryChainAbbrev, $WriteSDFile);
  }
}

# Process specified compound abbreviations containing any wild cards
# and count the number of valid abbreviations along with optional generation
# of a SD file.
#
sub _ProcessSPCmpdAbbrevs {
  my($CmpdAbbrevsRef, $AllowArbitraryChainAbbrev, $WriteSDFile, $SDFileRef) = @_;
  my($Abbrev, $AbbrevCount, $HeadGroupAbbrev, $ChainsAbbrev, $NewAbbrev, $Sn1Abbrev, $Sn2Abbrev, $AllowSubstituents, $AllowRings, @AbbrevList, @ChainsAbbrevList);

  $AbbrevCount = 0;

  ABBREV: for $Abbrev (@{$CmpdAbbrevsRef}) {
    $Sn1Abbrev = "0:0"; $Sn2Abbrev = "0:0";

    if (!($Abbrev =~ /\(/ && $Abbrev =~ /\)/)) {
      warn "Warning: Ignored SP compound abbreviation $Abbrev due to incorrect format. Valid format: HeadGroupAbbrev(Chain1Abbrev/Chain2Abbrev)\n";
      next ABBREV;
    }
    if (!ValidateSPAbbrev($Abbrev)) {
      next ABBREV;
    }

    ($HeadGroupAbbrev, $ChainsAbbrev) = ParseSPAbbrev($Abbrev);

    @ChainsAbbrevList = ();
    @ChainsAbbrevList = split /\//, $ChainsAbbrev;
    if (@ChainsAbbrevList != 2) {
      warn "Warning: Ignored SP compound Abbreviation $Abbrev due to incorrect format: Must contain chain abbreviations for both sn1 and sn2 chains only.\n";
      next ABBREV;
    }
    ($Sn1Abbrev, $Sn2Abbrev) = @ChainsAbbrevList;
    ($AllowSubstituents, $AllowRings) = (1, 0);
    if (!(SPChainAbbrev::IsChainAbbrevOkay($Sn1Abbrev, 'Sn1', $AllowSubstituents, $AllowRings, $AllowArbitraryChainAbbrev) && SPChainAbbrev::IsChainAbbrevOkay($Sn2Abbrev, 'Sn2', $AllowSubstituents, $AllowRings, $AllowArbitraryChainAbbrev))) {
      warn "Warning: Ignoring SP compound abbreviation $Abbrev)\n";
      next ABBREV;
    }
    if (!IsSPChainsAbbrevSupported($ChainsAbbrev, 1)) {
      warn "Warning: Ignoring SP compound abbreviation $Abbrev\n";
      next ABBREV;
    }
    my($WildCardInHeadGroup) = ($HeadGroupAbbrev =~ /\*/ ) ? 1 : 0;

    if (!($WildCardInHeadGroup || ChainAbbrev::IsWildCardInChainAbbrev($Sn1Abbrev) || ChainAbbrev::IsWildCardInChainAbbrev($Sn2Abbrev))) {
      my($TemplateHeadGroup) = GetSPTemplateID($HeadGroupAbbrev, $ChainsAbbrev);
      if (exists $SPTemplatesDataMap{$TemplateHeadGroup}) {
	$NewAbbrev = "$HeadGroupAbbrev($ChainsAbbrev)";
	$AbbrevCount++;
	if ($WriteSDFile) {
	  _GenerateAndWriteCmdDataString($SDFileRef, $NewAbbrev);
	}
      }
      else {
	warn "Warning: Ignored SP compound abbreviation $Abbrev : Abbreviation doesn't match any template\n";
      }
      next ABBREV;
    }

    # Arbitrary acyl chain abbreviation is not supported with wild cards...
    if ($AllowArbitraryChainAbbrev) {
      warn "Warning: Ignoring SP compound abbreviation $Abbrev : Allow arbitrary chain abbreviation option is not supported with wild cards in any part of the abbreviation\n";
      next ABBREV;
    }

    # Substituents are not supported with wild cards...
    if (ChainAbbrev::IsSubstituentInChainAbbrev($Sn1Abbrev) || ChainAbbrev::IsSubstituentInChainAbbrev($Sn2Abbrev)) {
      warn "Warning: Ignoring SP compound abbreviation $Abbrev : Substituent specifications are not supported with wild cards in any part of the abbreviation\n";
      next ABBREV;
    }

    # Get expanded abbreviation for each position...
    my($Sn1ExpandedAbbrevRef) = SPChainAbbrev::ExpandChainAbbrev($Sn1Abbrev, 'Sn1');
    my($Sn2ExpandedAbbrevRef) = SPChainAbbrev::ExpandChainAbbrev($Sn2Abbrev, 'Sn2');

    if ($Sn2Abbrev =~ /^(\*|\*:\*)$/i) {
      # Add 0:0 to Sn2Abbrev list containing wild cards...
      unshift(@{$Sn2ExpandedAbbrevRef}, "0:0")
    }

    # Enumerate various possibilities...
    my($ExpandedAbbrev, $ExpandedSn1Abbrev, $ExpandedSn2Abbrev, @WildHeadGroupCmpdAbbrevs);
    @WildHeadGroupCmpdAbbrevs = ();

    if ($WildCardInHeadGroup) {
      my($SupportedHeadGroupAbbrev);
      for $SupportedHeadGroupAbbrev (sort { $a cmp $b } keys %SPSupportedHeadGroupMap ) {
	for $ExpandedSn1Abbrev (@$Sn1ExpandedAbbrevRef) {
	  for $ExpandedSn2Abbrev (@$Sn2ExpandedAbbrevRef) {
	    $ExpandedAbbrev = $SupportedHeadGroupAbbrev . '(' . $ExpandedSn1Abbrev . '/' . $ExpandedSn2Abbrev . ')';
	    $AbbrevCount++;
	    if ($WriteSDFile) {
	      _GenerateAndWriteCmdDataString($SDFileRef, $ExpandedAbbrev);
	    }
	  }
	}
      }
    }
    else {
      for $ExpandedSn1Abbrev (@$Sn1ExpandedAbbrevRef) {
	for $ExpandedSn2Abbrev (@$Sn2ExpandedAbbrevRef) {
	  $ExpandedAbbrev = $HeadGroupAbbrev . '(' . $ExpandedSn1Abbrev . '/' . $ExpandedSn2Abbrev . ')';
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
  $CmpdAbbrevTemplateDataMapRef = SetupSPCmpdAbbrevTemplateDataMap($Abbrev);

  # Generate structure data for chains...
  my($Sn1AtomBlockLinesRef, $Sn1BondBlockLinesRef) = SPStr::GenerateSPChainStrData('Sn1', $CmpdAbbrevTemplateDataMapRef);
  my($Sn2AtomBlockLinesRef, $Sn2BondBlockLinesRef) = SPStr::GenerateSPChainStrData('Sn2', $CmpdAbbrevTemplateDataMapRef);
  my($Sn3AtomBlockLinesRef, $Sn3BondBlockLinesRef) = SPStr::GenerateSPChainStrData('Sn3', $CmpdAbbrevTemplateDataMapRef);

  # Generate data block lines including various desriptors...
  my($OntologyDataMapRef) = GenerateCmpdOntologyData($CmpdAbbrevTemplateDataMapRef);
  my($OntologyDataLinesRef) = GenerateCmpdOntologySDDataLines($CmpdAbbrevTemplateDataMapRef, $OntologyDataMapRef);


  # Write out first four SD file lines: Name, MiscInfo, Comments, Count
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
sub GenerateSPChainStrData {
  my($ChainType, $CmpdAbbrevTemplateDataMapRef) = @_;

  return ChainStr::GenerateChainStrData($ChainType, $CmpdAbbrevTemplateDataMapRef);
}

# Generate ontology data lines for SD file...
sub GenerateCmpdOntologySDDataLines {
  my($CmpdAbbrevTemplateDataMapRef, $OntologyDataMapRef);

  if (@_ == 2) {
    ($CmpdAbbrevTemplateDataMapRef, $OntologyDataMapRef) = @_;
  }
  else {
    ($CmpdAbbrevTemplateDataMapRef) = @_;
    $OntologyDataMapRef = GenerateCmpdOntologyData($CmpdAbbrevTemplateDataMapRef);
  }

  my(@OntologyDataLines) = ();

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
  my($Sn1ChainIndex, $Sn2ChainIndex, $Sn1Abbrev, $Sn2Abbrev, $Sn1ChainAbbrev, $Sn2ChainAbbrev, $Sn2SubstituentsAbbrev, $Sn2ChainName, $HeadGroupAbbrevID, $HeadGroupName, $HeadGroupNameBeforeBase, $HeadGroupAbbrevName, $HeadGroupAbbrev,  $Sn2ChainAbbrevToNameMap, $SphingosineBaseAtSn1, $Sn1SubstituentsName, $Sn2SubstituentsName, $ChainLength, $DoubleBondCount, $DoubleBondGeometry, $Substituents, $SystematicName, $Sn1Name, $Sn2Name);

  $Sn1ChainIndex = 0; $Sn2ChainIndex = 1;
  $Sn1Abbrev = $CmpdAbbrevTemplateDataMapRef->{SnChainAdd}[$Sn1ChainIndex] ? $CmpdAbbrevTemplateDataMapRef->{SnAbbrev}[$Sn1ChainIndex] : '0:0';
  $Sn2Abbrev = $CmpdAbbrevTemplateDataMapRef->{SnChainAdd}[$Sn2ChainIndex] ? $CmpdAbbrevTemplateDataMapRef->{SnAbbrev}[$Sn2ChainIndex] : '0:0';

  $Sn1ChainAbbrev = $Sn1Abbrev;
  $Sn2ChainAbbrev = $Sn2Abbrev;

  $Sn1SubstituentsName = ChainAbbrev::SetupChainSubstituentsName($CmpdAbbrevTemplateDataMapRef, $Sn1ChainIndex);

  $Sn2ChainAbbrevToNameMap = SPChainAbbrev::GetChainAbbrevToNameMap('Sn2');
  $Sn2ChainName = (exists $Sn2ChainAbbrevToNameMap->{$Sn2ChainAbbrev}) ? $Sn2ChainAbbrevToNameMap->{$Sn2ChainAbbrev} : ChainAbbrev::SetupChainName($CmpdAbbrevTemplateDataMapRef, $Sn2ChainIndex);
  if ($Sn2ChainName =~ /^\(/) {
    $Sn2ChainName =~ s/^\(//;
    $Sn2ChainName =~ s/\)$//;
  };
  $Sn2SubstituentsName = ChainAbbrev::SetupChainSubstituentsName($CmpdAbbrevTemplateDataMapRef, $Sn2ChainIndex);

  $HeadGroupAbbrevID = $CmpdAbbrevTemplateDataMapRef->{AbbrevID};
  $HeadGroupName = $CmpdAbbrevTemplateDataMapRef->{HeadGroupName};
  $HeadGroupNameBeforeBase = $CmpdAbbrevTemplateDataMapRef->{HeadGroupNameBeforeBase};

  $HeadGroupAbbrevName = $CmpdAbbrevTemplateDataMapRef->{HeadGroupAbbrev};

  # Setup systematic name...
  ($SystematicName, $Sn1Name, $Sn2Name) = ('') x 3;
  ($ChainLength, $DoubleBondCount, $DoubleBondGeometry) = ChainAbbrev::ParseChainAbbrev($Sn1Abbrev);

  my($Sn1SubstituentName, $Sn1HeadGroupName, $Sn1ChainNamePrefix, $Sn1BaseName);

  $Sn1SubstituentName = '';
  $Sn1HeadGroupName = '';
  $Sn1ChainNamePrefix = '';
  $Sn1BaseName = '';

  if ($Sn1SubstituentsName) {
    $Sn1SubstituentName = "${Sn1SubstituentsName}-";
  }

  if ($ChainLength != 18) {
    my($ChainNamePrefixMapRef);
    $ChainNamePrefixMapRef = ChainAbbrev::GetChainLenToNamePrefixMap();

    if (exists $ChainNamePrefixMapRef->{$ChainLength}) {
      $Sn1ChainNamePrefix = $ChainNamePrefixMapRef->{$ChainLength} . 'a';
    }
  }

  if ($HeadGroupAbbrevID =~ /Cert$/i) {
    $Sn1HeadGroupName = '4R-hydroxy-';
  }
  elsif ($HeadGroupAbbrevID =~ /Cerm$/i) {
    $Sn1HeadGroupName = '3-keto-';
  }

  if ($DoubleBondCount == 0) {
    $Sn1BaseName = 'sphinganine';
    $Sn1Name = "${Sn1SubstituentName}${Sn1HeadGroupName}${Sn1ChainNamePrefix}${Sn1BaseName}";
  }
  elsif ($DoubleBondCount == 1) {
    if ($DoubleBondGeometry =~ /^4E$/) {
      $Sn1BaseName = 'sphing-4-enine';
    }
    else {
      $Sn1BaseName = "sphing-${DoubleBondGeometry}-enine";
    }
    $Sn1Name = "${Sn1SubstituentName}${Sn1HeadGroupName}${Sn1ChainNamePrefix}${Sn1BaseName}";
  }
  else {
    my($CountNamePrefix, $CountNamePrefixMapRef);

    $CountNamePrefixMapRef = ChainAbbrev::GetCountToNamePrefixMap();
    $CountNamePrefix = exists $CountNamePrefixMapRef->{$DoubleBondCount} ? $CountNamePrefixMapRef->{$DoubleBondCount} : '';

    $Sn1BaseName = "sphinga${CountNamePrefix}enine";

    $Sn1Name = "${Sn1SubstituentName}${DoubleBondGeometry}-${Sn1HeadGroupName}${Sn1ChainNamePrefix}${Sn1BaseName}";
  }

  if ($HeadGroupName) {
    if ($HeadGroupNameBeforeBase) {
      $Sn1Name = "${HeadGroupName}-${Sn1Name}";
    }
    else {
      $Sn1Name = "${Sn1Name}-${HeadGroupName}";
    }
  }

  if ($Sn2Abbrev !~ /^0:0$/i) {
    if (IsEmpty($Sn2SubstituentsName)) {
      $Sn2Name = "N-(${Sn2ChainName})";
    }
    else {
      $Sn2Name = "N-(${Sn2SubstituentsName}-${Sn2ChainName})";
    }
  }

  $SystematicName = $Sn1Name;
  if ($Sn2Name) {
    $SystematicName = "${Sn2Name}-${Sn1Name}";
  }

  if ($HeadGroupAbbrevID =~ /Cert$/i) {
    $Sn1Abbrev = "t${Sn1Abbrev}";
  }
  elsif ($HeadGroupAbbrevID =~ /Cerm$/i) {
    $Sn1Abbrev = "m${Sn1Abbrev}";
  }
  else {
    $Sn1Abbrev = "d${Sn1Abbrev}";
  }
  $OntologyDataMapRef->{'Systematic Name'} = "$SystematicName";

  # Set up abbreviaton...
  if ($HeadGroupAbbrevID =~ /(Cert|Cerm)$/i) {
    $HeadGroupAbbrev = $HeadGroupAbbrevID;
    $HeadGroupAbbrev =~ s/[tm]$//;
  }
  else {
    $HeadGroupAbbrev = $HeadGroupAbbrevID;
  }

  $HeadGroupAbbrev = $HeadGroupAbbrevName ? $HeadGroupAbbrevName : $HeadGroupAbbrev;
  $OntologyDataMapRef->{Abbrev} = "$HeadGroupAbbrev($Sn1Abbrev/$Sn2Abbrev)";
}


# LM classification info...
sub _SetupLMClassificationInfo {
  my($CmpdAbbrevTemplateDataMapRef, $OntologyDataMapRef) = @_;
  my($LMCategory, $LMMainClass, $LMSubClass, $SphingosineSubClass, $SphinganineSubClass, $SPTemplateID, $HeadGroupAbbrevID, $Sn1Abbrev, $Sn2Abbrev);


  $LMCategory = $CmpdAbbrevTemplateDataMapRef->{LMCategory};
  $LMMainClass = '';
  $LMSubClass = '';

  $HeadGroupAbbrevID = $CmpdAbbrevTemplateDataMapRef->{AbbrevID};
  $Sn1Abbrev = $CmpdAbbrevTemplateDataMapRef->{SnChainAdd}[0] ? $CmpdAbbrevTemplateDataMapRef->{SnAbbrev}[0] : '0:0';
  $Sn2Abbrev = $CmpdAbbrevTemplateDataMapRef->{SnChainAdd}[1] ? $CmpdAbbrevTemplateDataMapRef->{SnAbbrev}[1] : '0:0';
  $SPTemplateID = GetSPTemplateID($HeadGroupAbbrevID, "${Sn1Abbrev}/${Sn2Abbrev}");

  ($LMMainClass, $LMSubClass) = _GetLMMainAndSubClasses($SPTemplateID, $Sn1Abbrev, $Sn2Abbrev);

  $OntologyDataMapRef->{'LM Category'} = $LMCategory;
  $OntologyDataMapRef->{'LM Main Class'} = $LMCategory . $LMMainClass;
  $OntologyDataMapRef->{'LM Sub Class'} = $LMCategory . $LMMainClass . $LMSubClass;

}

# Multiple bond counts...
sub _SetupChainLengthAndMultipleBondsCountInfo {
  my($CmpdAbbrevTemplateDataMapRef, $OntologyDataMapRef) = @_;

  my($Sn1Abbrev, $Sn2Abbrev, $HeadGroupAbbrev, $Sn1ChainLength, $Sn1DoubleBonds, $Sn2ChainLength, $Sn2DoubleBonds, $ChainLengthAbbrev, $DoubleBondCountAbbrev, $DoubleBondGeometryAbbrev);

  $Sn1Abbrev = $CmpdAbbrevTemplateDataMapRef->{SnChainAdd}[0] ? $CmpdAbbrevTemplateDataMapRef->{SnAbbrev}[0] : '0:0';
  $Sn2Abbrev = $CmpdAbbrevTemplateDataMapRef->{SnChainAdd}[1] ? $CmpdAbbrevTemplateDataMapRef->{SnAbbrev}[1] : '0:0';

  $HeadGroupAbbrev = $CmpdAbbrevTemplateDataMapRef->{AbbrevID};

  # Sn1 chain length and double bond count...
  ($ChainLengthAbbrev, $DoubleBondCountAbbrev, $DoubleBondGeometryAbbrev) = ChainAbbrev::ParseChainAbbrev($Sn1Abbrev);
  $Sn1ChainLength = $ChainLengthAbbrev ? $CmpdAbbrevTemplateDataMapRef->{SnChainLength}[0] : '0';
  $Sn1DoubleBonds = $DoubleBondCountAbbrev ? $DoubleBondCountAbbrev : '0';
  $OntologyDataMapRef->{'Sn1 Chain Length'} = $Sn1ChainLength;
  $OntologyDataMapRef->{'Sn1 Double Bonds'} = $Sn1DoubleBonds;

  # Sn2 chain length and double bond count...
  ($ChainLengthAbbrev, $DoubleBondCountAbbrev, $DoubleBondGeometryAbbrev) = ChainAbbrev::ParseChainAbbrev($Sn2Abbrev);
  $Sn2ChainLength = $ChainLengthAbbrev ? $CmpdAbbrevTemplateDataMapRef->{SnChainLength}[1] : '0';
  $Sn2DoubleBonds = $DoubleBondCountAbbrev ? $DoubleBondCountAbbrev : '0';
  $OntologyDataMapRef->{'Sn2 Chain Length'} = $Sn2ChainLength;
  $OntologyDataMapRef->{'Sn2 Double Bonds'} = $Sn2DoubleBonds;
}

# Return reference to SPTemplatesDataMap...
sub GetSPTemplatesData {
  return \%SPTemplatesDataMap;
}

# Return reference to SPSupportedHeadGroupMap...
sub GetSPSupportedHeadGroupMap {
  return \%SPSupportedHeadGroupMap;
}

# Get SP template ID...
sub GetSPTemplateID {
  my($HeadGroupAbbrev, $ChainsAbbrev) = @_;
  my($HeadGroupID);

  $HeadGroupID = "";
  if ($HeadGroupAbbrev) {
    my(@AbbrevWords) = split /\//, $ChainsAbbrev;
    my ($Sn1Abbrev, $Sn2Abbrev) = @AbbrevWords;

    if ($Sn2Abbrev eq "0:0") {
      $HeadGroupID = "${HeadGroupAbbrev}Sn1";
    }
    else {
      $HeadGroupID = "${HeadGroupAbbrev}Sn1Sn2";
    }
  }
  return $HeadGroupID;
}

# Does template exist to handle this abbreviation?
#
#
sub IsSPChainsAbbrevSupported {
  my($Abbrev, $PrintWarning) = @_;
  my($Sn1Abbrev, $Sn2Abbrev, @AbbrevList);

  $Sn1Abbrev = "0:0"; $Sn2Abbrev = "0:0";
  @AbbrevList = split /\//, $Abbrev;
  if (@AbbrevList != 2) {
    if ($PrintWarning) {
      warn "Warning: Ignoring SP compound abbreviation $Abbrev : Must contain chain abbreviation for sn1 and sn2 only\n";
    }
    return 0;
  }
  if (ChainAbbrev::IsAlkenylChainAbbrev($Sn2Abbrev)) {
    if ($PrintWarning) {
      warn "Warning: Ignoring SP compound abbreviation $Abbrev : alkenyl chain at sn2 position.\n";
    }
    return 0;
  }
  if (ChainAbbrev::IsAlkylChainAbbrev($Sn2Abbrev) && !ChainAbbrev::IsAlkylChainAbbrev($Sn1Abbrev)) {
    if ($PrintWarning) {
      warn "Warning: Ignoring SP compound abbreviation $Abbrev : alkyl chain at sn2 position without alkyl chain at sn1.\n";
    }
    return 0;
  }

  return 1;
}

# Parse abbreviation...
sub ParseSPAbbrev {
  my($Abbrev, $ModifySn1Abbrev, $ModifyHeadGroupID);

  $ModifySn1Abbrev = 1;
  $ModifyHeadGroupID = 1;
  if (@_ == 3) {
    ($Abbrev, $ModifySn1Abbrev, $ModifyHeadGroupID) = @_;
  }
  elsif (@_ == 2) {
    ($Abbrev, $ModifySn1Abbrev) = @_;
  }
  else {
    ($Abbrev) = @_;
  }

  my($HeadGroup, $ChainsAbbrev, $AbbrevModifier);
  if ($Abbrev =~ /\]$/) {
    ($HeadGroup, $ChainsAbbrev, $AbbrevModifier) = $Abbrev =~ /^(.+?)\((.+?)\)\[(.+?)\]$/;
    $AbbrevModifier = '[' . $AbbrevModifier . ']';
  }
  else {
    ($HeadGroup, $ChainsAbbrev) = $Abbrev =~ /^(.+?)\((.+?)\)$/;
    $AbbrevModifier = '';
  }

  if ($ModifyHeadGroupID) {
    if ($ChainsAbbrev =~ /^[tm]/) {
      # Append the 3-keto (m) and 4-hydroxy (t) desingation with the headgroup...
      $HeadGroup .= substr($ChainsAbbrev, 0, 1)
    }
  }

  # Take out starting "d", "t", and "m" from Sn1 abbreviation...
  if ($ModifySn1Abbrev) {
    $ChainsAbbrev =~ s/^[dtm]//;
  }

  return ($HeadGroup, $ChainsAbbrev, $AbbrevModifier);
}

#
# Check out the validity of SP abbreviation...
sub ValidateSPAbbrev {
  my($Abbrev) = @_;

  if (!($Abbrev =~ /\(/ && $Abbrev =~ /\)/)) {
    warn "Warning: Ignored SP compound abbreviation $Abbrev due to incorrect format. Valid format: HeadGroupAbbrev(Chain1Abbrev/Chain2Abbrev)\n";
    return 0;
  }

  # Make sure head group is there...
  if ($Abbrev =~ /^\(/) {
    warn "Warning: Ignored SP compound abbreviation $Abbrev due to incorrect format. Valid format: HeadGroupAbbrev(Chain1Abbrev/Chain2Abbrev)\n";
    return 0;
  }

  # Number of '(' and ')' parenthesis must match...
  my(@LeftParenthesis) = $Abbrev =~ /\(/g;
  my(@RightParenthesis) = $Abbrev =~ /\)/g;
  if ($#LeftParenthesis != $#RightParenthesis) {
    warn "Warning: Ignored SP compound abbreviation $Abbrev due to incorrect format: Number of ( and ) doesn't match. Valid format: HeadGroupAbbrev(Chain1Abbrev/Chain2Abbrev).\n";
    return 0;
  }

  @LeftParenthesis = $Abbrev =~ /\[/g;
  @RightParenthesis = $Abbrev =~ /\]/g;
  if ($#LeftParenthesis != $#RightParenthesis) {
    warn "Warning: Ignored SP compound abbreviation $Abbrev due to incorrect format: Number of [ and ] doesn't match. Valid format: HeadGroupAbbrev(Chain1Abbrev/Chain2Abbrev).\n";
    return 0;
  }

  my($ModifyHeadGroupID, $ModifySn1Abbrev);
  $ModifyHeadGroupID = 0;
  $ModifySn1Abbrev = 0;
  my($HeadGroup, $ChainsAbbrev, $AbbrevModifier) = ParseSPAbbrev($Abbrev, $ModifySn1Abbrev, $ModifyHeadGroupID);

  if ($HeadGroup !~ /\*/ ) {
    if ($HeadGroup =~ /(Cert|Cerm)$/i) {
      # These are used for internal tracking only..
      warn "Warning: Ignored SP compound abbreviation $Abbrev : Unknown head group $HeadGroup.\n";
      return 0;
    }
    if (!(exists $SPSupportedHeadGroupMap{$HeadGroup})) {
      warn "Warning: Ignored SP compound abbreviation $Abbrev : Unknown head group $HeadGroup.\n";
      return 0;
    }
  }

  my(@ChainsAbbrevList) = ();
  @ChainsAbbrevList = split /\//, $ChainsAbbrev;
  if (@ChainsAbbrevList != 2) {
    warn "Warning: Ignored SP compound abbreviation $Abbrev due to incorrect format: Didn't specify abbreviations for  both sn1 and sn2 chains.\n";
    return 0;
  }

  my($Sn1Abbrev, $Sn2Abbrev) = @ChainsAbbrevList;

  if ($Sn1Abbrev !~ /\*/ ) {
    if ($Sn1Abbrev =~ /^[tm]/) {
      if ($HeadGroup =~ /^(PICer|MIPCer|MIP2Cer)$/i) {
	# Keto not allowed...
	if ($Sn1Abbrev !~ /^t/) {
	  warn "Warning: Ignored SP compound abbreviation $Abbrev :  sn1 abbreviation value, $Sn1Abbrev. not valid for headgroup $HeadGroup. Allowed value: t chain prefix only\n";
	  return 0;
	}
      }
      elsif ($HeadGroup !~ /^Cer$/i) {
	warn "Warning: Ignored SP compound abbreviation $Abbrev : Headgroup value, $HeadGroup, is not valid for sn1 abbreviation $Sn1Abbrev. Allowed values: Cer, PICer, MIPCer, MIP2Cer\n";
	return 0;
      }
    }

    if ($Sn1Abbrev =~ /^t/) {
      # 4-hydroxy...
      # Double bond at position 4 is not allowed...
      my($ChainLength, $DoubleBondCount, $DoubleBondGeometry) = ChainAbbrev::ParseChainAbbrev($Sn1Abbrev);
      if ($DoubleBondGeometry =~ /^4E/) {
	warn "Warning: Ignoring SP compound abbreviation $Abbrev : sn1 abbreviation value, $Sn1Abbrev, is not valid: Double bond specification 4E is not allowed\n";
	return 0;
      }
    }
    elsif ($Sn1Abbrev =~ /^m/) {
      # 3-keto...
      if (SPChainAbbrev::IsSphingosineChainAbbrev($Sn1Abbrev)) {
	warn "Warning: Ignoring SP compound abbreviation $Abbrev : sn1 abbreviation value, $Sn1Abbrev, is not valid. Allowed value: only sphinganine abbreviations\n";
	return 0;
      }
      if ($Sn2Abbrev ne '0:0') {
	warn "Warning: Ignoring SP compound abbreviation $Abbrev : sn2 abbreviation value, $Sn2Abbrev, is not valid. Allowed value: 0:0\n";
	return 0;
      }
    }
    elsif ($Sn1Abbrev !~ /^[dtm]/) {
      # 3-keto(t) and 4-hydroxy(m)...
      warn "Warning: Ignoring SP compound abbreviation $Abbrev : sn1 abbreviation value, $Sn1Abbrev, must start with a letter d followed by supported chain length; additionally, t18:0 and m18:0 are also supported.\n";
      return 0;
    }
  }

  if ($Sn1Abbrev eq "d0:0") {
    warn "Warning: Ignoring SP compound abbreviation $Abbrev : sn1 abbreviation value of d0:0 is not allowed.\n";
    return 0;
  }
  if (ChainAbbrev::IsAlkylChainAbbrev($Sn1Abbrev)) {
    warn "Warning: Ignoring SP compound abbreviation $Abbrev : Alkyl chain format for sn1 abbreviation, $Sn1Abbrev, is not supported.\n";
    return 0;
  }
  if (ChainAbbrev::IsAlkenylChainAbbrev($Sn1Abbrev)) {
    warn "Warning: Ignoring SP compound abbreviation $Abbrev : Alkenyl chain format for sn1 abbreviation, $Sn1Abbrev, is not supported.\n";
    return 0;
  }

  if ($Sn2Abbrev =~ /^d/) {
    warn "Warning: Ignoring SP compound abbreviation $Abbrev : sn2 abbreviation value $Sn2Abbrev, starting with a letter d is not allowed.\n";
    return 0;
  }
  if (ChainAbbrev::IsAlkylChainAbbrev($Sn2Abbrev)) {
    warn "Warning: Ignoring SP compound abbreviation $Abbrev : Alkyl chain format for sn2 abbreviation, $Sn2Abbrev, is not supported.\n";
    return 0;
  }
  if (ChainAbbrev::IsAlkenylChainAbbrev($Sn2Abbrev)) {
    warn "Warning: Ignoring SP compound abbreviation $Abbrev : Alkenyl chain format for sn2 abbreviation, $Sn1Abbrev, is not supported.\n";
    return 0;
  }

  if ($AbbrevModifier) {
    warn "Warning: Ignored SP compound abbreviation $Abbrev: stereochemistry abbreviation modifier is not supported.\n";
    return 0;
  }

  return 1;
}

# Set up template data for a specific cmpd abbreviation and return a
# reference to the hash...
sub SetupSPCmpdAbbrevTemplateDataMap {
  my($AbbrevHeadGroup, $Abbrev, $ChainsAbbrev, $AbbrevModifier, $Sn1Abbrev, $Sn2Abbrev, $Sn3Abbrev, $TemplateID, %AbbrevTemplateDataMap);

  ($Abbrev) = @_;
  %AbbrevTemplateDataMap = ();

  ($AbbrevHeadGroup, $ChainsAbbrev, $AbbrevModifier) = ParseSPAbbrev($Abbrev);

  #$ChainsAbbrev =~ s/\)//g;
  ($Sn1Abbrev, $Sn2Abbrev) = split /\//, $ChainsAbbrev;
  $Sn3Abbrev = '';

  my(@SnChainAdd) = ();
  $TemplateID = GetSPTemplateID($AbbrevHeadGroup, $ChainsAbbrev);
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

  ChainStr::SetupTemplateDataMap('SP', \%AbbrevTemplateDataMap, $SPTemplatesDataMap{$TemplateID});

  return \%AbbrevTemplateDataMap;
}

# Assign LM MainClass and SubClass...
sub _GetLMMainAndSubClasses {
  my($SPTemplateID, $Sn1Abbrev, $Sn2Abbrev) = @_;
  my($LMMainClass, $LMSubClass) = ('') x 2;

  TEMPLATEID: {
      if ($SPTemplateID =~ /^CerSn1Sn2$/i) {$LMMainClass = '02'; $LMSubClass = (SPChainAbbrev::IsSphingosineChainAbbrev($Sn1Abbrev)) ? '01' : '02'; last TEMPLATEID; };
      if ($SPTemplateID =~ /^CerSn1$/i) {$LMMainClass = '01'; $LMSubClass = (SPChainAbbrev::IsSphingosineC18ChainAbbrev($Sn1Abbrev)) ? '01' : ((SPChainAbbrev::IsSphinganineC18ChainAbbrev($Sn1Abbrev)) ? '02' : '04'); last TEMPLATEID; };

      # 4-hydroxy...
      if ($SPTemplateID =~ /^CertSn1Sn2$/i) {$LMMainClass = '02'; $LMSubClass = (SPChainAbbrev::IsSphingosineChainAbbrev($Sn1Abbrev)) ? '01' : '03'; last TEMPLATEID; };
      if ($SPTemplateID =~ /^CertSn1$/i) {$LMMainClass = '01'; $LMSubClass = (SPChainAbbrev::IsSphinganineC18ChainAbbrev($Sn1Abbrev)) ? '03' : '04'; last TEMPLATEID; };

      # 3-keto...
      if ($SPTemplateID =~ /^CermSn1Sn2$/i) {$LMMainClass = '02'; $LMSubClass = (SPChainAbbrev::IsSphingosineChainAbbrev($Sn1Abbrev)) ? '01' : '02'; last TEMPLATEID; };
      if ($SPTemplateID =~ /^CermSn1$/i) {$LMMainClass = '01'; $LMSubClass = (SPChainAbbrev::IsSphinganineC18ChainAbbrev($Sn1Abbrev)) ? '02' : '04'; last TEMPLATEID; };

      if ($SPTemplateID =~ /^CerPSn1Sn2$/i) {$LMMainClass = '02'; $LMSubClass = '05'; last TEMPLATEID; };
      if ($SPTemplateID =~ /^CerPSn1$/i) {$LMMainClass = '01'; $LMSubClass = (SPChainAbbrev::IsSphinganineC18ChainAbbrev($Sn1Abbrev) || SPChainAbbrev::IsSphingosineC18ChainAbbrev($Sn1Abbrev)) ? '05' : '04'; last TEMPLATEID; };

      if ($SPTemplateID =~ /^SMSn1Sn2$/i) {$LMMainClass = '03'; $LMSubClass = '01'; last TEMPLATEID; };
      if ($SPTemplateID =~ /^SMSn1$/i) {$LMMainClass = '01'; $LMSubClass = (SPChainAbbrev::IsSphinganineC18ChainAbbrev($Sn1Abbrev) || SPChainAbbrev::IsSphingosineC18ChainAbbrev($Sn1Abbrev)) ? '06' : '04'; last TEMPLATEID; };

      if ($SPTemplateID =~ /^PECerSn1Sn2$/i) {$LMMainClass = '03'; $LMSubClass = '02'; last TEMPLATEID; };
      if ($SPTemplateID =~ /^PECerSn1$/i) {$LMMainClass = '01'; $LMSubClass = '04'; last TEMPLATEID; };

      if ($SPTemplateID =~ /^PICerSn1Sn2$/i) {$LMMainClass = '03'; $LMSubClass = '03'; last TEMPLATEID; };
      if ($SPTemplateID =~ /^PICerSn1$/i) {$LMMainClass = '01'; $LMSubClass = '04'; last TEMPLATEID; };

      if ($SPTemplateID =~ /^PICertSn1Sn2$/i) {$LMMainClass = '03'; $LMSubClass = '03'; last TEMPLATEID; };
      if ($SPTemplateID =~ /^PICertSn1$/i) {$LMMainClass = '01'; $LMSubClass = '04'; last TEMPLATEID; };

      if ($SPTemplateID =~ /^MIPCerSn1Sn2$/i) {$LMMainClass = '03'; $LMSubClass = '03'; last TEMPLATEID; };
      if ($SPTemplateID =~ /^MIPCerSn1$/i) {$LMMainClass = '01'; $LMSubClass = '04'; last TEMPLATEID; };

      if ($SPTemplateID =~ /^MIPCertSn1Sn2$/i) {$LMMainClass = '03'; $LMSubClass = '03'; last TEMPLATEID; };
      if ($SPTemplateID =~ /^MIPCertSn1$/i) {$LMMainClass = '01'; $LMSubClass = '04'; last TEMPLATEID; };

      if ($SPTemplateID =~ /^MIP2CerSn1Sn2$/i) {$LMMainClass = '03'; $LMSubClass = '03'; last TEMPLATEID; };
      if ($SPTemplateID =~ /^MIP2CerSn1$/i) {$LMMainClass = '01'; $LMSubClass = '04'; last TEMPLATEID; };

      if ($SPTemplateID =~ /^MIP2CertSn1Sn2$/i) {$LMMainClass = '03'; $LMSubClass = '03'; last TEMPLATEID; };
      if ($SPTemplateID =~ /^MIP2CertSn1$/i) {$LMMainClass = '01'; $LMSubClass = '04'; last TEMPLATEID; };

      if ($SPTemplateID =~ /^GlcCerSn1Sn2$/i) {$LMMainClass = '05'; $LMSubClass = '01'; last TEMPLATEID; };
      if ($SPTemplateID =~ /^GlcCerSn1$/i) {$LMMainClass = '01'; $LMSubClass = (SPChainAbbrev::IsSphinganineC18ChainAbbrev($Sn1Abbrev) || SPChainAbbrev::IsSphingosineC18ChainAbbrev($Sn1Abbrev)) ? '06' : '04'; last TEMPLATEID; };

      if ($SPTemplateID =~ /^GalCerSn1Sn2$/i) {$LMMainClass = '05'; $LMSubClass = '09'; last TEMPLATEID; };
      if ($SPTemplateID =~ /^GalCerSn1$/i) {$LMMainClass = '01'; $LMSubClass = (SPChainAbbrev::IsSphinganineC18ChainAbbrev($Sn1Abbrev) || SPChainAbbrev::IsSphingosineC18ChainAbbrev($Sn1Abbrev)) ? '06' : '04'; last TEMPLATEID; };

      if ($SPTemplateID =~ /^LacCerSn1Sn2$/i) {$LMMainClass = '05'; $LMSubClass = '01'; last TEMPLATEID; };
      if ($SPTemplateID =~ /^LacCerSn1$/i) {$LMMainClass = '01'; $LMSubClass = (SPChainAbbrev::IsSphinganineC18ChainAbbrev($Sn1Abbrev) || SPChainAbbrev::IsSphingosineC18ChainAbbrev($Sn1Abbrev)) ? '06' : '04'; last TEMPLATEID; };

      if ($SPTemplateID =~ /^GM3CerSn1Sn2$/i) {$LMMainClass = '06'; $LMSubClass = '01'; last TEMPLATEID; };
      if ($SPTemplateID =~ /^GM3CerSn1$/i) {$LMMainClass = '01'; $LMSubClass = (SPChainAbbrev::IsSphinganineC18ChainAbbrev($Sn1Abbrev) || SPChainAbbrev::IsSphingosineC18ChainAbbrev($Sn1Abbrev)) ? '06' : '04'; last TEMPLATEID; };

      if ($SPTemplateID =~ /^GM4CerSn1Sn2$/i) {$LMMainClass = '06'; $LMSubClass = '01'; last TEMPLATEID; };
      if ($SPTemplateID =~ /^GM4CerSn1$/i) {$LMMainClass = '01'; $LMSubClass = (SPChainAbbrev::IsSphinganineC18ChainAbbrev($Sn1Abbrev) || SPChainAbbrev::IsSphingosineC18ChainAbbrev($Sn1Abbrev)) ? '06' : '04'; last TEMPLATEID; };

      if ($SPTemplateID =~ /^GB3CerSn1Sn2$/i) {$LMMainClass = '05'; $LMSubClass = '02'; last TEMPLATEID; };
      if ($SPTemplateID =~ /^GB3CerSn1$/i) {$LMMainClass = '01'; $LMSubClass = (SPChainAbbrev::IsSphinganineC18ChainAbbrev($Sn1Abbrev) || SPChainAbbrev::IsSphingosineC18ChainAbbrev($Sn1Abbrev)) ? '06' : '04'; last TEMPLATEID; };

      if ($SPTemplateID =~ /^iGB3CerSn1Sn2$/i) {$LMMainClass = '05'; $LMSubClass = '06'; last TEMPLATEID; };
      if ($SPTemplateID =~ /^iGB3CerSn1$/i) {$LMMainClass = '01'; $LMSubClass = (SPChainAbbrev::IsSphinganineC18ChainAbbrev($Sn1Abbrev) || SPChainAbbrev::IsSphingosineC18ChainAbbrev($Sn1Abbrev)) ? '06' : '04'; last TEMPLATEID; };

      if ($SPTemplateID =~ /^asialo-GM2CerSn1Sn2$/i) {$LMMainClass = '05'; $LMSubClass = '03'; last TEMPLATEID; };
      if ($SPTemplateID =~ /^asialo-GM2CerSn1$/i) {$LMMainClass = '01'; $LMSubClass = (SPChainAbbrev::IsSphinganineC18ChainAbbrev($Sn1Abbrev) || SPChainAbbrev::IsSphingosineC18ChainAbbrev($Sn1Abbrev)) ? '06' : '04'; last TEMPLATEID; };

      if ($SPTemplateID =~ /^Lc3CerSn1Sn2$/i) {$LMMainClass = '05'; $LMSubClass = '04'; last TEMPLATEID; };
      if ($SPTemplateID =~ /^Lc3CerSn1$/i) {$LMMainClass = '01'; $LMSubClass = (SPChainAbbrev::IsSphinganineC18ChainAbbrev($Sn1Abbrev) || SPChainAbbrev::IsSphingosineC18ChainAbbrev($Sn1Abbrev)) ? '06' : '04'; last TEMPLATEID; };

      if ($SPTemplateID =~ /^Manb1-4GlcCerSn1Sn2$/i) {$LMMainClass = '05'; $LMSubClass = '01'; last TEMPLATEID; };
      if ($SPTemplateID =~ /^Manb1-4GlcCerSn1$/i) {$LMMainClass = '01'; $LMSubClass = (SPChainAbbrev::IsSphinganineC18ChainAbbrev($Sn1Abbrev) || SPChainAbbrev::IsSphingosineC18ChainAbbrev($Sn1Abbrev)) ? '06' : '04'; last TEMPLATEID; };

      if ($SPTemplateID =~ /^MolluCerSn1Sn2$/i) {$LMMainClass = '05'; $LMSubClass = '07'; last TEMPLATEID; };
      if ($SPTemplateID =~ /^MolluCerSn1$/i) {$LMMainClass = ''; $LMSubClass = (SPChainAbbrev::IsSphinganineC18ChainAbbrev($Sn1Abbrev) || SPChainAbbrev::IsSphingosineC18ChainAbbrev($Sn1Abbrev)) ? '06' : '04'; last TEMPLATEID; };

      if ($SPTemplateID =~ /^ArthroCerSn1Sn2$/i) {$LMMainClass = '05'; $LMSubClass = '08'; last TEMPLATEID; };
      if ($SPTemplateID =~ /^ArthroCerSn1$/i) {$LMMainClass = '01'; $LMSubClass = (SPChainAbbrev::IsSphinganineC18ChainAbbrev($Sn1Abbrev) || SPChainAbbrev::IsSphingosineC18ChainAbbrev($Sn1Abbrev)) ? '06' : '04'; last TEMPLATEID; };

      ($LMMainClass, $LMSubClass) = ('') x 2;
  }
  return ($LMMainClass, $LMSubClass);
}

# Initialize SP data...
sub _InitializeData {
  _InitializeStrTemplates();
  _InitializeSupportedHeadGroupsData();
}

# Initialize structure template data for these supported templates:
#
sub _InitializeStrTemplates {
  %SPTemplatesDataMap = ();

  # Cer templates...
  my($CerSn1Sn2TemplateString)=<<ENDTEMPLATESTRING;
Cer sn1 acyl and sn2 acyl template structure
  LipdMAPS02060709152D

 11 10  0  0  0  0  0  0  0  0999 V2000
    0.2232    0.4310    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4911    0.8423    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2056    0.4310    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6361   -0.2833    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1897   -0.2833    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.9378    0.8435    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.6523    0.4310    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9380   -0.7055    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9380   -1.5316    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1022    1.5158    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8892    1.5316    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  1  5  1  6      
  1  4  1  1      
  6  1  1  0      
  7  6  1  0      
  8  9  2  0      
  8  5  1  0      
  2 10  1  1      
  2 11  1  6      
M  END

ENDTEMPLATESTRING

  my($CerSn1TemplateString)=<<ENDTEMPLATESTRING;
Cer sn1 acyl template structure
  LipdMAPS02070609152D

  9  8  0  0  0  0  0  0  0  0999 V2000
   -0.0001   -0.1931    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7144    0.2181    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4289   -0.1931    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4127   -0.9074    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4130   -0.9074    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.7144    0.2193    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4289   -0.1931    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3255    0.8916    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1125    0.9074    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  1  5  1  6      
  1  4  1  1      
  6  1  1  0      
  7  6  1  0      
  2  8  1  1      
  2  9  1  6      
M  END

ENDTEMPLATESTRING

  # Cert (4-hydroxy) templates...
  my($CertSn1Sn2TemplateString)=<<ENDTEMPLATESTRING;
Cert sn1 acyl and sn2 acyl template structure
  LipdMAPS02060709152D

 14 13  0  0  0  0  0  0  0  0999 V2000
    0.7199    0.3747    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0056    0.7859    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7089    0.3747    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4289    0.7859    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1492    0.3747    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.1327   -0.3396    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.3070   -0.3396    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    1.4344    0.7871    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.1492    0.3747    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.3945    1.4594    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3925    1.4752    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7050   -0.1369    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5670   -0.8615    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5670   -1.4752    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  1  0      
  1  7  1  6      
  1  6  1  1      
  8  1  1  0      
  9  8  1  0      
  2 10  1  1      
  2 11  1  6      
  3 12  1  6      
  7 13  1  0      
 13 14  2  0      
M  END

ENDTEMPLATESTRING

  my($CertSn1TemplateString)=<<ENDTEMPLATESTRING;
Cert sn1 acyl template structure
  LipdMAPS02070609152D

 12 11  0  0  0  0  0  0  0  0999 V2000
    0.7199   -0.1931    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0056    0.2181    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7089   -0.1931    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4289    0.2181    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1491   -0.1931    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.1327   -0.9074    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.3070   -0.9074    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    1.4344    0.2193    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.1491   -0.1931    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.3945    0.8916    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3925    0.9074    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7050   -0.7460    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  1  0      
  1  7  1  6      
  1  6  1  1      
  8  1  1  0      
  9  8  1  0      
  2 10  1  1      
  2 11  1  6      
  3 12  1  6      
M  END

ENDTEMPLATESTRING

  # Cerm (3-keto) templates...
  my($CermSn1Sn2TemplateString)=<<ENDTEMPLATESTRING;
Cerm sn1 acyl and sn2 acyl template structure
  LipdMAPS02060709152D

 11 10  0  0  0  0  0  0  0  0999 V2000
    0.2395    0.3963    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4748    0.8075    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1893    0.3963    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6523   -0.3181    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1734   -0.3181    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.9540    0.8087    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.6686    0.3963    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4748    1.4134    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9477   -0.7652    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6686   -0.3490    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9477   -1.4134    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  1  5  1  6      
  1  4  1  1      
  6  1  1  0      
  7  6  1  0      
  2  8  2  0      
  5  9  1  0      
  9 10  1  0      
  9 11  2  0      
M  END

ENDTEMPLATESTRING

  my($CermSn1TemplateString)=<<ENDTEMPLATESTRING;
Cerm sn1 acyl template structure
  LipdMAPS02070609152D

  8  7  0  0  0  0  0  0  0  0999 V2000
   -0.0001   -0.1513    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7144    0.2599    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4289   -0.1513    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4127   -0.8657    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4130   -0.8657    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.7144    0.2611    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4289   -0.1513    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7144    0.8657    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  1  5  1  6      
  1  4  1  1      
  6  1  1  0      
  7  6  1  0      
  2  8  2  0      
M  END

ENDTEMPLATESTRING

  # PCer templates...
  my($CerPSn1Sn2TemplateString)=<<ENDTEMPLATESTRING;
Cer sn1 acyl and sn2 acyl template structure
  LipdMAPS02060709152D

 15 14  0  0  0  0  0  0  0  0999 V2000
   -0.6611    0.4310    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3755    0.8423    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0899    0.4310    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2482   -0.2833    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0741   -0.2833    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.0534    0.8435    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7679    0.4310    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.5366    0.4138    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.7845    0.7200    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    1.4205    0.0892    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.7845    1.4723    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8223   -0.7055    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8223   -1.5316    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9866    1.5158    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7735    1.5316    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
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
 12  5  1  0      
  9  7  1  0      
  2 14  1  1      
  2 15  1  6      
M  END

ENDTEMPLATESTRING

  my($CerPSn1TemplateString)=<<ENDTEMPLATESTRING;
Cer sn1 acyl template structure
  LipdMAPS02070609152D

 13 12  0  0  0  0  0  0  0  0999 V2000
   -0.8845   -0.1932    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5989    0.2182    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.3134   -0.1932    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4716   -0.9075    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2975   -0.9075    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1700    0.2194    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5446   -0.1932    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.3134   -0.2104    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.5612    0.0959    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    1.1972   -0.5350    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.5612    0.8482    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2100    0.8917    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9970    0.9075    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
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
  2 12  1  1      
  2 13  1  6      
M  END

ENDTEMPLATESTRING

  # PECer templates...
  my($PECerSn1Sn2TemplateString)=<<ENDTEMPLATESTRING;
PECer sn1 acyl and sn2 acyl template structure
  LipdMAPS02060709152D

 18 17  0  0  0  0  0  0  0  0999 V2000
   -1.7328    0.4310    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4473    0.8423    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1617    0.4310    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3199   -0.2833    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1459   -0.2833    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0183    0.8435    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3038    0.4310    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.4650    0.4138    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.7129    0.7200    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    0.3489    0.0892    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.7129    1.4724    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.8941   -0.7055    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.8941   -1.5317    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0584    1.5159    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.8453    1.5317    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.1795    0.0013    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.8940    0.4138    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.6085    0.0013    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
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
 12  5  1  0      
  9  7  1  0      
  2 14  1  1      
  2 15  1  6      
 16  8  1  0      
 17 16  1  0      
 18 17  1  0      
M  END

ENDTEMPLATESTRING

  my($PECerSn1TemplateString)=<<ENDTEMPLATESTRING;
PECer sn1 acyl template structure
  LipdMAPS02070609152D

 16 15  0  0  0  0  0  0  0  0999 V2000
   -1.9561   -0.1932    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6706    0.2181    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.3850   -0.1932    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5433   -0.9075    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.3692   -0.9075    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2417    0.2193    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5272   -0.1932    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.2416   -0.2104    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.4895    0.0958    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    0.1255   -0.5350    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.4895    0.8482    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2817    0.8917    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0686    0.9075    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.9560   -0.6229    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.6705   -0.2104    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.3850   -0.6229    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
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
  2 12  1  1      
  2 13  1  6      
 14  8  1  0      
 15 14  1  0      
 16 15  1  0      
M  END

ENDTEMPLATESTRING

  # PICer templates...
  my($PICerSn1Sn2TemplateString)=<<ENDTEMPLATESTRING;
PICer sn1 acyl and sn2 acyl template structure
  LipdMAPS02060709152D

 26 26  0  0  0  0  0  0  0  0999 V2000
   -3.0591    0.4310    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.7735    0.8423    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.4880    0.4310    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6462   -0.2833    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -3.4721   -0.2833    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -2.3445    0.8435    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6300    0.4310    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.1388    0.4138    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6133    0.7200    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9774    0.0892    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6133    1.4724    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.2204   -0.7055    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.2204   -1.5317    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.3846    1.5159    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.1715    1.5317    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.9337    0.7404    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.6228    1.0906    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.9420   -0.0941    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.2544    0.2593    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.5686   -0.0941    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.2462    1.0906    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.8261    0.8772    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.2504    1.3049    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.4147    1.3372    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.4819    0.2053    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.9347    0.9168    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
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
 12  5  1  0      
  9  7  1  0      
  2 14  1  1      
  2 15  1  6      
 16 17  1  0      
 16 21  1  0      
 17 18  1  0      
 20 21  1  0      
 20 19  1  1      
 18 19  1  1      
 17 22  1  0      
 19 23  1  0      
 16 24  1  0      
 20 25  1  0      
 21 26  1  0      
  8 18  1  0      
M  END

ENDTEMPLATESTRING

  # PICerSn1 templates...
  my($PICerSn1TemplateString)=<<ENDTEMPLATESTRING;
PICer sn1 acyl template structure
  LipdMAPS02060709152D

 24 24  0  0  0  0  0  0  0  0999 V2000
   -3.2877   -0.1933    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.0026    0.2183    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.7175   -0.1933    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.8746   -0.9081    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -3.7010   -0.9081    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5728    0.2195    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8577   -0.1933    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0878   -0.2105    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8405    0.0960    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2047   -0.5353    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8405    0.8487    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.6134    0.8923    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.4009    0.9081    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.7168    0.1085    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4060    0.4587    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7253   -0.7259    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0376   -0.3725    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.3516   -0.7259    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.0291    0.4587    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6094    0.2453    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.0336    0.6730    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.1977    0.7053    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.2648   -0.4265    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.7175    0.2849    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
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
  2 12  1  1      
  2 13  1  6      
 14 15  1  0      
 14 19  1  0      
 15 16  1  0      
 18 19  1  0      
 18 17  1  1      
 16 17  1  1      
 15 20  1  0      
 17 21  1  0      
 14 22  1  0      
 18 23  1  0      
 19 24  1  0      
 16  8  1  0      
M  END

ENDTEMPLATESTRING

  # PICert templates...
  my($PICertSn1Sn2TemplateString)=<<ENDTEMPLATESTRING;
PICert sn1 acyl and sn2 acyl template structure
  LipdMAPS02060709152D

 29 29  0  0  0  0  0  0  0  0999 V2000
   -2.5627    0.3747    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.2770    0.7859    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.9916    0.3747    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.7116    0.7859    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.4319    0.3747    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1499   -0.3396    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.9756   -0.3396    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8482    0.7871    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1334    0.3747    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.8881    1.4594    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.6751    1.4752    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -3.9877   -0.1369    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.8497   -0.8615    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.8497   -1.4752    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.6306    0.3596    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1215    0.6658    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4854    0.0350    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1215    1.4179    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.4329    0.6783    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.1232    1.0282    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4430   -0.1554    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.7542    0.1977    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.0671   -0.1554    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.7441    1.0282    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.3272    0.8150    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.7502    1.2423    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.9134    1.2746    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.9796    0.1437    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.4319    0.8546    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  1  0      
  1  7  1  6      
  1  6  1  1      
  8  1  1  0      
  9  8  1  0      
  2 10  1  1      
  2 11  1  6      
  3 12  1  6      
  7 13  1  0      
 13 14  2  0      
 16 15  1  0      
 16 17  1  0      
 16 18  2  0      
 19 20  1  0      
 19 24  1  0      
 20 21  1  0      
 23 24  1  0      
 23 22  1  1      
 21 22  1  1      
 20 25  1  0      
 22 26  1  0      
 19 27  1  0      
 23 28  1  0      
 24 29  1  0      
 21 15  1  0      
  9 16  1  0      
M  END

ENDTEMPLATESTRING

  my($PICertSn1TemplateString)=<<ENDTEMPLATESTRING;
PICert sn1 acyl template structure
  LipdMAPS02060709152D

 27 27  0  0  0  0  0  0  0  0999 V2000
   -2.5608   -0.1931    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.2752    0.2181    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.9897   -0.1931    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.7097    0.2181    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.4299   -0.1931    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1480   -0.9074    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.9737   -0.9074    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8463    0.2193    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1316   -0.1931    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.8862    0.8916    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.6733    0.9074    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -3.9858   -0.7460    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.6289   -0.2103    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1231    0.0959    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4870   -0.5349    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1231    0.8479    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.4310    0.1084    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.1214    0.4582    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4413   -0.7253    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.7524   -0.3722    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.0652   -0.7253    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.7422    0.4582    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.3255    0.2451    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.7484    0.6723    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.9115    0.7046    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.9776   -0.4262    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.4299    0.2847    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  1  0      
  1  7  1  6      
  1  6  1  1      
  8  1  1  0      
  9  8  1  0      
  2 10  1  1      
  2 11  1  6      
  3 12  1  6      
 14 13  1  0      
 14 15  1  0      
 14 16  2  0      
 17 18  1  0      
 17 22  1  0      
 18 19  1  0      
 21 22  1  0      
 21 20  1  1      
 19 20  1  1      
 18 23  1  0      
 20 24  1  0      
 17 25  1  0      
 21 26  1  0      
 22 27  1  0      
 19 13  1  0      
  9 14  1  0      
M  END

ENDTEMPLATESTRING

  # SM templates...
  my($SMSn1Sn2TemplateString)=<<ENDTEMPLATESTRING;
SM sn1 acyl and sn2 acyl template structure
  LipdMAPS02060709152D

 21 20  0  0  0  0  0  0  0  0999 V2000
   -2.0900    0.4310    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.8043    0.8422    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.5187    0.4310    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6771   -0.2833    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5029   -0.2833    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3755    0.8435    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6610    0.4310    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.1076    0.4138    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.8219    0.0013    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.5364    0.4138    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.2508    0.0013    0.0000 N   0  3  0  0  0  0  0  0  0  0  0  0
    3.9654    0.4138    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.2508   -0.8237    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.9654   -0.4112    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3555    0.7200    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0085    0.0892    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.3555    1.4722    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.2511   -0.7055    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.2511   -1.5315    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4154    1.5157    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.2023    1.5315    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
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
 18  5  1  0      
 15  7  1  0      
  2 20  1  1      
  2 21  1  6      
M  CHG  2  11   1  16  -1
M  END

ENDTEMPLATESTRING

  my($SMSn1TemplateString)=<<ENDTEMPLATESTRING;
SM sn1 acyl template structure
  LipdMAPS02070609152D

 19 18  0  0  0  0  0  0  0  0999 V2000
   -2.3133    0.0771    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0276    0.4883    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.7421    0.0771    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9005   -0.6372    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.7262   -0.6372    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5988    0.4896    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8844    0.0771    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.8842    0.0599    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.5985   -0.3526    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.3130    0.0599    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0275   -0.3526    0.0000 N   0  3  0  0  0  0  0  0  0  0  0  0
    3.7421    0.0599    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0275   -1.1776    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.7421   -0.7651    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.1322    0.3661    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2319   -0.2647    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.1322    1.1183    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6388    1.1618    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.4257    1.1776    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
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
  2 18  1  1      
  2 19  1  6      
M  CHG  2  11   1  16  -1
M  END

ENDTEMPLATESTRING

  # GlcCer templates...
  my($GlcCerSn1Sn2TemplateString)=<<ENDTEMPLATESTRING;
GlcCer sn1 acyl and sn2 acyl template structure
  LipdMAPS02060709152D

 22 22  0  0  0  0  0  0  0  0999 V2000
   -2.2212    0.4310    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.9355    0.8423    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.6500    0.4310    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8083   -0.2833    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6341   -0.2833    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5067    0.8435    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7922    0.4310    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.3824   -0.7055    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.3824   -1.5316    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5466    1.5158    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.3336    1.5316    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.4173    1.1519    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.7496   -0.0047    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4574    0.3338    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.1732   -0.0335    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.8410    1.1233    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.1331    0.7850    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.0967    0.9700    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.0842    0.9207    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.4120    0.1727    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.8294   -0.0683    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.7165    1.3683    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  1  5  1  6      
  1  4  1  1      
  6  1  1  0      
  7  6  1  0      
  8  9  2  0      
  8  5  1  0      
  2 10  1  1      
  2 11  1  6      
 15 16  1  0      
 16 17  1  0      
 17 12  1  0      
 12 18  1  0      
 12 13  1  1      
 15 14  1  1      
 13 14  1  1      
 16 19  1  0      
 13 20  1  0      
 20 21  1  0      
 17 22  1  0      
  7 15  1  0      
M  END

ENDTEMPLATESTRING

  my($GlcCerSn1TemplateString)=<<ENDTEMPLATESTRING;
GlcCer sn1 acyl template structure
  LipdMAPS02070609152D

 20 20  0  0  0  0  0  0  0  0999 V2000
   -2.4446   -0.1931    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1589    0.2181    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.8733   -0.1931    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0317   -0.9074    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.8575   -0.9074    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7300    0.2193    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0155   -0.1931    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.7700    0.8916    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.5569    0.9074    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.1940    0.5278    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.5263   -0.6289    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2341   -0.2903    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0501   -0.6576    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6176    0.4992    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.9098    0.1609    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.8733    0.3458    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1392    0.2965    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.1886   -0.4514    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.6060   -0.6924    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.4932    0.7442    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  1  5  1  6      
  1  4  1  1      
  6  1  1  0      
  7  6  1  0      
  2  8  1  1      
  2  9  1  6      
 13 14  1  0      
 14 15  1  0      
 15 10  1  0      
 10 16  1  0      
 10 11  1  1      
 13 12  1  1      
 11 12  1  1      
 14 17  1  0      
 11 18  1  0      
 18 19  1  0      
 15 20  1  0      
  7 13  1  0      
M  END

ENDTEMPLATESTRING

  # GalCer templates...
  my($GalCerSn1Sn2TemplateString)=<<ENDTEMPLATESTRING;
GalCer sn1 acyl and sn2 acyl template structure
  LipdMAPS02060709152D

 22 22  0  0  0  0  0  0  0  0999 V2000
   -2.0875    0.2693    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.8018    0.6806    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.5162    0.2693    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6746   -0.4450    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5004   -0.4450    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3730    0.6818    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6585    0.2693    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.2487   -0.8672    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.2487   -1.6932    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4129    1.3541    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1999    1.3699    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.5508    0.9902    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.8832   -0.1664    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5910    0.1721    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.3068   -0.1952    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.9746    0.9616    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.2667    0.6233    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.2178    0.7590    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.5455    0.0110    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.9629   -0.2300    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.8501    1.2066    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.5508    1.6932    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  1  5  1  6      
  1  4  1  1      
  6  1  1  0      
  7  6  1  0      
  8  9  2  0      
  8  5  1  0      
  2 10  1  1      
  2 11  1  6      
 15 16  1  0      
 16 17  1  0      
 17 12  1  0      
 12 13  1  1      
 15 14  1  1      
 13 14  1  1      
 16 18  1  0      
 13 19  1  0      
 19 20  1  0      
 17 21  1  0      
  7 15  1  0      
 12 22  1  0      
M  END

ENDTEMPLATESTRING

  my($GalCerSn1TemplateString)=<<ENDTEMPLATESTRING;
GalCer sn1 acyl template structure
  LipdMAPS02070609152D

 20 20  0  0  0  0  0  0  0  0999 V2000
   -2.3108   -0.3060    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0251    0.1051    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.7394   -0.3060    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8979   -1.0203    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.7237   -1.0203    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5963    0.1063    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8818   -0.3060    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6362    0.7786    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.4230    0.7944    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.3274    0.4148    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.6598   -0.7418    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.3677   -0.4032    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.0835   -0.7705    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7512    0.3862    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0433    0.0479    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0055    0.1835    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.3220   -0.5643    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.7394   -0.8053    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.6267    0.6312    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.3274    1.0203    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  1  5  1  6      
  1  4  1  1      
  6  1  1  0      
  7  6  1  0      
  2  8  1  1      
  2  9  1  6      
 13 14  1  0      
 14 15  1  0      
 15 10  1  0      
 10 11  1  1      
 13 12  1  1      
 11 12  1  1      
 14 16  1  0      
 11 17  1  0      
 17 18  1  0      
 15 19  1  0      
  7 13  1  0      
 10 20  1  0      
M  END

ENDTEMPLATESTRING

  # LacCer templates...
  my($LacCerSn1Sn2TemplateString)=<<ENDTEMPLATESTRING;
LacCer sn1 acyl and sn2 acyl template structure
  LipdMAPS02060709152D

 33 34  0  0  0  0  0  0  0  0999 V2000
   -4.5828   -0.0598    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.2972    0.3515    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.0116   -0.0598    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.1700   -0.7740    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -4.9958   -0.7740    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -3.8684    0.3527    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1540   -0.0598    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.7441   -1.1962    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.7441   -2.0222    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.9083    1.0249    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.6953    1.0407    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.0551    0.6611    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3875   -0.4955    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9046   -0.1570    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1887   -0.5243    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5209    0.6325    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2289    0.2942    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.7345    0.4792    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2777    0.4299    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.0498   -0.3181    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4672   -0.5591    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.3544    0.8774    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    6.0462    1.1972    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.3785    0.0408    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.0863    0.3793    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.8023    0.0120    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.4700    1.1686    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.7620    0.8304    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.7133    0.9660    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    6.0409    0.2182    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.4583   -0.0228    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.3454    1.4136    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    6.0462    2.0222    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  1  5  1  6      
  1  4  1  1      
  6  1  1  0      
  7  6  1  0      
  8  9  2  0      
  8  5  1  0      
  2 10  1  1      
  2 11  1  6      
 15 16  1  0      
 16 17  1  0      
 17 12  1  0      
 12 18  1  0      
 12 13  1  1      
 15 14  1  1      
 13 14  1  1      
 16 19  1  0      
 13 20  1  0      
 20 21  1  0      
 17 22  1  0      
  7 15  1  0      
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
 18 26  1  0      
M  END

ENDTEMPLATESTRING

  my($LacCerSn1TemplateString)=<<ENDTEMPLATESTRING;
LacCer sn1 acyl template structure
  LipdMAPS02070609152D

 31 32  0  0  0  0  0  0  0  0999 V2000
   -4.8061   -0.6839    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.5206   -0.2726    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.2350   -0.6839    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.3933   -1.3981    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -5.2191   -1.3981    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -4.0918   -0.2714    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.3773   -0.6839    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.1316    0.4008    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.9186    0.4166    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.8318    0.0370    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.1642   -1.1196    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1279   -0.7811    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4120   -1.1484    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7442    0.0084    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4522   -0.3299    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5112   -0.1449    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5010   -0.1942    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.8265   -0.9422    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2439   -1.1832    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.1311    0.2533    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.8229    0.5731    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.1552   -0.5833    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.8630   -0.2448    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.5790   -0.6121    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.2467    0.5445    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.5387    0.2063    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.4900    0.3419    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.8176   -0.4059    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.2350   -0.6469    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.1221    0.7895    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.8229    1.3981    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  1  5  1  6      
  1  4  1  1      
  6  1  1  0      
  7  6  1  0      
  2  8  1  1      
  2  9  1  6      
 13 14  1  0      
 14 15  1  0      
 15 10  1  0      
 10 16  1  0      
 10 11  1  1      
 13 12  1  1      
 11 12  1  1      
 14 17  1  0      
 11 18  1  0      
 18 19  1  0      
 15 20  1  0      
  7 13  1  0      
 24 25  1  0      
 25 26  1  0      
 26 21  1  0      
 21 22  1  1      
 24 23  1  1      
 22 23  1  1      
 25 27  1  0      
 22 28  1  0      
 28 29  1  0      
 26 30  1  0      
 21 31  1  0      
 16 24  1  0      
M  END

ENDTEMPLATESTRING

  # GM3Cer templates...
  my($GM3CerSn1Sn2TemplateString)=<<ENDTEMPLATESTRING;
GM3Cer sn1 acyl and sn2 acyl template structure
  LipdMAPS02060709152D

 53 55  0  0  0  0  0  0  0  0999 V2000
   -7.0194   -1.2312    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.7337   -0.8199    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -8.4481   -1.2312    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.6066   -1.9454    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -7.4323   -1.9454    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -6.3050   -0.8187    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.5906   -1.2312    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -8.1806   -2.3676    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -8.1806   -3.1936    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -7.3448   -0.1464    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -8.1318   -0.1306    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3814   -0.5103    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0490   -1.6669    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.3411   -1.3284    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.6252   -1.6957    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.9574   -0.5389    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6654   -0.8772    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7020   -0.6922    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.7142   -0.7415    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3867   -1.4895    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9693   -1.7305    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0821   -0.2939    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.6096    0.0259    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.9419   -1.1306    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.6498   -0.7921    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.3658   -1.1594    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0335   -0.0027    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.3255   -0.3409    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.2768   -0.2053    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.6043   -0.9532    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.0217   -1.1942    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.9088    0.2423    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.6096    0.8509    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.9821    2.6373    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.7026    1.6819    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.6937    2.3516    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.8768    2.1650    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.5213    1.8658    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.4113    1.3958    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.4669    1.9926    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    8.1966    1.6077    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.8948    2.0472    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.2206    0.9823    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    6.6775    2.4793    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    6.0222    1.4019    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.7429    0.8723    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.0618    0.3657    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.7825   -0.1639    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    6.5825    1.2745    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.1781    0.6976    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.6937    2.9001    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.2020    3.1936    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.2373    3.1636    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  1  5  1  6      
  1  4  1  1      
  6  1  1  0      
  7  6  1  0      
  8  9  2  0      
  8  5  1  0      
  2 10  1  1      
  2 11  1  6      
 15 16  1  0      
 16 17  1  0      
 17 12  1  0      
 12 18  1  0      
 12 13  1  1      
 15 14  1  1      
 13 14  1  1      
 16 19  1  0      
 13 20  1  0      
 20 21  1  0      
 17 22  1  0      
  7 15  1  0      
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
 18 26  1  0      
 34 35  1  0      
 35 38  1  1      
 36 39  1  1      
 36 37  1  0      
 34 37  1  0      
 38 39  1  1      
 35 40  1  0      
 40 41  1  0      
 41 42  1  0      
 41 43  2  0      
 34 44  1  0      
 38 45  1  0      
 46 45  1  0      
 47 46  1  0      
 48 47  1  0      
 45 49  1  6      
 46 50  1  6      
 32 36  1  0      
 36 51  1  0      
 51 52  1  0      
 51 53  2  0      
M  END

ENDTEMPLATESTRING

  my($GM3CerSn1TemplateString)=<<ENDTEMPLATESTRING;
GM3Cer sn1 acyl template structure
  LipdMAPS02070609152D

 51 53  0  0  0  0  0  0  0  0999 V2000
   -7.2431   -1.8554    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.9574   -1.4441    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -8.6718   -1.8554    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.8302   -2.5696    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -7.6560   -2.5696    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -6.5286   -1.4429    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.8142   -1.8554    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -7.5685   -0.7705    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -8.3555   -0.7547    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6048   -1.1344    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2724   -2.2911    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.5646   -1.9526    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.8488   -2.3199    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.1809   -1.1630    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.8889   -1.5014    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9254   -1.3164    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.9378   -1.3657    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6101   -2.1137    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1927   -2.3547    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.3055   -0.9180    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.3864   -0.5982    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.7187   -1.7548    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4265   -1.4163    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.1425   -1.7836    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.8102   -0.6268    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.1022   -0.9650    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0535   -0.8294    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.3811   -1.5774    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.7985   -1.8184    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.6856   -0.3818    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.3864    0.2268    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.7590    2.0133    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.4795    1.0578    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.4705    1.7276    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.6536    1.5410    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.2982    1.2418    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.1881    0.7717    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.2439    1.3686    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    7.9736    0.9836    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.6718    1.4232    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.9976    0.3582    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    6.4544    1.8553    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.7991    0.7778    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.5198    0.2482    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.8387   -0.2584    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.5594   -0.7880    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    6.3594    0.6504    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.9550    0.0735    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.4705    2.2761    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.9788    2.5696    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.0141    2.5396    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  1  5  1  6      
  1  4  1  1      
  6  1  1  0      
  7  6  1  0      
  2  8  1  1      
  2  9  1  6      
 13 14  1  0      
 14 15  1  0      
 15 10  1  0      
 10 16  1  0      
 10 11  1  1      
 13 12  1  1      
 11 12  1  1      
 14 17  1  0      
 11 18  1  0      
 18 19  1  0      
 15 20  1  0      
  7 13  1  0      
 24 25  1  0      
 25 26  1  0      
 26 21  1  0      
 21 22  1  1      
 24 23  1  1      
 22 23  1  1      
 25 27  1  0      
 22 28  1  0      
 28 29  1  0      
 26 30  1  0      
 21 31  1  0      
 16 24  1  0      
 32 33  1  0      
 33 36  1  1      
 34 37  1  1      
 34 35  1  0      
 32 35  1  0      
 36 37  1  1      
 33 38  1  0      
 38 39  1  0      
 39 40  1  0      
 39 41  2  0      
 32 42  1  0      
 36 43  1  0      
 44 43  1  0      
 45 44  1  0      
 46 45  1  0      
 43 47  1  6      
 44 48  1  6      
 30 34  1  0      
 34 49  1  0      
 49 50  1  0      
 49 51  2  0      
M  END

ENDTEMPLATESTRING

  # GM4Cer templates...
  my($GM4CerSn1Sn2TemplateString)=<<ENDTEMPLATESTRING;
GM4Cer sn1 acyl and sn2 acyl template structure
  LipdMAPS02060709152D

 42 43  0  0  0  0  0  0  0  0999 V2000
   -4.4831   -0.8704    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.1975   -0.4590    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.9120   -0.8704    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.0703   -1.5846    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -4.8961   -1.5846    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -3.7687   -0.4578    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0542   -0.8704    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.6445   -2.0070    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.6445   -2.8330    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.8086    0.2145    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.5957    0.2303    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.1555   -0.1494    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4878   -1.3061    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8045   -0.9676    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0887   -1.3349    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4208   -0.1780    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1286   -0.5163    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1777   -0.3806    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.1502   -1.1287    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5676   -1.3697    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.4547    0.0670    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.4458    2.2767    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.1663    1.3211    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.1571    1.9909    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.3403    1.8043    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.9849    1.5050    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.8747    1.0350    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.9307    1.6318    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    5.6605    1.2469    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.3588    1.6864    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.6845    0.6214    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.1412    2.1187    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.4859    1.0411    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.2065    0.5114    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.5255    0.0048    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.2461   -0.5248    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.0462    0.9137    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.6416    0.3367    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.1571    2.5395    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.6654    2.8330    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.7007    2.8030    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.1555    0.5286    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  1  5  1  6      
  1  4  1  1      
  6  1  1  0      
  7  6  1  0      
  8  9  2  0      
  8  5  1  0      
  2 10  1  1      
  2 11  1  6      
 15 16  1  0      
 16 17  1  0      
 17 12  1  0      
 12 13  1  1      
 15 14  1  1      
 13 14  1  1      
 16 18  1  0      
 13 19  1  0      
 19 20  1  0      
 17 21  1  0      
  7 15  1  0      
 22 23  1  0      
 23 26  1  1      
 24 27  1  1      
 24 25  1  0      
 22 25  1  0      
 26 27  1  1      
 23 28  1  0      
 28 29  1  0      
 29 30  1  0      
 29 31  2  0      
 22 32  1  0      
 26 33  1  0      
 34 33  1  0      
 35 34  1  0      
 36 35  1  0      
 33 37  1  6      
 34 38  1  6      
 24 39  1  0      
 39 40  1  0      
 39 41  2  0      
 21 24  1  0      
 12 42  1  0      
M  END

ENDTEMPLATESTRING

  my($GM4CerSn1TemplateString)=<<ENDTEMPLATESTRING;
GM4Cer sn1 acyl template structure
  LipdMAPS02070609152D

 40 41  0  0  0  0  0  0  0  0999 V2000
   -4.7064   -1.4946    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.4208   -1.0832    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.1353   -1.4946    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.2936   -2.2088    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -5.1194   -2.2088    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -3.9920   -1.0820    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.2775   -1.4946    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.0319   -0.4097    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.8190   -0.3939    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.9321   -0.7736    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.2644   -1.9303    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0279   -1.5918    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.3121   -1.9591    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6442   -0.8022    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3520   -1.1405    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4011   -1.0048    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.9268   -1.7529    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.3442   -1.9939    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.2313   -0.5572    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.2223    1.6525    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.9428    0.6969    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.9337    1.3667    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.1169    1.1801    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.7614    0.8808    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.6513    0.4108    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.7072    1.0076    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    5.4370    0.6227    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.1353    1.0622    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.4610   -0.0028    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.9177    1.4945    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.2624    0.4169    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.9830   -0.1128    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.3020   -0.6194    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0226   -1.1490    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.8227    0.2895    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.4182   -0.2875    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.9337    1.9153    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4420    2.2088    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.4773    2.1788    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.9321   -0.0956    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  1  5  1  6      
  1  4  1  1      
  6  1  1  0      
  7  6  1  0      
  2  8  1  1      
  2  9  1  6      
 13 14  1  0      
 14 15  1  0      
 15 10  1  0      
 10 11  1  1      
 13 12  1  1      
 11 12  1  1      
 14 16  1  0      
 11 17  1  0      
 17 18  1  0      
 15 19  1  0      
  7 13  1  0      
 20 21  1  0      
 21 24  1  1      
 22 25  1  1      
 22 23  1  0      
 20 23  1  0      
 24 25  1  1      
 21 26  1  0      
 26 27  1  0      
 27 28  1  0      
 27 29  2  0      
 20 30  1  0      
 24 31  1  0      
 32 31  1  0      
 33 32  1  0      
 34 33  1  0      
 31 35  1  6      
 32 36  1  6      
 22 37  1  0      
 37 38  1  0      
 37 39  2  0      
 19 22  1  0      
 10 40  1  0      
M  END

ENDTEMPLATESTRING

  # GB3Cer templates...
  my($GB3CerSn1Sn2TemplateString)=<<ENDTEMPLATESTRING;
GB3Cer sn1 acyl and sn2 acyl template structure
  LipdMAPS02060709152D

 44 46  0  0  0  0  0  0  0  0999 V2000
   -6.2225   -1.5242    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.9369   -1.1129    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.6513   -1.5242    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.8097   -2.2384    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -6.6355   -2.2384    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -5.5081   -1.1117    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.7937   -1.5242    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -7.3838   -2.6606    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.3838   -3.4866    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -6.5480   -0.4395    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -7.3350   -0.4237    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5846   -0.8033    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2522   -1.9599    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5443   -1.6214    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.8284   -1.9887    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1606   -0.8319    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8686   -1.1702    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0948   -0.9852    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.9174   -1.0345    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5899   -1.7825    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1725   -2.0235    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2853   -0.5870    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.4065   -0.2672    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.7388   -1.4236    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.4466   -1.0851    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.1626   -1.4524    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.8303   -0.2958    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.1223   -0.6340    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0736   -0.4984    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.4012   -1.2462    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.8186   -1.4872    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.4065    0.5578    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.6859    2.6616    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.0182    1.5052    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.7260    1.8437    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.4420    1.4764    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.1097    2.6330    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.4017    2.2948    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.3530    2.4304    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.6806    1.6826    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.0980    1.4416    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    6.9851    2.8780    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.6859    3.4866    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.7108   -0.0455    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  1  5  1  6      
  1  4  1  1      
  6  1  1  0      
  7  6  1  0      
  8  9  2  0      
  8  5  1  0      
  2 10  1  1      
  2 11  1  6      
 15 16  1  0      
 16 17  1  0      
 17 12  1  0      
 12 18  1  0      
 12 13  1  1      
 15 14  1  1      
 13 14  1  1      
 16 19  1  0      
 13 20  1  0      
 20 21  1  0      
 17 22  1  0      
  7 15  1  0      
 26 27  1  0      
 27 28  1  0      
 28 23  1  0      
 23 24  1  1      
 26 25  1  1      
 24 25  1  1      
 27 29  1  0      
 24 30  1  0      
 30 31  1  0      
 23 32  1  0      
 18 26  1  0      
 36 37  1  0      
 37 38  1  0      
 38 33  1  0      
 33 34  1  1      
 36 35  1  1      
 34 35  1  1      
 37 39  1  0      
 34 40  1  0      
 40 41  1  0      
 38 42  1  0      
 33 43  1  0      
 28 44  1  0      
 32 36  1  0      
M  END

ENDTEMPLATESTRING

  my($GB3CerSn1TemplateString)=<<ENDTEMPLATESTRING;
GB3Cer sn1 acyl template structure
  LipdMAPS02070609152D

 42 44  0  0  0  0  0  0  0  0999 V2000
   -6.4380   -2.1457    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.1516   -1.7349    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.8651   -2.1457    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.0257   -2.8590    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -6.8505   -2.8590    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -5.7245   -1.7337    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.0110   -2.1457    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -6.7631   -1.0623    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -7.5492   -1.0465    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8070   -1.4257    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4738   -2.5809    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.7643   -2.2428    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.0468   -2.6096    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.3799   -1.4542    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0894   -1.7921    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1284   -1.6074    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.1357   -1.6566    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8123   -2.4037    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3954   -2.6444    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5068   -1.2096    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.1781   -0.8902    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.5112   -2.0452    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.2206   -1.7071    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.9381   -2.0740    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.6050   -0.9188    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.8954   -1.2566    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.8492   -1.1211    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.1728   -1.8680    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.5897   -2.1087    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.1781   -0.0662    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.4535    2.0350    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.7866    0.8800    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.4960    1.2181    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.2135    0.8513    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.8804    2.0065    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.1709    1.6687    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.1246    1.8041    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.4482    1.0572    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.8651    0.8165    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    6.7536    2.2512    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.4535    2.8590    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.4832   -0.6688    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  1  5  1  6      
  1  4  1  1      
  6  1  1  0      
  7  6  1  0      
  2  8  1  1      
  2  9  1  6      
 13 14  1  0      
 14 15  1  0      
 15 10  1  0      
 10 16  1  0      
 10 11  1  1      
 13 12  1  1      
 11 12  1  1      
 14 17  1  0      
 11 18  1  0      
 18 19  1  0      
 15 20  1  0      
  7 13  1  0      
 24 25  1  0      
 25 26  1  0      
 26 21  1  0      
 21 22  1  1      
 24 23  1  1      
 22 23  1  1      
 25 27  1  0      
 22 28  1  0      
 28 29  1  0      
 21 30  1  0      
 16 24  1  0      
 34 35  1  0      
 35 36  1  0      
 36 31  1  0      
 31 32  1  1      
 34 33  1  1      
 32 33  1  1      
 35 37  1  0      
 32 38  1  0      
 38 39  1  0      
 36 40  1  0      
 31 41  1  0      
 26 42  1  0      
 30 34  1  0      
M  END

ENDTEMPLATESTRING

  # iGB3Cer templates...
  my($iGB3CerSn1Sn2TemplateString)=<<ENDTEMPLATESTRING;
iGB3Cer sn1 acyl and sn2 acyl template structure
  LipdMAPS02060709152D

 44 46  0  0  0  0  0  0  0  0999 V2000
   -5.9440   -1.3489    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.6584   -0.9376    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.3728   -1.3489    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.5312   -2.0631    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -6.3570   -2.0631    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -5.2296   -0.9364    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.5152   -1.3489    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -7.1053   -2.4853    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.1053   -3.3113    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -6.2695   -0.2642    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -7.0565   -0.2484    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3062   -0.6280    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9738   -1.7846    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2659   -1.4461    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.5500   -1.8134    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.8821   -0.6566    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5902   -0.9949    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3732   -0.8099    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.6389   -0.8592    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3115   -1.6072    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.1059   -1.8482    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0069   -0.4117    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.6849   -0.0919    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.0172   -1.2483    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.7250   -0.9098    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.4410   -1.2771    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.1087   -0.1205    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.4007   -0.4587    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.3520   -0.3231    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.6796   -1.0709    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.0970   -1.3119    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.6849    0.7331    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.4074    2.4863    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.7397    1.3299    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.4475    1.6684    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.1635    1.3011    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.8312    2.4577    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.1232    2.1195    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.0745    2.2551    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.4021    1.5073    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.8195    1.2663    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    6.7066    2.7027    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.4074    3.3113    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.1383    0.4139    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  1  5  1  6      
  1  4  1  1      
  6  1  1  0      
  7  6  1  0      
  8  9  2  0      
  8  5  1  0      
  2 10  1  1      
  2 11  1  6      
 15 16  1  0      
 16 17  1  0      
 17 12  1  0      
 12 18  1  0      
 12 13  1  1      
 15 14  1  1      
 13 14  1  1      
 16 19  1  0      
 13 20  1  0      
 20 21  1  0      
 17 22  1  0      
  7 15  1  0      
 26 27  1  0      
 27 28  1  0      
 28 23  1  0      
 23 24  1  1      
 26 25  1  1      
 24 25  1  1      
 27 29  1  0      
 24 30  1  0      
 30 31  1  0      
 23 32  1  0      
 18 26  1  0      
 36 37  1  0      
 37 38  1  0      
 38 33  1  0      
 33 34  1  1      
 36 35  1  1      
 34 35  1  1      
 37 39  1  0      
 34 40  1  0      
 40 41  1  0      
 38 42  1  0      
 33 43  1  0      
 36 44  1  0      
 28 44  1  0      
M  END

ENDTEMPLATESTRING

  my($iGB3CerSn1TemplateString)=<<ENDTEMPLATESTRING;
iGB3Cer sn1 acyl template structure
  LipdMAPS02070609152D

 42 44  0  0  0  0  0  0  0  0999 V2000
   -6.1599   -1.9706    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.8734   -1.5598    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.5870   -1.9706    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.7476   -2.6839    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -6.5724   -2.6839    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -5.4463   -1.5586    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.7328   -1.9706    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -6.4850   -0.8872    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -7.2710   -0.8714    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5289   -1.2506    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1957   -2.4058    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4862   -2.0677    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.7688   -2.4345    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1017   -1.2791    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8114   -1.6170    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.1497   -1.4323    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.8576   -1.4815    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5342   -2.2286    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1173   -2.4693    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2288   -1.0345    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.4561   -0.7151    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.7893   -1.8701    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.4986   -1.5320    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.2162   -1.8989    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.8831   -0.7437    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.1735   -1.0815    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.1273   -0.9461    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.4509   -1.6929    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.8678   -1.9337    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.4561    0.1089    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.1754    1.8599    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.5085    0.7049    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.2178    1.0430    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.9354    0.6762    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.6023    1.8314    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.8927    1.4936    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.8465    1.6290    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.1701    0.8821    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.5870    0.6414    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    6.4754    2.0761    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.1754    2.6839    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.9102   -0.2099    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  1  5  1  6      
  1  4  1  1      
  6  1  1  0      
  7  6  1  0      
  2  8  1  1      
  2  9  1  6      
 13 14  1  0      
 14 15  1  0      
 15 10  1  0      
 10 16  1  0      
 10 11  1  1      
 13 12  1  1      
 11 12  1  1      
 14 17  1  0      
 11 18  1  0      
 18 19  1  0      
 15 20  1  0      
  7 13  1  0      
 24 25  1  0      
 25 26  1  0      
 26 21  1  0      
 21 22  1  1      
 24 23  1  1      
 22 23  1  1      
 25 27  1  0      
 22 28  1  0      
 28 29  1  0      
 21 30  1  0      
 16 24  1  0      
 34 35  1  0      
 35 36  1  0      
 36 31  1  0      
 31 32  1  1      
 34 33  1  1      
 32 33  1  1      
 35 37  1  0      
 32 38  1  0      
 38 39  1  0      
 36 40  1  0      
 31 41  1  0      
 34 42  1  0      
 26 42  1  0      
M  END

ENDTEMPLATESTRING

  # asialo-GM2Cer templates...
  my($asialoGM2CerSn1Sn2TemplateString)=<<ENDTEMPLATESTRING;
asialo-GM2Cer sn1 acyl and sn2 acyl template structure
  LipdMAPS02060709152D

 47 49  0  0  0  0  0  0  0  0999 V2000
   -6.2150   -1.5224    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.9285   -1.1116    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.6420   -1.5224    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.8027   -2.2357    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -6.6275   -2.2357    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -5.5014   -1.1104    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.7879   -1.5224    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -7.3749   -2.6574    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.3749   -3.4824    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -6.5401   -0.4390    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -7.3261   -0.4232    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5839   -0.8023    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2507   -1.9575    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5412   -1.6194    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.8238   -1.9863    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1568   -0.8309    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8663   -1.1688    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0947   -0.9840    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.9127   -1.0332    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5892   -1.7803    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1723   -2.0211    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2837   -0.5863    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.4012   -0.2669    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.7343   -1.4219    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.4436   -1.0838    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.1612   -1.4506    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.8281   -0.2954    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.1185   -0.6332    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0723   -0.4978    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.3959   -1.2447    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.8128   -1.4854    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.4012    0.5571    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.6766    2.6584    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.0097    1.5034    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.7191    1.8415    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.4366    1.4746    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.1035    2.6298    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.3940    2.2920    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.3477    2.4275    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    7.6713    1.6806    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.0882    1.4399    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    6.9766    2.8745    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.6766    3.4824    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.7063   -0.0454    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.5967    2.8610    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.9480    2.4865    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.5967    3.3667    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  1  5  1  6      
  1  4  1  1      
  6  1  1  0      
  7  6  1  0      
  8  9  2  0      
  8  5  1  0      
  2 10  1  1      
  2 11  1  6      
 15 16  1  0      
 16 17  1  0      
 17 12  1  0      
 12 18  1  0      
 12 13  1  1      
 15 14  1  1      
 13 14  1  1      
 16 19  1  0      
 13 20  1  0      
 20 21  1  0      
 17 22  1  0      
  7 15  1  0      
 26 27  1  0      
 27 28  1  0      
 28 23  1  0      
 23 24  1  1      
 26 25  1  1      
 24 25  1  1      
 27 29  1  0      
 24 30  1  0      
 30 31  1  0      
 23 32  1  0      
 18 26  1  0      
 36 37  1  0      
 37 38  1  0      
 38 33  1  0      
 33 34  1  1      
 36 35  1  1      
 34 35  1  1      
 37 39  1  0      
 34 40  1  0      
 40 41  1  0      
 38 42  1  0      
 33 43  1  0      
 28 44  1  0      
 32 36  1  0      
 39 45  1  0      
 45 46  1  0      
 45 47  2  0      
M  END

ENDTEMPLATESTRING

  my($asialoGM2CerSn1TemplateString)=<<ENDTEMPLATESTRING;
asialo-GM2Cer sn1 acyl template structure
  LipdMAPS02070609152D

 45 47  0  0  0  0  0  0  0  0999 V2000
   -6.4381   -2.1457    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.1516   -1.7349    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.8651   -2.1457    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.0258   -2.8590    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -6.8506   -2.8590    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -5.7245   -1.7337    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.0110   -2.1457    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -6.7632   -1.0623    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -7.5492   -1.0465    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8070   -1.4256    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4738   -2.5808    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.7643   -2.2427    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.0469   -2.6096    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.3799   -1.4542    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0894   -1.7921    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1284   -1.6073    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.1358   -1.6565    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8123   -2.4036    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3954   -2.6444    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5068   -1.2096    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.1781   -0.8902    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.5112   -2.0452    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.2205   -1.7071    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.9381   -2.0739    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.6050   -0.9187    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.8954   -1.2565    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.8492   -1.1211    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.1728   -1.8680    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.5897   -2.1088    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.1781   -0.0662    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.4535    2.0351    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.7866    0.8801    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.4960    1.2182    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.2135    0.8513    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.8804    2.0065    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.1709    1.6687    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.1246    1.8042    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    7.4482    1.0573    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.8651    0.8166    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    6.7535    2.2512    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.4535    2.8590    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.4832   -0.6687    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.3736    2.2377    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.7249    1.8631    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.3736    2.7434    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  1  5  1  6      
  1  4  1  1      
  6  1  1  0      
  7  6  1  0      
  2  8  1  1      
  2  9  1  6      
 13 14  1  0      
 14 15  1  0      
 15 10  1  0      
 10 16  1  0      
 10 11  1  1      
 13 12  1  1      
 11 12  1  1      
 14 17  1  0      
 11 18  1  0      
 18 19  1  0      
 15 20  1  0      
  7 13  1  0      
 24 25  1  0      
 25 26  1  0      
 26 21  1  0      
 21 22  1  1      
 24 23  1  1      
 22 23  1  1      
 25 27  1  0      
 22 28  1  0      
 28 29  1  0      
 21 30  1  0      
 16 24  1  0      
 34 35  1  0      
 35 36  1  0      
 36 31  1  0      
 31 32  1  1      
 34 33  1  1      
 32 33  1  1      
 35 37  1  0      
 32 38  1  0      
 38 39  1  0      
 36 40  1  0      
 31 41  1  0      
 26 42  1  0      
 30 34  1  0      
 37 43  1  0      
 43 44  1  0      
 43 45  2  0      
M  END

ENDTEMPLATESTRING

  # Lc3Cer templates...
  my($Lc3CerSn1Sn2TemplateString)=<<ENDTEMPLATESTRING;
Lc3Cer sn1 acyl and sn2 acyl template structure
  LipdMAPS02060709152D

 47 49  0  0  0  0  0  0  0  0999 V2000
   -6.4276   -1.4645    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.1412   -1.0537    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.8547   -1.4645    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.0153   -2.1779    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -6.8401   -2.1779    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -5.7141   -1.0525    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.0006   -1.4645    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -7.5875   -2.5996    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.5875   -3.4246    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -6.7528   -0.3811    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -7.5388   -0.3654    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7966   -0.7445    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4634   -1.8997    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.7539   -1.5616    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.0364   -1.9285    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.3694   -0.7731    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0790   -1.1110    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1180   -0.9262    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.1253   -0.9754    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8019   -1.7225    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3850   -1.9632    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4964   -0.5285    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.1885   -0.2091    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.5216   -1.3641    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.2310   -1.0260    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.9485   -1.3928    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.6154   -0.2376    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.9058   -0.5754    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.8596   -0.4400    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.1832   -1.1869    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.6001   -1.4276    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.1885    0.6149    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.4639    2.7162    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.7970    1.5612    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.5064    1.8993    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.2240    1.5324    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.8908    2.6876    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.1813    2.3498    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.1351    2.4853    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    7.4586    1.7384    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.8755    1.4977    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    6.7640    2.9323    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.4936    0.0124    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.3841    2.9189    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.7354    2.5443    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.3841    3.4246    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    8.3009    2.4919    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  1  5  1  6      
  1  4  1  1      
  6  1  1  0      
  7  6  1  0      
  8  9  2  0      
  8  5  1  0      
  2 10  1  1      
  2 11  1  6      
 15 16  1  0      
 16 17  1  0      
 17 12  1  0      
 12 18  1  0      
 12 13  1  1      
 15 14  1  1      
 13 14  1  1      
 16 19  1  0      
 13 20  1  0      
 20 21  1  0      
 17 22  1  0      
  7 15  1  0      
 26 27  1  0      
 27 28  1  0      
 28 23  1  0      
 23 24  1  1      
 26 25  1  1      
 24 25  1  1      
 27 29  1  0      
 24 30  1  0      
 30 31  1  0      
 23 32  1  0      
 18 26  1  0      
 36 37  1  0      
 37 38  1  0      
 38 33  1  0      
 33 34  1  1      
 36 35  1  1      
 34 35  1  1      
 37 39  1  0      
 34 40  1  0      
 40 41  1  0      
 38 42  1  0      
 28 43  1  0      
 32 36  1  0      
 39 44  1  0      
 44 45  1  0      
 44 46  2  0      
 33 47  1  0      
M  END

ENDTEMPLATESTRING

  my($Lc3CerSn1TemplateString)=<<ENDTEMPLATESTRING;
Lc3Cer sn1 acyl template structure
  LipdMAPS02070609152D

 45 47  0  0  0  0  0  0  0  0999 V2000
   -6.6005   -2.0721    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.3087   -1.6644    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -8.0168   -2.0721    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.1913   -2.7801    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -7.0099   -2.7801    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -5.8924   -1.6632    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.1843   -2.0721    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -6.9232   -0.9969    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -7.7033   -0.9813    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0120   -1.3575    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6738   -2.5040    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.9545   -2.1685    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.2273   -2.5326    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.5654   -1.3859    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2847   -1.7213    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3385   -1.5379    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.3156   -1.5867    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0173   -2.3281    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6035   -2.5670    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7065   -1.1432    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.9355   -0.8262    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.2736   -1.9724    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.9927   -1.6369    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.7199   -2.0009    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.3818   -0.8545    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.6624   -1.1897    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6317   -1.0553    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.9302   -1.7966    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.3440   -2.0355    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.9355   -0.0084    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.1861    2.0771    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.5243    0.9308    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.2434    1.2663    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.9707    0.9022    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.6325    2.0487    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.9132    1.7134    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.8825    1.8479    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    7.1809    1.1066    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.5946    0.8678    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    6.4915    2.2915    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.2458   -0.6063    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.1371    2.2782    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.4933    1.9064    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.1371    2.7801    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    8.0168    1.8544    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  1  5  1  6      
  1  4  1  1      
  6  1  1  0      
  7  6  1  0      
  2  8  1  1      
  2  9  1  6      
 13 14  1  0      
 14 15  1  0      
 15 10  1  0      
 10 16  1  0      
 10 11  1  1      
 13 12  1  1      
 11 12  1  1      
 14 17  1  0      
 11 18  1  0      
 18 19  1  0      
 15 20  1  0      
  7 13  1  0      
 24 25  1  0      
 25 26  1  0      
 26 21  1  0      
 21 22  1  1      
 24 23  1  1      
 22 23  1  1      
 25 27  1  0      
 22 28  1  0      
 28 29  1  0      
 21 30  1  0      
 16 24  1  0      
 34 35  1  0      
 35 36  1  0      
 36 31  1  0      
 31 32  1  1      
 34 33  1  1      
 32 33  1  1      
 35 37  1  0      
 32 38  1  0      
 38 39  1  0      
 36 40  1  0      
 26 41  1  0      
 30 34  1  0      
 37 42  1  0      
 42 43  1  0      
 42 44  2  0      
 31 45  1  0      
M  END

ENDTEMPLATESTRING

  # Manb1-4GlcCer templates...
  my($Manb14GlcCerSn1Sn2TemplateString)=<<ENDTEMPLATESTRING;
Manb1-4GlcCer sn1 acyl and sn2 acyl template structure
  LipdMAPS02060709152D

 33 34  0  0  0  0  0  0  0  0999 V2000
   -4.7850    0.0431    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.4994    0.4544    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.2138    0.0431    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.3722   -0.6711    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -5.1980   -0.6711    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -4.0706    0.4556    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.3562    0.0431    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.9463   -1.0933    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.9463   -1.9193    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.1105    1.1278    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.8975    1.1436    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.8528    0.7640    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.1852   -0.3926    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1069   -0.0541    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.3910   -0.4214    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7232    0.7354    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4312    0.3971    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5322    0.5821    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4800    0.5328    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.8475   -0.2152    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2649   -0.4562    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.1521    0.9803    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.8438    1.3001    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.1761    0.1437    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.8839    0.4822    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.6000    0.1149    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.2676    1.2715    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.5596    0.9333    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.8385    0.3211    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.2559    0.0801    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.1430    1.5165    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    6.6605    1.0813    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.2676    1.9193    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  1  5  1  6      
  1  4  1  1      
  6  1  1  0      
  7  6  1  0      
  8  9  2  0      
  8  5  1  0      
  2 10  1  1      
  2 11  1  6      
 15 16  1  0      
 16 17  1  0      
 17 12  1  0      
 12 18  1  0      
 12 13  1  1      
 15 14  1  1      
 13 14  1  1      
 16 19  1  0      
 13 20  1  0      
 20 21  1  0      
 17 22  1  0      
  7 15  1  0      
 26 27  1  0      
 27 28  1  0      
 28 23  1  0      
 23 24  1  1      
 26 25  1  1      
 24 25  1  1      
 24 29  1  0      
 29 30  1  0      
 28 31  1  0      
 18 26  1  0      
 23 32  1  0      
 27 33  1  0      
M  END

ENDTEMPLATESTRING

  my($Manb14GlcCerSn1TemplateString)=<<ENDTEMPLATESTRING;
Manb1-4GlcCer sn1 acyl template structure
  LipdMAPS02070609152D

 31 32  0  0  0  0  0  0  0  0999 V2000
   -5.0083   -0.5810    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.7226   -0.1697    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.4370   -0.5810    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.5955   -1.2952    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -5.4212   -1.2952    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -4.2939   -0.1685    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.5795   -0.5810    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.3337    0.5037    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -6.1207    0.5195    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.6294    0.1399    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0381   -1.0167    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3302   -0.6782    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6143   -1.0455    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9465    0.1113    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6545   -0.2270    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.3088   -0.0420    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.7033   -0.0913    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.6241   -0.8393    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0415   -1.0803    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0712    0.3562    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.6203    0.6760    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.9527   -0.4804    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.6605   -0.1419    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.3766   -0.5092    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0442    0.6474    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.3362    0.3092    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.6150   -0.3030    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.0324   -0.5440    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.9196    0.8924    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    6.4370    0.4572    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.0442    1.2952    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  1  5  1  6      
  1  4  1  1      
  6  1  1  0      
  7  6  1  0      
  2  8  1  1      
  2  9  1  6      
 13 14  1  0      
 14 15  1  0      
 15 10  1  0      
 10 16  1  0      
 10 11  1  1      
 13 12  1  1      
 11 12  1  1      
 14 17  1  0      
 11 18  1  0      
 18 19  1  0      
 15 20  1  0      
  7 13  1  0      
 24 25  1  0      
 25 26  1  0      
 26 21  1  0      
 21 22  1  1      
 24 23  1  1      
 22 23  1  1      
 22 27  1  0      
 27 28  1  0      
 26 29  1  0      
 16 24  1  0      
 21 30  1  0      
 25 31  1  0      
M  END

ENDTEMPLATESTRING

  # MolluCer templates...
  my($MolluCerSn1Sn2TemplateString)=<<ENDTEMPLATESTRING;
MolluCer sn1 acyl and sn2 acyl template structure
  LipdMAPS02060709152D

 44 46  0  0  0  0  0  0  0  0999 V2000
   -6.0637   -1.1428    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.7780   -0.7315    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.4924   -1.1428    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.6509   -1.8570    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -6.4767   -1.8570    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -5.3493   -0.7303    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.6349   -1.1428    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -7.2249   -2.2792    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.2249   -3.1052    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -6.3892   -0.0582    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -7.1761   -0.0424    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4260   -0.4220    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0936   -1.5785    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.3856   -1.2400    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.6697   -1.6073    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0019   -0.4506    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7099   -0.7888    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.2534   -0.6038    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.7587   -0.6531    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4313   -1.4011    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0139   -1.6421    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1267   -0.2057    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.5649    0.1141    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.8973   -1.0422    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.6051   -0.7037    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.3212   -1.0710    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.9888    0.0855    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.2808   -0.2527    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.5596   -0.8648    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.9770   -1.1058    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.8642    0.3305    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.3816   -0.1047    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.9888    0.7333    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.1224    2.4860    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.4548    1.3296    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.1626    1.6681    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.8787    1.3008    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.5463    2.4574    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.8383    2.1192    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.1171    1.5070    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.5345    1.2660    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    6.4217    2.7024    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.9391    2.2672    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.5463    3.1052    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  1  5  1  6      
  1  4  1  1      
  6  1  1  0      
  7  6  1  0      
  8  9  2  0      
  8  5  1  0      
  2 10  1  1      
  2 11  1  6      
 15 16  1  0      
 16 17  1  0      
 17 12  1  0      
 12 18  1  0      
 12 13  1  1      
 15 14  1  1      
 13 14  1  1      
 16 19  1  0      
 13 20  1  0      
 20 21  1  0      
 17 22  1  0      
  7 15  1  0      
 26 27  1  0      
 27 28  1  0      
 28 23  1  0      
 23 24  1  1      
 26 25  1  1      
 24 25  1  1      
 24 29  1  0      
 29 30  1  0      
 28 31  1  0      
 18 26  1  0      
 23 32  1  0      
 27 33  1  0      
 37 38  1  0      
 38 39  1  0      
 39 34  1  0      
 34 35  1  1      
 37 36  1  1      
 35 36  1  1      
 35 40  1  0      
 40 41  1  0      
 39 42  1  0      
 34 43  1  0      
 38 44  1  0      
 37 31  1  0      
M  END

ENDTEMPLATESTRING

  my($MolluCerSn1TemplateString)=<<ENDTEMPLATESTRING;
MolluCer sn1 acyl template structure
  LipdMAPS02070609152D

 42 44  0  0  0  0  0  0  0  0999 V2000
   -6.1346   -1.7241    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.8316   -1.3227    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.5287   -1.7241    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.7318   -2.4209    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -6.5376   -2.4209    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -5.4375   -1.3216    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.7405   -1.7241    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -6.4522   -0.6658    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -7.2200   -0.6503    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6336   -1.0207    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2850   -2.1492    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5457   -1.8189    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.7987   -2.1773    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1470   -1.0486    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8864   -1.3786    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0293   -1.1981    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.8855   -1.2462    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6388   -1.9761    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2315   -2.2113    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3173   -0.8097    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.2363   -0.4976    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.5849   -1.6259    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.3240   -1.2956    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.0712   -1.6540    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.7226   -0.5255    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.9833   -0.8555    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.2311   -1.4528    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.6384   -1.6880    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.5526   -0.2865    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.0332   -0.7111    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.7226    0.1066    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    6.7318    1.8168    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.0804    0.6884    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.8195    1.0187    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.5667    0.6603    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.2181    1.7888    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.4788    1.4588    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.7266    0.8615    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.1339    0.6263    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    6.0481    2.0279    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.5287    1.6033    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.2181    2.4209    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  1  5  1  6      
  1  4  1  1      
  6  1  1  0      
  7  6  1  0      
  2  8  1  1      
  2  9  1  6      
 13 14  1  0      
 14 15  1  0      
 15 10  1  0      
 10 16  1  0      
 10 11  1  1      
 13 12  1  1      
 11 12  1  1      
 14 17  1  0      
 11 18  1  0      
 18 19  1  0      
 15 20  1  0      
  7 13  1  0      
 24 25  1  0      
 25 26  1  0      
 26 21  1  0      
 21 22  1  1      
 24 23  1  1      
 22 23  1  1      
 22 27  1  0      
 27 28  1  0      
 26 29  1  0      
 16 24  1  0      
 21 30  1  0      
 25 31  1  0      
 35 36  1  0      
 36 37  1  0      
 37 32  1  0      
 32 33  1  1      
 35 34  1  1      
 33 34  1  1      
 33 38  1  0      
 38 39  1  0      
 37 40  1  0      
 32 41  1  0      
 36 42  1  0      
 35 29  1  0      
M  END

ENDTEMPLATESTRING

  # ArthroCer templates...
  my($ArthroCerSn1Sn2TemplateString)=<<ENDTEMPLATESTRING;
ArthroCer sn1 acyl and sn2 acyl template structure
  LipdMAPS02060709152D

 47 49  0  0  0  0  0  0  0  0999 V2000
   -7.2599   -0.1893    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.9743    0.2219    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -8.6887   -0.1893    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.8471   -0.9035    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -7.6729   -0.9035    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -6.5455    0.2231    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.8311   -0.1893    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -8.4212   -1.3257    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -8.4212   -2.1517    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -7.5854    0.8953    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -8.3724    0.9111    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6222    0.5315    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2898   -0.6250    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.5819   -0.2865    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.8660   -0.6538    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.1982    0.5029    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.9062    0.1646    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9428    0.3496    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.9550    0.3003    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6275   -0.4476    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2101   -0.6886    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.3229    0.7478    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.3687    1.0676    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.7010   -0.0887    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4088    0.2497    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.1250   -0.1175    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7925    1.0390    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0845    0.7008    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.3634    0.0887    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.7808   -0.1523    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.6679    1.2840    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.1854    0.8488    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.7925    1.6868    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    8.3187    1.5214    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.6510    0.3650    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.3588    0.7035    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.0750    0.3362    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.7425    1.4928    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.0345    1.1546    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.3134    0.5424    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.7308    0.3014    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.6179    1.7378    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    9.1354    1.3026    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.9846    1.2897    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    4.4127    1.6198    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.8634    1.3027    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.4127    2.1517    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  1  5  1  6      
  1  4  1  1      
  6  1  1  0      
  7  6  1  0      
  8  9  2  0      
  8  5  1  0      
  2 10  1  1      
  2 11  1  6      
 15 16  1  0      
 16 17  1  0      
 17 12  1  0      
 12 18  1  0      
 12 13  1  1      
 15 14  1  1      
 13 14  1  1      
 16 19  1  0      
 13 20  1  0      
 20 21  1  0      
 17 22  1  0      
  7 15  1  0      
 26 27  1  0      
 27 28  1  0      
 28 23  1  0      
 23 24  1  1      
 26 25  1  1      
 24 25  1  1      
 24 29  1  0      
 29 30  1  0      
 28 31  1  0      
 18 26  1  0      
 23 32  1  0      
 27 33  1  0      
 37 38  1  0      
 38 39  1  0      
 39 34  1  0      
 34 35  1  1      
 37 36  1  1      
 35 36  1  1      
 35 40  1  0      
 40 41  1  0      
 39 42  1  0      
 34 43  1  0      
 32 37  1  0      
 38 44  1  0      
 44 45  1  0      
 45 46  1  0      
 45 47  2  0      
M  END

ENDTEMPLATESTRING

  my($ArthroCerSn1TemplateString)=<<ENDTEMPLATESTRING;
ArthroCer sn1 acyl template structure
  LipdMAPS02070609152D

 45 47  0  0  0  0  0  0  0  0999 V2000
   -7.4742   -0.8124    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -8.1877   -0.4017    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -8.9013   -0.8124    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.0619   -1.5258    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -7.8867   -1.5258    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -6.7607   -0.4005    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.0471   -0.8124    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -7.7993    0.2709    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -8.5853    0.2867    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8433   -0.0925    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5101   -1.2476    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.8006   -0.9095    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.0832   -1.2764    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.4162   -0.1211    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1258   -0.4589    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1647   -0.2742    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.1721   -0.3234    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8486   -1.0704    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4317   -1.3111    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5432    0.1236    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.1415    0.4430    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.4747   -0.7119    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.1840   -0.3739    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0982   -0.7407    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5685    0.4144    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.8589    0.0766    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.1362   -0.5348    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.5531   -0.7755    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.4416    0.6591    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.9573    0.2244    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.5685    1.0614    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    8.0855    0.8962    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.4187   -0.2588    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.1280    0.0793    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.8458   -0.2876    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.5125    0.8676    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.8029    0.5299    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.0803   -0.0816    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.4972   -0.3223    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.3856    1.1124    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    8.9013    0.6777    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.7555    0.6648    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    4.1843    0.9945    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.6356    0.6778    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.1843    1.5258    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  1  5  1  6      
  1  4  1  1      
  6  1  1  0      
  7  6  1  0      
  2  8  1  1      
  2  9  1  6      
 13 14  1  0      
 14 15  1  0      
 15 10  1  0      
 10 16  1  0      
 10 11  1  1      
 13 12  1  1      
 11 12  1  1      
 14 17  1  0      
 11 18  1  0      
 18 19  1  0      
 15 20  1  0      
  7 13  1  0      
 24 25  1  0      
 25 26  1  0      
 26 21  1  0      
 21 22  1  1      
 24 23  1  1      
 22 23  1  1      
 22 27  1  0      
 27 28  1  0      
 26 29  1  0      
 16 24  1  0      
 21 30  1  0      
 25 31  1  0      
 35 36  1  0      
 36 37  1  0      
 37 32  1  0      
 32 33  1  1      
 35 34  1  1      
 33 34  1  1      
 33 38  1  0      
 38 39  1  0      
 37 40  1  0      
 32 41  1  0      
 30 35  1  0      
 36 42  1  0      
 42 43  1  0      
 43 44  1  0      
 43 45  2  0      
M  END

ENDTEMPLATESTRING

  # MIPCer templates...
  my($MIPCerSn1Sn2TemplateString)=<<ENDTEMPLATESTRING;
MIPCer sn1 acyl and sn2 acyl template structure
  LipdMAPS08061209252D

 37 38  0  0  0  0  0  0  0  0999 V2000
   -3.8973   -0.8873    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.6116   -0.4761    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.3260   -0.8873    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.4844   -1.6015    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -4.3102   -1.6015    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1828   -0.4749    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4684   -0.8873    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6998   -0.9045    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4518   -0.5984    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8159   -1.2291    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4518    0.1539    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.0584   -2.0237    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.0584   -2.8497    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.2227    0.1974    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.0095    0.2132    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.0947   -0.5780    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7840   -0.2278    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.1033   -1.4123    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4155   -1.0590    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.7295   -1.4123    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.4070   -0.2278    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0126   -0.4412    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.4115   -0.0135    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.5756    0.0188    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.6427   -1.1130    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.0954   -0.4016    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.5954    2.2670    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.9405    1.1322    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.6727    1.4641    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.4127    1.1040    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0678    2.2389    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.3357    1.9071    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.2622    2.0885    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.9081    2.4795    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.7265    1.3426    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0678    2.8497    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.3260    0.9968    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
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
 12  5  1  0      
  9  7  1  0      
  2 14  1  1      
  2 15  1  6      
 16 17  1  0      
 16 21  1  0      
 17 18  1  0      
 20 21  1  0      
 20 19  1  1      
 18 19  1  1      
 17 22  1  0      
 19 23  1  0      
 16 24  1  0      
 20 25  1  0      
 21 26  1  0      
  8 18  1  0      
 30 31  1  0      
 31 32  1  0      
 32 27  1  0      
 27 33  1  0      
 27 28  1  1      
 30 29  1  1      
 28 29  1  1      
 32 34  1  0      
 28 35  1  0      
 31 36  1  0      
 35 37  1  0      
 30 23  1  0      
M  END

ENDTEMPLATESTRING

  my($MIPCerSn1TemplateString)=<<ENDTEMPLATESTRING;
MIPCer sn1 acyl template structure
  LipdMAPS08061209252D

 35 36  0  0  0  0  0  0  0  0999 V2000
   -3.8926   -1.5096    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.6060   -1.0989    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.3196   -1.5096    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.4802   -2.2229    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -4.3050   -2.2229    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1789   -1.0977    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4654   -1.5096    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6990   -1.5267    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4500   -1.2210    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8137   -1.8510    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4500   -0.4696    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.2176   -0.4262    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.0034   -0.4104    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.0922   -1.2006    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7831   -0.8509    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.1032   -2.0339    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4138   -1.6811    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.7262   -2.0339    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.4029   -0.8509    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0126   -1.0640    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.4098   -0.6368    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.5725   -0.6046    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.6383   -1.7350    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.0904   -1.0245    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.5898    1.6409    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.9357    0.5075    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.6695    0.8390    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.4110    0.4793    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0653    1.6128    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.3317    1.2814    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.2558    1.4626    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.9034    1.8532    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.7208    0.7176    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0653    2.2229    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.3196    0.3722    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
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
  2 12  1  1      
  2 13  1  6      
 14 15  1  0      
 14 19  1  0      
 15 16  1  0      
 18 19  1  0      
 18 17  1  1      
 16 17  1  1      
 15 20  1  0      
 17 21  1  0      
 14 22  1  0      
 18 23  1  0      
 19 24  1  0      
  8 16  1  0      
 28 29  1  0      
 29 30  1  0      
 30 25  1  0      
 25 31  1  0      
 25 26  1  1      
 28 27  1  1      
 26 27  1  1      
 30 32  1  0      
 26 33  1  0      
 29 34  1  0      
 33 35  1  0      
 28 21  1  0      
M  END

ENDTEMPLATESTRING

  # MIPCert templates...
  my($MIPCertSn1Sn2TemplateString)=<<ENDTEMPLATESTRING;
MIPCert sn1 acyl and sn2 acyl template structure
  LipdMAPS08061209252D

 38 39  0  0  0  0  0  0  0  0999 V2000
   -3.8926   -0.8862    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.6060   -0.4755    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.3196   -0.8862    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.4802   -1.5996    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -4.3050   -1.5996    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1789   -0.4743    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4654   -0.8862    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6990   -0.9034    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4500   -0.5977    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8137   -1.2276    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4500    0.1537    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.0523   -2.0212    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.0523   -2.8462    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.2176    0.1972    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.0034    0.2129    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.0922   -0.5773    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7831   -0.2275    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.1032   -1.4106    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4138   -1.0577    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.7262   -1.4106    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.4029   -0.2275    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0126   -0.4407    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.4098   -0.0135    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.5725    0.0188    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.6383   -1.1117    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.0904   -0.4011    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.5898    2.2643    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.9357    1.1308    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.6695    1.4623    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.4110    1.1027    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0653    2.2362    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.3317    1.9048    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.2558    2.0860    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.9034    2.4765    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.7208    1.3410    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0653    2.8462    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.3196    0.9956    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.3196   -1.4019    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
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
 12  5  1  0      
  9  7  1  0      
  2 14  1  1      
  2 15  1  6      
 16 17  1  0      
 16 21  1  0      
 17 18  1  0      
 20 21  1  0      
 20 19  1  1      
 18 19  1  1      
 17 22  1  0      
 19 23  1  0      
 16 24  1  0      
 20 25  1  0      
 21 26  1  0      
  8 18  1  0      
 30 31  1  0      
 31 32  1  0      
 32 27  1  0      
 27 33  1  0      
 27 28  1  1      
 30 29  1  1      
 28 29  1  1      
 32 34  1  0      
 28 35  1  0      
 31 36  1  0      
 35 37  1  0      
 30 23  1  0      
  3 38  1  6      
M  END

ENDTEMPLATESTRING

  my($MIPCertSn1TemplateString)=<<ENDTEMPLATESTRING;
MIPCert sn1 acyl  template structure
  LipdMAPS08061209252D

 36 37  0  0  0  0  0  0  0  0999 V2000
   -3.8972   -1.5113    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.6114   -1.1001    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.3259   -1.5113    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.4843   -2.2255    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -4.3101   -2.2255    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1826   -1.0989    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4683   -1.5113    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6998   -1.5285    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4517   -1.2224    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8158   -1.8531    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4517   -0.4702    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.2226   -0.4266    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.0093   -0.4109    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.0947   -1.2020    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7840   -0.8518    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.1033   -2.0363    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4155   -1.6830    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.7294   -2.0363    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.4069   -0.8518    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0126   -1.0653    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.4115   -0.6375    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.5755   -0.6052    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.6426   -1.7370    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.0952   -1.0256    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.5952    1.6429    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.9403    0.5081    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.6726    0.8400    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.4127    0.4800    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0677    1.6148    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.3356    1.2830    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.2620    1.4644    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.9080    1.8554    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.7263    0.7185    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0677    2.2255    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.3259    0.3727    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.3259   -2.2255    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
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
  2 12  1  1      
  2 13  1  6      
 14 15  1  0      
 14 19  1  0      
 15 16  1  0      
 18 19  1  0      
 18 17  1  1      
 16 17  1  1      
 15 20  1  0      
 17 21  1  0      
 14 22  1  0      
 18 23  1  0      
 19 24  1  0      
  8 16  1  0      
 28 29  1  0      
 29 30  1  0      
 30 25  1  0      
 25 31  1  0      
 25 26  1  1      
 28 27  1  1      
 26 27  1  1      
 30 32  1  0      
 26 33  1  0      
 29 34  1  0      
 33 35  1  0      
 28 21  1  0      
  3 36  1  6      
M  END

ENDTEMPLATESTRING

  # MIP2Cer templates...
  my($MIP2CerSn1Sn2TemplateString)=<<ENDTEMPLATESTRING;
MIP2Cer sn1 acyl and sn2 acyl template structure
  LipdMAPS08061209252D

 52 54  0  0  0  0  0  0  0  0999 V2000
   -7.2943   -0.8564    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -8.0086   -0.4451    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -8.7230   -0.8564    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.8815   -1.5706    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -7.7073   -1.5706    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -6.5798   -0.4439    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.8654   -0.8564    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.0969   -0.8736    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.8489   -0.5674    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -5.2129   -1.1981    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.8489    0.1849    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -8.4555   -1.9927    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -8.4555   -2.8188    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -7.6198    0.2284    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -8.4066    0.2442    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3023   -0.5470    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6131   -0.1968    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.2938   -1.3814    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9815   -1.0280    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6675   -1.3814    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0100   -0.1968    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.4096   -0.4102    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9855    0.0174    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8214    0.0497    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.2457   -1.0820    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.6984   -0.3706    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.2021    2.2361    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5472    1.1014    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7207    1.4333    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9805    1.0731    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3255    2.2081    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0576    1.8763    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.8689    2.0576    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.5147    2.4485    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.3331    1.3117    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3255    2.8188    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.9325    0.9659    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.9448    0.9948    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.1960    1.2995    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    2.9616    0.6080    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.1960    2.0481    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    6.7336    1.3117    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.4301    1.6601    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.7530    0.4822    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.0580    0.8336    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.3646    0.4822    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.0381    1.6601    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.6379    1.4483    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    6.0540    1.8732    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.2116    1.9053    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    8.2727    0.7799    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    8.7230    1.4874    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
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
 12  5  1  0      
  9  7  1  0      
  2 14  1  1      
  2 15  1  6      
 16 17  1  0      
 16 21  1  0      
 17 18  1  0      
 20 21  1  0      
 20 19  1  1      
 18 19  1  1      
 17 22  1  0      
 19 23  1  0      
 16 24  1  0      
 20 25  1  0      
 21 26  1  0      
  8 18  1  0      
 30 31  1  0      
 31 32  1  0      
 32 27  1  0      
 27 33  1  0      
 27 28  1  1      
 30 29  1  1      
 28 29  1  1      
 32 34  1  0      
 28 35  1  0      
 31 36  1  0      
 35 37  1  0      
 39 38  1  0      
 39 40  1  0      
 39 41  2  0      
 42 43  1  0      
 42 47  1  0      
 43 44  1  0      
 46 47  1  0      
 46 45  1  1      
 44 45  1  1      
 43 48  1  0      
 45 49  1  0      
 42 50  1  0      
 46 51  1  0      
 47 52  1  0      
 44 38  1  0      
 37 39  1  0      
 30 23  1  0      
M  END

ENDTEMPLATESTRING

  my($MIP2CerSn1TemplateString)=<<ENDTEMPLATESTRING;
MIP2Cer sn1 acyl template structure
  LipdMAPS08061209252D

 50 52  0  0  0  0  0  0  0  0999 V2000
   -6.3334   -1.2855    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.9536   -0.9283    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.5738   -1.2855    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.9749   -1.9056    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -6.6920   -1.9056    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -5.7130   -0.9273    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.0927   -1.2855    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.5572   -1.3004    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.2101   -1.0345    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -4.5262   -1.5821    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.2101   -0.3813    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -6.6160   -0.3436    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -7.2991   -0.3299    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1307   -1.0168    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2689   -0.7128    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.8599   -1.7413    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7205   -1.4345    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5796   -1.7413    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0087   -0.7128    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.9604   -0.8980    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7239   -0.5268    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7132   -0.4987    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.2133   -1.4813    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.6064   -0.8637    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.0437    1.3996    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4751    0.4144    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6258    0.7026    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7196    0.3898    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1509    1.3753    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0500    1.0872    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.6227    1.2447    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.4469    1.5841    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.1575    0.5970    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1509    1.9056    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.6779    0.2968    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.4251    0.3219    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.7750    0.5864    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    2.5714   -0.0140    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.7750    1.2364    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.8465    0.5970    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.7147    0.8995    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.1268   -0.1232    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.2599    0.1819    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.3944   -0.1232    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.9792    0.8995    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.0269    0.7156    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.2565    1.0845    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    6.2616    1.1124    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.1829    0.1353    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.5738    0.7496    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
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
  2 12  1  1      
  2 13  1  6      
 14 15  1  0      
 14 19  1  0      
 15 16  1  0      
 18 19  1  0      
 18 17  1  1      
 16 17  1  1      
 15 20  1  0      
 17 21  1  0      
 14 22  1  0      
 18 23  1  0      
 19 24  1  0      
  8 16  1  0      
 28 29  1  0      
 29 30  1  0      
 30 25  1  0      
 25 31  1  0      
 25 26  1  1      
 28 27  1  1      
 26 27  1  1      
 30 32  1  0      
 26 33  1  0      
 29 34  1  0      
 33 35  1  0      
 37 36  1  0      
 37 38  1  0      
 37 39  2  0      
 40 41  1  0      
 40 45  1  0      
 41 42  1  0      
 44 45  1  0      
 44 43  1  1      
 42 43  1  1      
 41 46  1  0      
 43 47  1  0      
 40 48  1  0      
 44 49  1  0      
 45 50  1  0      
 42 36  1  0      
 35 37  1  0      
 28 21  1  0      
M  END

ENDTEMPLATESTRING

  # MIP2Cert templates...
  my($MIP2CertSn1Sn2TemplateString)=<<ENDTEMPLATESTRING;
MIP2Cert sn1 acyl and sn2 acyl template structure
  LipdMAPS08061209252D

 53 55  0  0  0  0  0  0  0  0999 V2000
   -6.3334   -0.7436    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.9536   -0.3865    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.5738   -0.7436    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.9749   -1.3637    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -6.6920   -1.3637    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -5.7130   -0.3854    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.0927   -0.7436    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.5572   -0.7585    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.2101   -0.4927    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -4.5262   -1.0403    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.2101    0.1605    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -7.3416   -1.7302    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.3416   -2.4475    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -6.6160    0.1983    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -7.2991    0.2120    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1307   -0.4749    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2689   -0.1709    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.8599   -1.1994    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7205   -0.8926    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5796   -1.1994    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0087   -0.1709    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.9604   -0.3562    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7239    0.0151    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7132    0.0432    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.2133   -0.9395    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.6064   -0.3218    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.0437    1.9415    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4751    0.9563    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6258    1.2445    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7196    0.9317    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1509    1.9172    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0500    1.6291    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.6227    1.7865    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.4469    2.1259    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.1575    1.1389    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1509    2.4475    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.6779    0.8387    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.4251    0.8637    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.7750    1.1283    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    2.5714    0.5279    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.7750    1.7783    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.8465    1.1389    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.7147    1.4414    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.1268    0.4187    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.2599    0.7238    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.3944    0.4187    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.9792    1.4414    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.0269    1.2575    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.2565    1.6264    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    6.2616    1.6543    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.1829    0.6772    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.5738    1.2915    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -7.5738   -1.2386    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
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
 12  5  1  0      
  9  7  1  0      
  2 14  1  1      
  2 15  1  6      
 16 17  1  0      
 16 21  1  0      
 17 18  1  0      
 20 21  1  0      
 20 19  1  1      
 18 19  1  1      
 17 22  1  0      
 19 23  1  0      
 16 24  1  0      
 20 25  1  0      
 21 26  1  0      
  8 18  1  0      
 30 31  1  0      
 31 32  1  0      
 32 27  1  0      
 27 33  1  0      
 27 28  1  1      
 30 29  1  1      
 28 29  1  1      
 32 34  1  0      
 28 35  1  0      
 31 36  1  0      
 35 37  1  0      
 39 38  1  0      
 39 40  1  0      
 39 41  2  0      
 42 43  1  0      
 42 47  1  0      
 43 44  1  0      
 46 47  1  0      
 46 45  1  1      
 44 45  1  1      
 43 48  1  0      
 45 49  1  0      
 42 50  1  0      
 46 51  1  0      
 47 52  1  0      
 44 38  1  0      
 37 39  1  0      
 30 23  1  0      
  3 53  1  6      
M  END

ENDTEMPLATESTRING

  my($MIP2CertSn1TemplateString)=<<ENDTEMPLATESTRING;
MIP2Cert sn1 acyl template structure
  LipdMAPS08061209252D

 51 53  0  0  0  0  0  0  0  0999 V2000
   -7.0054   -1.4219    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.6914   -1.0269    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -8.3775   -1.4219    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.6089   -2.1078    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -7.4021   -2.1078    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -6.3192   -1.0257    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.6331   -1.4219    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.9347   -1.4384    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.6568   -1.1444    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
   -5.0065   -1.7501    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.6568   -0.4219    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -7.3180   -0.3801    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -8.0736   -0.3649    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2507   -1.1247    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5097   -0.7884    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1634   -1.9261    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9031   -1.5867    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6411   -1.9261    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0096   -0.7884    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.2745   -0.9934    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9068   -0.5827    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7889   -0.5516    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.2359   -1.6386    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.6707   -0.9553    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.1544    1.5481    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5255    0.4584    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6922    0.7772    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9021    0.4312    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2730    1.5212    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0553    1.2026    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.7949    1.3767    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.4943    1.7521    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.2803    0.6603    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2730    2.1078    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.8559    0.3283    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.7885    0.3559    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.0695    0.6486    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    2.8443   -0.0155    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.0695    1.3676    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    6.4669    0.6603    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.2150    0.9949    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.5647   -0.1363    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.8180    0.2012    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.0729   -0.1363    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.7198    0.9949    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.4542    0.7915    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.8143    1.1996    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    6.9260    1.2304    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.9451    0.1497    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    8.3775    0.8291    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -8.3775   -2.1078    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
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
  2 12  1  1      
  2 13  1  6      
 14 15  1  0      
 14 19  1  0      
 15 16  1  0      
 18 19  1  0      
 18 17  1  1      
 16 17  1  1      
 15 20  1  0      
 17 21  1  0      
 14 22  1  0      
 18 23  1  0      
 19 24  1  0      
  8 16  1  0      
 28 29  1  0      
 29 30  1  0      
 30 25  1  0      
 25 31  1  0      
 25 26  1  1      
 28 27  1  1      
 26 27  1  1      
 30 32  1  0      
 26 33  1  0      
 29 34  1  0      
 33 35  1  0      
 37 36  1  0      
 37 38  1  0      
 37 39  2  0      
 40 41  1  0      
 40 45  1  0      
 41 42  1  0      
 44 45  1  0      
 44 43  1  1      
 42 43  1  1      
 41 46  1  0      
 43 47  1  0      
 40 48  1  0      
 44 49  1  0      
 45 50  1  0      
 42 36  1  0      
 35 37  1  0      
 28 21  1  0      
  3 51  1  6      
M  END

ENDTEMPLATESTRING


# Format: ID => AbbrevID|HeadGroupName|HeadGroupNameBeforeBase|HeadGroupAbbrev|Y1Sn1|Y2Sn1|Y1Sn2|Y2Sn2|Y1Sn3|Y2Sn3|Sn1ChainAtomNum|Sn2ChainAtomNum|Sn3ChainAtomNum|Sn1ChainCarbons|Sn2ChainCarbons|Sn3ChainCarbons|Sn2CAtomNum|Sn2NAtomNum|Sn2HAtomNum|LMCategory|LMMainClass|LMSubClass|TemplateString
#
# Note:
#  . Based on the sphingoid base abbreviation at sn1, d18:1(4E) or d18:0, head group could
#    start with shing-4-enine or sphinganine. However, both options are not valid for some
#    subclasses: For Cert, Cerm, and PICert shing-4-enine is not allowed.
#
#  . LMMainClass and LMSubClass fields are intentionally left empty and based on the type of
#    sphingoid base, N-acyl chain, and head group, these are assigned dynamically.
#
#  . Templates Cert (4-hydroxy) and Cerm (3-keto) sphingoid bases.
#
  %SPTemplatesDataMap = (
			       "CerSn1Sn2" => "Cer||||0.4310|0.8422|-0.7055|-0.2928|0|0|3|8|0|4|1|0|1|5|4|SP|||$CerSn1Sn2TemplateString",
			       "CerSn1" => "Cer||||-0.1931|0.2181|0|0|0|0|3|0|0|4|0|0|1|5|4|SP|||$CerSn1TemplateString",

			       "CertSn1Sn2" => "Cert|||Cer|0.3747|0.7859|-0.8615|-0.4745|0|0|5|13|0|6|1|0|1|7|6|SP|||$CertSn1Sn2TemplateString",
			       "CertSn1" => "Cert|||Cer|-0.1931|0.2181|0|0|0|0|5|0|0|6|0|0|1|7|6|SP|||$CertSn1TemplateString",

			       "CermSn1Sn2" => "Cerm|||Cer|0.3963|0.8075|-0.7652|-0.3490|0|0|3|10|0|4|2|0|1|5|4|SP|||$CermSn1Sn2TemplateString",
			       "CermSn1" => "Cerm|||Cer|-0.1513|0.2599|0|0|0|0|3|0|0|4|0|0|1|5|4|SP|||$CermSn1TemplateString",

			       "CerPSn1Sn2" => "CerP|1-phosphate|||0.4310|0.8422|-0.7055|-0.2928|0|0|3|12|0|4|1|0|1|5|4|SP|||$CerPSn1Sn2TemplateString",
			       "CerPSn1" => "CerP|1-phosphate|||-0.1932|0.2182|0|0|0|0|3|0|0|4|0|0|1|5|4|SP|||$CerPSn1TemplateString",

			       "SMSn1Sn2" => "SM|1-phosphocholine|||0.4310|0.8422|-0.7055|-0.2928|0|0|3|18|0|4|1|0|1|5|4|SP|||$SMSn1Sn2TemplateString",
			       "SMSn1" => "SM|1-phosphocholine|||0.0771|0.4883|0|0|0|0|3|0|0|4|0|0|1|5|4|SP|||$SMSn1TemplateString",

			       "PECerSn1Sn2" => "PECer|1-phosphoethanolamine||PE-Cer|0.4310|0.8423|-0.7055|-0.2928|0|0|3|12|0|4|1|0|1|5|4|SP|||$PECerSn1Sn2TemplateString",
			       "PECerSn1" => "PECer|1-phosphoethanolamine||PE-Cer|-0.1932|0.2181|0|0|0|0|3|0|0|4|0|0|1|5|4|SP|||$PECerSn1TemplateString",

			       "PICerSn1Sn2" => "PICer|1-phospho-(1'-myo-inositol)||PI-Cer|0.4310|0.8422|-0.7055|-0.2928|0|0|3|12|0|4|1|0|1|5|4|SP|||$PICerSn1Sn2TemplateString",
			       "PICerSn1" => "PICer|1-phospho-(1'-myo-inositol)||PI-Cer|-0.1932|0.2182|0|0|0|0|3|0|0|4|0|0|1|5|4|SP|||$PICerSn1TemplateString",

			       "PICertSn1Sn2" => "PICert|1-phospho-(1'-myo-inositol)||PI-Cer|0.3747|0.7859|-0.8615|-0.4745|0|0|5|13|0|6|1|0|1|7|6|SP|||$PICertSn1Sn2TemplateString",
			       "PICertSn1" => "PICert|1-phospho-(1'-myo-inositol)||PI-Cer|-0.1931|0.2181|0|0|0|0|5|0|0|6|0|0|1|7|6|SP|||$PICertSn1TemplateString",

			       "GlcCerSn1Sn2" => "GlcCer|1-beta-glucosyl|1||0.4310|0.8422|-0.7055|-0.2928|0|0|3|8|0|4|1|0|1|5|4|SP|||$GlcCerSn1Sn2TemplateString",
			       "GlcCerSn1" => "GlcCer|1-beta-glucosyl|1||-0.1931|0.2181|0|0|0|0|3|0|0|4|0|0|1|5|4|SP|||$GlcCerSn1TemplateString",

			       "GalCerSn1Sn2" => "GalCer|1-beta-galactosyl|1||0.2693|0.6806|-0.8672|-0.4545|0|0|3|8|0|4|1|0|1|5|4|SP|||$GalCerSn1Sn2TemplateString",
			       "GalCerSn1" => "GalCer|1-beta-galactosyl|1||-0.3060|0.1051|0|0|0|0|3|0|0|4|0|0|1|5|4|SP|||$GalCerSn1TemplateString",

			       "LacCerSn1Sn2" => "LacCer|1-beta-lactosyl|1|Galb1-4Glcb-Cer|-0.0598|0.3515|-1.1963|-0.7836|0|0|3|8|0|4|1|0|1|5|4|SP|||$LacCerSn1Sn2TemplateString",
			       "LacCerSn1" => "LacCer|1-beta-lactosyl|1|Galb1-4Glcb-Cer|-0.6839|-0.2726|0|0|0|0|3|0|0|4|0|0|1|5|4|SP|||$LacCerSn1TemplateString",

			       "GM3CerSn1Sn2" => "GM3Cer|1-beta-GM3|1|NeuAca1-3Galb1-4Glcb-Cer|-1.2312|-0.8199|-2.3676|-1.9549|0|0|3|8|0|4|1|0|1|5|4|SP|||$GM3CerSn1Sn2TemplateString",
			       "GM3CerSn1" => "GM3Cer|1-beta-GM3|1|NeuAca1-3Galb1-4Glcb-Cer|-1.8554|-1.4441|0|0|0|0|3|0|0|4|0|0|1|5|4|SP|||$GM3CerSn1TemplateString",

			       "GM4CerSn1Sn2" => "GM4Cer|1-beta-GM4|1|NeuAca1-3Galb-Cer|-0.8703|-0.4590|-2.0068|-1.5940|0|0|3|8|0|4|1|0|1|5|4|SP|||$GM4CerSn1Sn2TemplateString",
			       "GM4CerSn1" => "GM4Cer|1-beta-GM4|1|NeuAca1-3Galb-Cer|-1.4946|-1.0832|0|0|0|0|3|0|0|4|0|0|1|5|4|SP|||$GM4CerSn1TemplateString",

			       "GB3CerSn1Sn2" => "GB3Cer|1-beta-GB3|1|Gala1-4Galb1-4Glcb-Cer|-1.5242|-1.1129|-2.6606|-2.2479|0|0|3|8|0|4|1|0|1|5|4|SP|||$GB3CerSn1Sn2TemplateString",
			       "GB3CerSn1" => "GB3Cer|1-beta-GB3|1|Gala1-4Galb1-4Glcb-Cer|-2.1457|-1.7349|0|0|0|0|3|0|0|4|0|0|1|5|4|SP|||$GB3CerSn1TemplateString",

			       "iGB3CerSn1Sn2" => "iGB3Cer|1-beta-iGB3|1|Gala1-3Galb1-4Glcb-Cer|-1.3489|-0.9376|-2.4853|-2.0726|0|0|3|8|0|4|1|0|1|5|4|SP|||$iGB3CerSn1Sn2TemplateString",
			       "iGB3CerSn1" => "iGB3Cer|1-beta-iGB3|1|Gala1-3Galb1-4Glcb-Cer|-1.9706|-1.5598|0|0|0|0|3|0|0|4|0|0|1|5|4|SP|||$iGB3CerSn1TemplateString",

			       "asialo-GM2CerSn1Sn2" => "asialo-GM2Cer|1-beta-asialo-GM2|1|GalNAcb1-4Galb1-4Glcb-Cer|-1.5224|-1.1116|-2.6574|-2.2452|0|0|3|8|0|4|1|0|1|5|4|SP|||$asialoGM2CerSn1Sn2TemplateString",
			       "asialo-GM2CerSn1" => "asialo-GM2Cer|1-beta-asialo-GM2|1|GalNAcb1-4Galb1-4Glcb-Cer|-2.1457|-1.7349|0|0|0|0|3|0|0|4|0|0|1|5|4|SP|||$asialoGM2CerSn1TemplateString",

			       "Lc3CerSn1Sn2" => "Lc3Cer|1-beta-Lc3|1|GlcNAcb1-4Galb1-4Glcb-Cer|-1.4645|-1.0537|-2.5996|-2.1874|0|0|3|8|0|4|1|0|1|5|4|SP|||$Lc3CerSn1Sn2TemplateString",
			       "Lc3CerSn1" => "Lc3Cer|1-beta-Lc3|1|GlcNAcb1-4Galb1-4Glcb-Cer|-2.0721|-1.6644|0|0|0|0|3|0|0|4|0|0|1|5|4|SP|||$Lc3CerSn1TemplateString",

			       "Manb1-4GlcCerSn1Sn2" => "Manb1-4GlcCer|1-beta-Manb4Glc|1|Manb1-4Glcb-Cer|0.0431|0.4544|-1.0933|-0.6806|0|0|3|8|0|4|1|0|1|5|4|SP|||$Manb14GlcCerSn1Sn2TemplateString",
			       "Manb1-4GlcCerSn1" => "Manb1-4GlcCer|1-beta-Manb4Glc|1|Manb1-4Glcb-Cer|-0.5810|-0.1697|0|0|0|0|3|0|0|4|0|0|1|5|4|SP|||$Manb14GlcCerSn1TemplateString",

			       "MolluCerSn1Sn2" => "MolluCer|1-beta-Mollu|1|Manb1-3Manb1-4Glcb-Cer|-1.1428|-0.7315|-2.2792|-1.8665|0|0|3|8|0|4|1|0|1|5|4|SP|||$MolluCerSn1Sn2TemplateString",
			       "MolluCerSn1" => "MolluCer|1-beta-Mollu|1|Manb1-3Manb1-4Glcb-Cer|-1.7241|-1.3227|0|0|0|0|3|0|0|4|0|0|1|5|4|SP|||$MolluCerSn1TemplateString",

			       "ArthroCerSn1Sn2" => "ArthroCer|1-beta-Arthro|1|GlcNAcb1-4Manb1-4Glc-Cer|-0.1893|0.2219|-1.3257| -0.9130|0|0|3|8|0|4|1|0|1|5|4|SP|||$ArthroCerSn1Sn2TemplateString",
			       "ArthroCerSn1" => "ArthroCer|1-beta-Arthro|1|GlcNAcb1-4Manb1-4Glc-Cer|-0.8124|-0.4017|0|0|0|0|3|0|0|4|0|0|1|5|4|SP|||$ArthroCerSn1TemplateString",

			       "MIPCerSn1Sn2" => "MIPCer|1-O-[D-mannopyranosyl-alpha1-2-myo-inositol-1-phosphate]||MIPC|-0.8873|-0.4761|-2.0237|-1.6110|0|0|3|12|0|4|1|0|1|5|4|SP|||$MIPCerSn1Sn2TemplateString",
			       "MIPCerSn1" => "MIPCer|1-O-[D-mannopyranosyl-alpha1-2-myo-inositol-1-phosphate]||MIPC|-1.5096|-1.0989|0|0|0|0|3|0|0|4|0|0|1|5|4|SP|||$MIPCerSn1TemplateString",

			       "MIPCertSn1Sn2" => "MIPCert|1-O-[D-mannopyranosyl-alpha1-2-myo-inositol-1-phosphate]||MIPC|-0.8862|-0.4755|-2.0212|-1.6085|0|0|3|12|0|4|1|0|1|5|4|SP|||$MIPCertSn1Sn2TemplateString",
			       "MIPCertSn1" => "MIPCert|1-O-[D-mannopyranosyl-alpha1-2-myo-inositol-1-phosphate]||MIPC|-1.5113|-1.1001|0|0|0|0|3|0|0|4|0|0|1|5|4|SP|||$MIPCertSn1TemplateString",

			       "MIP2CerSn1Sn2" => "MIP2Cer|1-O-[myo-inositol-1-phosphoryl-6-D-mannopyranosyl-alpha1-2-myo-inositol-1-phosphate]||M(IP)2C|-0.8564|-0.4451|-1.9927|-1.5800|0|0|3|12|0|4|1|0|1|5|4|SP|||$MIP2CerSn1Sn2TemplateString",
			       "MIP2CerSn1" => "MIP2Cer|1-O-[myo-inositol-1-phosphoryl-6-D-mannopyranosyl-alpha1-2-myo-inositol-1-phosphate]||M(IP)2C|-1.2855|-0.9283|0|0|0|0|3|0|0|4|0|0|1|5|4|SP|||$MIP2CerSn1TemplateString",

			       "MIP2CertSn1Sn2" => "MIP2Cert|1-O-[myo-inositol-1-phosphoryl-6-D-mannopyranosyl-alpha1-2-myo-inositol-1-phosphate]||M(IP)2C|-0.7436|-0.3865|-1.7302|-1.3275|0|0|3|12|0|4|1|0|1|5|4|SP|||$MIP2CertSn1Sn2TemplateString",
			       "MIP2CertSn1" => "MIP2Cert|1-O-[myo-inositol-1-phosphoryl-6-D-mannopyranosyl-alpha1-2-myo-inositol-1-phosphate]||M(IP)2C|-1.4219|-1.0269|0|0|0|0|3|0|0|4|0|0|1|5|4|SP|||$MIP2CertSn1TemplateString",
			       );
}
# Template format: ID => AbbrevID|HeadGroupName|HeadGroupNameBeforeBase|HeadGroupAbbrev|Y1Sn1|Y2Sn1|Y1Sn2|Y2Sn2|Y1Sn3|Y2Sn3|Sn1ChainAtomNum|Sn2ChainAtomNum|Sn3ChainAtomNum|Sn1ChainCarbons|Sn2ChainCarbons|Sn3ChainCarbons|Sn2CAtomNum|Sn2OAtomNum|Sn2HAtomNum|LMCategory|LMMainClass|LMSubClass|TemplateString

# Initialize supported head groups...
sub _InitializeSupportedHeadGroupsData {
  my($SPType, $SPHeadGroup);
  %SPSupportedHeadGroupMap = ();

  for $SPType (keys %SPTemplatesDataMap) {
    ($SPHeadGroup) = split /\|/, $SPTemplatesDataMap{$SPType};
    if (!(exists $SPSupportedHeadGroupMap{$SPHeadGroup})) {
      $SPSupportedHeadGroupMap{$SPHeadGroup} = $SPHeadGroup;
    }
  }
}


1;

__END__

=head1 NAME

SPStr - Sphingolipids (SP) structure generation methods

=head1 SYNOPSIS

use SPStr;

use SPStr qw(:all);

=head1 DESCRIPTION

SPStr module provides these methods:

    GenerateCmpdOntologyData - Generate ontology data
    GenerateCmpdOntologySDDataLines - Generate ontology data lines for SD file
    GenerateSPChainStrData - Generate chain structure data
    GenerateSDFile - Generate SD file
    GetSPTemplatesData - Get templates data
    GetSPSupportedHeadGroupMap - Get supported headgroups data
    GetSPTemplateID - Get templates ID
    IsSPChainsAbbrevSupported - Is it a supported SP abbreviation
    ParseSPAbbrev - Parse SP abbreviation
    ProcessSPCmpdAbbrevs - Process SP abbreviation
    SetupSPCmpdAbbrevTemplateDataMap - Setup template structure data map
    ValidateSPAbbrev - Validate SP abbreviation

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

=item B<GenerateSPChainStrData>

    ($AtomLinesArrayRef, $BondLinesArrayRef) =
       GenerateSPChainStrData($ChainType, $CmpdDataRef);

Return array references containing atom and bond data lines for SD file. Appropriate atom
and bond data lines are generated using chain type and abbreviation template data.

=item B<GenerateSDFile>

    GenerateSDFile($SDFileName, $CmdAbbrevsRef);

Generate a SD file for compound abbreviations. Structure data for specified abbreviation
is generated sequentially and written to SD file.

=item B<GetSPTemplatesData>

    $TemplatesDataRef = GetSPTemplatesData();

Return a reference to a hash containing SP templates data

=item B<GetSPSupportedHeadGroupMap>

    $SupportedHeadGroupDataRef = GetSPSupportedHeadGroupMap();

Return a reference to a hash containing supported head groups data.

=item B<GetSPTemplateID>

    $HeadGroupID = GetSPTemplateID($HeadGroupAbbrev, $ChainsAbbrev);

Return a supported template ID for compound abbreviation.

=item B<IsSPChainsAbbrevSupported>

    $Status = IsSPChainsAbbrevSupported($Abbrev, $PrintWarning);

Return 1 or 0 based on whether SP abbreviated is supported. For unsupported SP abbreviations,
a warning is printed unless PrintWarning flag is set.

=item B<ParseSPAbbrev>

    ($HeadGroup, $ChainsAbbrev, $AbbrevModifier) =
       ParseSPAbbrev($Abbrev);

Parse SP abbreviation and return these values: HeadGroup, ChainsAbbrev,
AbbrevModifier.

=item B<ProcessSPCmpdAbbrevs>

    ProcessSPCmpdAbbrevs($CmpdAbbrevsRef, $AllowArbitraryChainAbbrev,
                         $WriteSDFile, $SDFileName);

Process specified SP abbreviations to generate structures and write them out either
a SD file or simply report number of valid abbreviations.

=item B<SetupSPCmpdAbbrevTemplateDataMap>

    $AbbrevTemplateDataMapRef =
       SetupSPCmpdAbbrevTemplateDataMap($Abbrev);

Return a reference to a hash containing template data for compound abbreviation. The
template data is used to generate SD file for compound abbreviation.

=item B<ValidateSPAbbrev>

    $Status = ValidateSPAbbrev($Abbrev);

Return 1 or 0 based on whether a SP abbreviation is valid.

=back

=head1 AUTHOR

Manish Sud

=head1 CONTRIBUTOR

Eoin Fahy

=head1 SEE ALSO

ChainStr.pm

=head1 COPYRIGHT

Copyright (C) 2006-2012. The Regents of the University of California. All Rights Reserved.

=head1 LICENSE

Modified BSD License

=cut
