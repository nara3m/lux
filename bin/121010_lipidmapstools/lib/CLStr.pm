package CLStr;
#
# File: CLStr.pm
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
@EXPORT_OK = qw(GenerateCmpdOntologyData GenerateCmpdOntologySDDataLines GenerateCLChainStrData GenerateSDFile GetCLTemplatesData GetCLSupportedHeadGroupMap GetCLTemplateID IsCLChainsAbbrevSupported ParseCLAbbrev ProcessCLCmpdAbbrevs SetupCLChainsAbbrev SetupCLCmpdAbbrevTemplateDataMap ValidateCLAbbrev);
%EXPORT_TAGS = (
		all  => [@EXPORT, @EXPORT_OK],
		);

# Initialize Cl data...
my(%CLTemplatesDataMap, %CLSupportedHeadGroupMap);
_InitializeData();

# Process specified compound abbreviations containing any wild cards
# and count the number of valid abbreviations along with optional generation
# of a SD file.
#
sub ProcessCLCmpdAbbrevs {
  my($CmpdAbbrevsRef, $AllowArbitraryChainAbbrev, $WriteSDFile, $SDFileName) = @_;

  $AllowArbitraryChainAbbrev = defined $AllowArbitraryChainAbbrev ? $AllowArbitraryChainAbbrev : 0;

  $WriteSDFile = defined $WriteSDFile ? $WriteSDFile : 0;
  if ($WriteSDFile &&  IsEmpty($SDFileName)) {
    warn "Warning: CLStr::ProcessCLCmpdAbbrevs: No SD file name specified. Ingoring structure generation.\n";
    return;
  }

  if ($WriteSDFile) {
    print "Generating SD file $SDFileName...\n";

    open SDFILE, ">$SDFileName" or die "Error: Couldn't open $SDFileName: $! \n";

    _ProcessCLCmpdAbbrevs($CmpdAbbrevsRef, $AllowArbitraryChainAbbrev, $WriteSDFile, \*SDFILE);

    close SDFILE;
  }
  else {
    _ProcessCLCmpdAbbrevs($CmpdAbbrevsRef, $AllowArbitraryChainAbbrev, $WriteSDFile);
  }
}

# Process specified compound abbreviations containing any wild cards
# and count the number of valid abbreviations along with optional generation
# of a SD file.
#
sub _ProcessCLCmpdAbbrevs {
  my($CmpdAbbrevsRef, $AllowArbitraryChainAbbrev, $WriteSDFile, $SDFileRef) = @_;
  my($Abbrev, $AbbrevCount, $HeadGroupAbbrev, $CLChainsAbbrev, $ChainsAbbrev, $AbbrevModifier, $NewAbbrev, $Sn1Abbrev, $Sn2Abbrev, $Sn3Abbrev, $Sn4Abbrev, $AllowSubstituents, $AllowRings, @AbbrevList, @ChainsAbbrevList);

  $AbbrevCount = 0;

 ABBREV: for $Abbrev (@{$CmpdAbbrevsRef}) {
    ($Sn1Abbrev, $Sn2Abbrev, $Sn3Abbrev, $Sn4Abbrev) = ('0:0') x 4;

    if (!($Abbrev =~ /\(/ && $Abbrev =~ /\)/)) {
      warn "Warning: Ignored CL compound abbreviation $Abbrev due to incorrect format. Valid format: HeadGroupAbbrev(1'-[Chain1Abbrev/Chain2Abbrev],3'-[Chain3Abbrev/Chain4Abbrev])\n";
      next ABBREV;
    }
    if (!ValidateCLAbbrev($Abbrev)) {
      next ABBREV;
    }

    ($HeadGroupAbbrev, $ChainsAbbrev, $AbbrevModifier) = ParseCLAbbrev($Abbrev);

    # Headgroup is always CL...
    if ($HeadGroupAbbrev =~ /\*/) {
      $HeadGroupAbbrev = "CL";
    }

    @ChainsAbbrevList = ();
    @ChainsAbbrevList = split /\//, $ChainsAbbrev;
    if (@ChainsAbbrevList != 4) {
      warn "Warning: Ignored CL compound Abbreviation $Abbrev due to incorrect format: Must contain chain abbreviations for sn1 and sn2 chain positions at 1' and 3' glycerol positions...\n";
      next ABBREV;
    }
    ($Sn1Abbrev, $Sn2Abbrev, $Sn3Abbrev, $Sn4Abbrev) = @ChainsAbbrevList;
    $CLChainsAbbrev = SetupCLChainsAbbrev($Sn1Abbrev, $Sn2Abbrev, $Sn3Abbrev, $Sn4Abbrev);

    ($AllowSubstituents, $AllowRings) = (1, 0);
    if (!(ChainAbbrev::IsChainAbbrevOkay($Sn1Abbrev, $AllowSubstituents, $AllowRings, $AllowArbitraryChainAbbrev) && ChainAbbrev::IsChainAbbrevOkay($Sn2Abbrev, $AllowSubstituents, $AllowRings, $AllowArbitraryChainAbbrev) && ChainAbbrev::IsChainAbbrevOkay($Sn3Abbrev, $AllowSubstituents, $AllowRings, $AllowArbitraryChainAbbrev) && ChainAbbrev::IsChainAbbrevOkay($Sn4Abbrev, $AllowSubstituents, $AllowRings, $AllowArbitraryChainAbbrev))) {
      warn "Warning: Ignoring CL compound abbreviation $HeadGroupAbbrev(${CLChainsAbbrev})\n";
      next ABBREV;
    }
    if (!IsCLChainsAbbrevSupported($ChainsAbbrev, 1)) {
      warn "Warning: Ignoring CL compound abbreviation $HeadGroupAbbrev(${CLChainsAbbrev})\n";
      next ABBREV;
    }

    if (!(ChainAbbrev::IsWildCardInChainAbbrev($Sn1Abbrev) || ChainAbbrev::IsWildCardInChainAbbrev($Sn2Abbrev) || ChainAbbrev::IsWildCardInChainAbbrev($Sn3Abbrev) || ChainAbbrev::IsWildCardInChainAbbrev($Sn4Abbrev))) {
      my($TemplateHeadGroup) = GetCLTemplateID($HeadGroupAbbrev, $ChainsAbbrev);
      if (exists $CLTemplatesDataMap{$TemplateHeadGroup}) {
	$NewAbbrev = "$HeadGroupAbbrev($CLChainsAbbrev)$AbbrevModifier";
	$AbbrevCount++;
	if ($WriteSDFile) {
	  _GenerateAndWriteCmdDataString($SDFileRef, $NewAbbrev);
	}
      }
      else {
	warn "Warning: Ignored CL compound abbreviation $Abbrev : Abbreviation doesn't match any template\n";
      }
      next ABBREV;
    }

    # Arbitrary acyl chain abbreviation is not supported with wild cards...
    if ($AllowArbitraryChainAbbrev) {
      warn "Warning: Ignoring CL compound abbreviation $Abbrev : Allow arbitrary chain abbreviation option is not supported with wild cards in any part of the abbreviation\n";
      next ABBREV;
    }

    # Substituents are not supported with wild cards...
    if (ChainAbbrev::IsSubstituentInChainAbbrev($Sn1Abbrev) || ChainAbbrev::IsSubstituentInChainAbbrev($Sn2Abbrev) || ChainAbbrev::IsSubstituentInChainAbbrev($Sn3Abbrev) || ChainAbbrev::IsSubstituentInChainAbbrev($Sn4Abbrev)) {
      warn "Warning: Ignoring CL compound abbreviation $Abbrev : Substituent specifications are not supported with wild cards in any part of the abbreviation\n";
      next ABBREV;
    }

    # Get expanded abbreviation for each position...
    my($Sn1ExpandedAbbrevRef) = ChainAbbrev::ExpandChainAbbrev($Sn1Abbrev);
    my($Sn2ExpandedAbbrevRef) = ChainAbbrev::ExpandChainAbbrev($Sn2Abbrev);
    my($Sn3ExpandedAbbrevRef) = ChainAbbrev::ExpandChainAbbrev($Sn3Abbrev);
    my($Sn4ExpandedAbbrevRef) = ChainAbbrev::ExpandChainAbbrev($Sn4Abbrev);

    if ($Sn2Abbrev =~ /^(\*|\*:\*)$/i) {
      # Add 0:0 to Sn2Abbrev list containing wild cards...
      unshift(@{$Sn2ExpandedAbbrevRef}, "0:0")
    }

    if ($Sn4Abbrev =~ /^(\*|\*:\*)$/i) {
      # Add 0:0 to Sn4Abbrev list containing wild cards...
      unshift(@{$Sn4ExpandedAbbrevRef}, "0:0")
    }

    # Enumerate various possibilities...
    my($ExpandedAbbrev, $ExpandedSn1Abbrev, $ExpandedSn2Abbrev, $ExpandedSn3Abbrev, $ExpandedSn4Abbrev);

    SN1ABBREV: for $ExpandedSn1Abbrev (@$Sn1ExpandedAbbrevRef) {
      SN2ABBREV: for $ExpandedSn2Abbrev (@$Sn2ExpandedAbbrevRef) {
	if (ChainAbbrev::IsAlkylChainAbbrev($ExpandedSn2Abbrev) || ChainAbbrev::IsAlkenylChainAbbrev($ExpandedSn2Abbrev)) {
	  next SN2ABBREV;
	}
	if ($ExpandedSn2Abbrev =~ /^0:0$/i) {
	  if (!(ChainAbbrev::IsAlkylChainAbbrev($ExpandedSn1Abbrev) || ChainAbbrev::IsAlkenylChainAbbrev($ExpandedSn1Abbrev))) {
	    # Sn2 0:0 is only allowed for Sn1 alkyl...
	    next SN2ABBREV;
	  }
	}
	SN3ABBREV: for $ExpandedSn3Abbrev (@$Sn3ExpandedAbbrevRef) {
	  if (ChainAbbrev::IsAlkylChainAbbrev($ExpandedSn3Abbrev) || ChainAbbrev::IsAlkenylChainAbbrev($ExpandedSn3Abbrev)) {
	    next SN3ABBREV;
	  }
	  SN4ABBREV: for $ExpandedSn4Abbrev (@$Sn4ExpandedAbbrevRef) {
	    if (ChainAbbrev::IsAlkylChainAbbrev($ExpandedSn4Abbrev) || ChainAbbrev::IsAlkenylChainAbbrev($ExpandedSn4Abbrev)) {
	      next SN4ABBREV;
	    }
	    $CLChainsAbbrev = SetupCLChainsAbbrev($ExpandedSn1Abbrev, $ExpandedSn2Abbrev, $ExpandedSn3Abbrev, $ExpandedSn4Abbrev);
	    $ExpandedAbbrev = $HeadGroupAbbrev . '(' . $CLChainsAbbrev. ')' . $AbbrevModifier;
	    $AbbrevCount++;
	    if ($WriteSDFile) {
	      _GenerateAndWriteCmdDataString($SDFileRef, $ExpandedAbbrev);
	    }
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
  $CmpdAbbrevTemplateDataMapRef = SetupCLCmpdAbbrevTemplateDataMap($Abbrev);

  # Generate structure data for chains...
  my($Sn1AtomBlockLinesRef, $Sn1BondBlockLinesRef) = CLStr::GenerateCLChainStrData('Sn1', $CmpdAbbrevTemplateDataMapRef);
  my($Sn2AtomBlockLinesRef, $Sn2BondBlockLinesRef) = CLStr::GenerateCLChainStrData('Sn2', $CmpdAbbrevTemplateDataMapRef);
  my($Sn3AtomBlockLinesRef, $Sn3BondBlockLinesRef) = CLStr::GenerateCLChainStrData('Sn3', $CmpdAbbrevTemplateDataMapRef);
  my($Sn4AtomBlockLinesRef, $Sn4BondBlockLinesRef) = CLStr::GenerateCLChainStrData('Sn4', $CmpdAbbrevTemplateDataMapRef);

  # Generate data block lines including various desriptors...
  my($OntologyDataLinesRef) = GenerateCmpdOntologySDDataLines($CmpdAbbrevTemplateDataMapRef);

  # Setup first 4 SD file lines string: Name, MiscInfo, Comments, Count
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

# Generate atom and bond block lines for a chain
sub GenerateCLChainStrData {
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
  my($Sn1Abbrev, $Sn2Abbrev, $Sn3Abbrev, $Sn4Abbrev, $Sn1Name, $Sn2Name, $Sn3Name, $Sn4Name, $HeadGroupName, $ChainName, $ChainName1, $ChainName2, $NamePart1, $NamePart2, $NameSuffix, $HeadGroupAbbrev, $ChainAbbrevToNameMap);

  $Sn1Abbrev = $CmpdAbbrevTemplateDataMapRef->{SnChainAdd}[0] ? $CmpdAbbrevTemplateDataMapRef->{SnAbbrev}[0] : '0:0';
  $Sn2Abbrev = $CmpdAbbrevTemplateDataMapRef->{SnChainAdd}[1] ? $CmpdAbbrevTemplateDataMapRef->{SnAbbrev}[1] : '0:0';
  $Sn3Abbrev = $CmpdAbbrevTemplateDataMapRef->{SnChainAdd}[2] ? $CmpdAbbrevTemplateDataMapRef->{SnAbbrev}[2] : '0:0';
  $Sn4Abbrev = $CmpdAbbrevTemplateDataMapRef->{SnChainAdd}[3] ? $CmpdAbbrevTemplateDataMapRef->{SnAbbrev}[3] : '0:0';

  $Sn1Name = ChainAbbrev::SetupChainNameWithSubstituents($CmpdAbbrevTemplateDataMapRef, 0);
  $Sn2Name = ChainAbbrev::SetupChainNameWithSubstituents($CmpdAbbrevTemplateDataMapRef, 1);
  $Sn3Name = ChainAbbrev::SetupChainNameWithSubstituents($CmpdAbbrevTemplateDataMapRef, 2);
  $Sn4Name = ChainAbbrev::SetupChainNameWithSubstituents($CmpdAbbrevTemplateDataMapRef, 3);

  $HeadGroupAbbrev = $CmpdAbbrevTemplateDataMapRef->{AbbrevID};
  $HeadGroupName = $CmpdAbbrevTemplateDataMapRef->{HeadGroupName};

  my($AbbrevModifier) = $CmpdAbbrevTemplateDataMapRef->{AbbrevModifier};
  if ($AbbrevModifier =~ /^rac$/i) {
    $NameSuffix = '-rac-glycerol';
  }
  elsif ($AbbrevModifier =~ /^U$/i) {
    $NameSuffix = '-glycerol';
  }
  else {
    $NameSuffix = '-sn-glycerol';
  }

  my($SystematicName) = '';
  if ($Sn1Abbrev eq $Sn2Abbrev && $Sn2Abbrev eq $Sn3Abbrev && $Sn3Abbrev eq $Sn4Abbrev) {
    $ChainName = ($Sn1Name =~ /^\([1-9]/) ? "di-${Sn1Name}" : "di${Sn1Name}";
    $SystematicName = "1',3'-Bis-[1,2-${ChainName}-sn-glycero-3-phospho]${NameSuffix}"
  }
  elsif ($Sn1Abbrev eq $Sn3Abbrev && $Sn2Abbrev eq $Sn4Abbrev) {
    # Sn1 and Sn3, Sn2 and Sn4 are same...
    $ChainName1 = $Sn1Name ? "1-${Sn1Name}-" : '';
    $ChainName2 = $Sn2Name ? "2-${Sn2Name}-" : ($Sn1Name ? '' : '-');
    $SystematicName = "1',3'-Bis-[${ChainName1}${ChainName2}sn-glycero-3-phospho]${NameSuffix}"
  }
  elsif ($Sn1Abbrev eq $Sn2Abbrev && $Sn3Abbrev eq $Sn4Abbrev) {
    # Sn1 and Sn2, Sn3 and Sn4 are same...
    $ChainName1 = ($Sn1Name =~ /^\([1-9]/) ? "di-${Sn1Name}" : "di${Sn1Name}";
    $ChainName2 = ($Sn2Name =~ /^\([1-9]/) ? "di-${Sn2Name}" : "di${Sn2Name}";

    $NamePart1 = $Sn1Name ? "1'-[1,2-${ChainName1}-sn-glycero-3-phospho]," : '';
    $NamePart2 = $Sn1Name ? "3'-[1,2-${ChainName2}-sn-glycero-3-phospho]" : '';
    $SystematicName = "${NamePart1}${NamePart2}${NameSuffix}"
  }
  elsif ($Sn1Abbrev eq $Sn2Abbrev) {
    # Sn1 and Sn2 are same...
    $ChainName = ($Sn1Name =~ /^\([1-9]/) ? "di-${Sn1Name}" : "di${Sn1Name}";
    $NamePart1 = $Sn1Name ? "1'-[1,2-${ChainName}-sn-glycero-3-phospho]," : '';

    $ChainName1 = $Sn3Name ? "1-${Sn3Name}-" : '';
    $ChainName2 = $Sn4Name ? "2-${Sn4Name}-" : ($Sn3Name ? '' : '-');
    $NamePart2 = $Sn1Name ? "3'-[${ChainName1}${ChainName2}sn-glycero-3-phospho]" : '';

    $SystematicName = "${NamePart1}${NamePart2}${NameSuffix}"
  }
  elsif ($Sn3Abbrev eq $Sn4Abbrev) {
    # Sn3 and Sn4 are same...
    $ChainName1 = $Sn1Name ? "1-${Sn1Name}-" : '';
    $ChainName2 = $Sn2Name ? "2-${Sn2Name}-" : ($Sn1Name ? '' : '-');
    $NamePart1 = $Sn1Name ? "1'-[${ChainName1}${ChainName2}sn-glycero-3-phospho]," : '';

    $ChainName = ($Sn3Name =~ /^\([1-9]/) ? "di-${Sn3Name}" : "di${Sn3Name}";
    $NamePart2 = $Sn1Name ? "3'-[1,2-${ChainName}-sn-glycero-3-phospho]" : '';

    $SystematicName = "${NamePart1}${NamePart2}${NameSuffix}"
  }
  else {
    # All four chains are different...
    $ChainName1 = $Sn1Name ? "1-${Sn1Name}-" : '';
    $ChainName2 = $Sn2Name ? "2-${Sn2Name}-" : ($Sn1Name ? '' : '-');
    $NamePart1 = $Sn1Name ? "1'-[${ChainName1}${ChainName2}sn-glycero-3-phospho]," : '';

    $ChainName1 = $Sn3Name ? "1-${Sn3Name}-" : '';
    $ChainName2 = $Sn4Name ? "2-${Sn4Name}-" : ($Sn3Name ? '' : '-');
    $NamePart2 = $Sn1Name ? "3'-[${ChainName1}${ChainName2}sn-glycero-3-phospho]" : '';
    $SystematicName = "${NamePart1}${NamePart2}${NameSuffix}"
  }

  if (IsEmpty($AbbrevModifier)) {
    $OntologyDataMapRef->{Abbrev} = "$HeadGroupAbbrev(1'-[$Sn1Abbrev/$Sn2Abbrev],3'-[$Sn3Abbrev/$Sn4Abbrev])";
  }
  else {
    $OntologyDataMapRef->{Abbrev} = "$HeadGroupAbbrev(1'-[$Sn1Abbrev/$Sn2Abbrev],3'-[$Sn3Abbrev/$Sn4Abbrev])[${AbbrevModifier}]";
  }
  $OntologyDataMapRef->{'Systematic Name'} = "$SystematicName";
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

  @ChainIndices = (0, 1, 2, 3);

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

# Return reference to CLTemplatesDataMap...
sub GetCLTemplatesData {
  return \%CLTemplatesDataMap;
}

# Return reference to CLSupportedHeadGroupMap...
sub GetCLSupportedHeadGroupMap {
  return \%CLSupportedHeadGroupMap;
}

# Get CL template ID...
#
# Various possibilites:
#
#    . Diacylglycerophosphoglycerophosphodiradylglycerols
#    . Diacylglycerophosphoglycerophosphomonoradylglycerols
#    . 1-alkyl,2-acylglycerophosphoglycerophosphodiradylglycerols
#    . 1-alkyl,2-acylglycerophosphoglycerophosphomonoradylglycerols
#    . 1Z-alkenyl,2-acylglycerophosphoglycerophosphodiradylglycerols
#    . 1Z-alkenyl,2-acylglycerophosphoglycerophosphomonoradylglycerols
#    . Monoacylglycerophosphoglycerophosphomonoradylglycerols
#    . 1-alkyl glycerophosphoglycerophosphodiradylglycerols
#    . 1-alkyl glycerophosphoglycerophosphomonoradylglycerols
#    . 1Z-alkenylglycerophosphoglycerophosphodiradylglycerols
#    . 1Z-alkenylglycerophosphoglycerophosphomonoradylglycerols
#
sub GetCLTemplateID {
  my($HeadGroupAbbrev, $ChainsAbbrev) = @_;
  my($HeadGroupID);

  $HeadGroupID = "";
  if ($HeadGroupAbbrev) {
    my(@AbbrevWords) = split /\//, $ChainsAbbrev;
    my ($Sn1Abbrev, $Sn2Abbrev, $Sn3Abbrev, $Sn4Abbrev) = @AbbrevWords;
    my($Sn3ID, $Sn4ID);

    $Sn3ID = 'Sn3'; $Sn4ID = 'Sn4';
    if ($Sn4Abbrev eq "0:0") {
      # For mono radyl positions...
      $Sn4ID = '';
    }

    if ($Sn2Abbrev eq "0:0") {
      $HeadGroupID = (ChainAbbrev::IsAlkylChainAbbrev($Sn1Abbrev)) ? "Sn1Alkyl" : ((ChainAbbrev::IsAlkenylChainAbbrev($Sn1Abbrev)) ? "Sn1Alkenyl" : "Sn1" ) ;
      $HeadGroupID = $HeadGroupAbbrev . $HeadGroupID . $Sn3ID . $Sn4ID;
    }
    else {
      $HeadGroupID = (ChainAbbrev::IsAlkylChainAbbrev($Sn1Abbrev) && ChainAbbrev::IsAlkylChainAbbrev($Sn2Abbrev)) ? "Sn1AlkylSn2Alkyl" : (ChainAbbrev::IsAlkylChainAbbrev($Sn1Abbrev)) ? "Sn1AlkylSn2" : ((ChainAbbrev::IsAlkenylChainAbbrev($Sn1Abbrev)) ? "Sn1AlkenylSn2" : "Sn1Sn2" ) ;
      $HeadGroupID = $HeadGroupAbbrev . $HeadGroupID . $Sn3ID . $Sn4ID;
    }
  }
  return $HeadGroupID;
}

# Does template exist to handle this abbreviation?
#
# Based on proposed LM subclasses, here is what's not allowed for CL.
#
# ForCL: alkenyl at sn2; alkyl at sn2 without alkyl at sn1; no acyl at sn3.
#
#
sub IsCLChainsAbbrevSupported {
  my($Abbrev, $PrintWarning) = @_;
  my($Sn1Abbrev, $Sn2Abbrev, $Sn3Abbrev, $Sn4Abbrev, @AbbrevList);

  ($Sn1Abbrev, $Sn2Abbrev, $Sn3Abbrev, $Sn4Abbrev) = ('0:0') x 4;
  @AbbrevList = split /\//, $Abbrev;
  if (@AbbrevList != 4) {
    if ($PrintWarning) {
      warn "Warning: Ignoring CL compound abbreviation $Abbrev : Must contain chain abbreviation for sn1, sn2, sn3 and sn4...\n";
    }
    return 0;
  }
  ($Sn1Abbrev, $Sn2Abbrev, $Sn3Abbrev, $Sn4Abbrev) = @AbbrevList;
  if (ChainAbbrev::IsAlkylChainAbbrev($Sn2Abbrev)) {
    if ($PrintWarning) {
      warn "Warning: Ignoring CL compound abbreviation $Abbrev : alkyl chain at sn2 position.\n";
    }
    return 0;
  }
  if (ChainAbbrev::IsAlkenylChainAbbrev($Sn2Abbrev)) {
    if ($PrintWarning) {
      warn "Warning: Ignoring CL compound abbreviation $Abbrev : alkenyl chain at sn2 position.\n";
    }
    return 0;
  }
  if (ChainAbbrev::IsAlkylChainAbbrev($Sn3Abbrev)) {
    if ($PrintWarning) {
      warn "Warning: Ignoring CL compound abbreviation $Abbrev : alkyl chain at sn3 position.\n";
    }
    return 0;
  }
  if (ChainAbbrev::IsAlkenylChainAbbrev($Sn3Abbrev)) {
    if ($PrintWarning) {
      warn "Warning: Ignoring CL compound abbreviation $Abbrev : alkenyl chain at sn3 position.\n";
    }
    return 0;
  }
  if (ChainAbbrev::IsAlkylChainAbbrev($Sn4Abbrev)) {
    if ($PrintWarning) {
      warn "Warning: Ignoring CL compound abbreviation $Abbrev : alkyl chain at sn4 position.\n";
    }
    return 0;
  }
  if (ChainAbbrev::IsAlkenylChainAbbrev($Sn4Abbrev)) {
    if ($PrintWarning) {
      warn "Warning: Ignoring CL compound abbreviation $Abbrev : alkenyl chain at sn4 position.\n";
    }
    return 0;
  }

  return 1;
}

# Parse abbreviation...
#
# Format: HeadGroupAbbrev(1'-[Chain1Abbrev/Chain2Abbrev],3'-[Chain3Abbrev/Chain4Abbrev])[AbbrevModifier]
#
sub ParseCLAbbrev {
  my($Abbrev) = @_;

  my($HeadGroup, $ChainsAbbrev, $Sn1ChainAbbev, $Sn2ChainAbbev, $Sn3ChainAbbev, $Sn4ChainAbbev, $AbbrevModifier);

  if ($Abbrev =~ /\]$/) {
    ($HeadGroup, $Sn1ChainAbbev, $Sn2ChainAbbev, $Sn3ChainAbbev, $Sn4ChainAbbev, $AbbrevModifier) = $Abbrev =~ /^(.*?)\(1\'-\[(.*?)\/(.*?)\]\,3\'-\[(.*?)\/(.*?)\]\)\[(.*?)$/;
    $AbbrevModifier = '[' . $AbbrevModifier;
  }
  else {
    ($HeadGroup, $Sn1ChainAbbev, $Sn2ChainAbbev, $Sn3ChainAbbev, $Sn4ChainAbbev) = $Abbrev =~ /^(.*?)\(1\'-\[(.*?)\/(.*?)\]\,3\'-\[(.*?)\/(.*?)\]\)$/;
    $AbbrevModifier = '';
  }
  $ChainsAbbrev = "$Sn1ChainAbbev/$Sn2ChainAbbev/$Sn3ChainAbbev/$Sn4ChainAbbev";

  return ($HeadGroup, $ChainsAbbrev, $AbbrevModifier);
}
#
# Check out the validity of CL abbreviation...
sub ValidateCLAbbrev {
  my($Abbrev) = @_;

  if (!($Abbrev =~ /\(/ && $Abbrev =~ /\)/)) {
    warn "Warning: Ignored CL compound abbreviation $Abbrev due to incorrect format. Valid format: HeadGroupAbbrev(1'-[Chain1Abbrev/Chain2Abbrev],3'-[Chain3Abbrev/Chain4Abbrev])\n";
    return 0;
  }

  # Make sure head group is these...
  if ($Abbrev =~ /^\(/) {
    warn "Warning: Ignored CL compound abbreviation $Abbrev due to incorrect format. Valid format: HeadGroupAbbrev(1'-[Chain1Abbrev/Chain2Abbrev],3'-[Chain3Abbrev/Chain4Abbrev])\n";
    return 0;
  }

  # Make sure 1'- abbreviation is these...
  if ($Abbrev !~ /\(1\'-\[/) {
    warn "Warning: Ignored CL compound abbreviation $Abbrev due to incorrect format. Valid format: HeadGroupAbbrev(1'-[Chain1Abbrev/Chain2Abbrev],3'-[Chain3Abbrev/Chain4Abbrev])\n";
    return 0;
  }

  # Make sure 3'- abbreviation is these...
  if ($Abbrev !~ /\]\,3\'-\[/) {
    warn "Warning: Ignored CL compound abbreviation $Abbrev due to incorrect format. Valid format: HeadGroupAbbrev(1'-[Chain1Abbrev/Chain2Abbrev],3'-[Chain3Abbrev/Chain4Abbrev])\n";
    return 0;
  }

  # Number of '(' and ')' parenthesis must match...
  my(@LeftParenthesis) = $Abbrev =~ /\(/g;
  my(@RightParenthesis) = $Abbrev =~ /\)/g;
  if ($#LeftParenthesis != $#RightParenthesis) {
    warn "Warning: Ignored CL compound abbreviation $Abbrev due to incorrect format: Number of ( and ) doesn't match. Valid format: HeadGroupAbbrev(1'-[Chain1Abbrev/Chain2Abbrev],3'-[Chain3Abbrev/Chain4Abbrev])\n";
    return 0;
  }

  @LeftParenthesis = $Abbrev =~ /\[/g;
  @RightParenthesis = $Abbrev =~ /\]/g;
  if ($#LeftParenthesis != $#RightParenthesis) {
    warn "Warning: Ignored CL compound abbreviation $Abbrev due to incorrect format: Number of [ and ] doesn't match. Valid format: Valid format: HeadGroupAbbrev(1'-[Chain1Abbrev/Chain2Abbrev],3'-[Chain3Abbrev/Chain4Abbrev])\n";
    return 0;
  }

  my($HeadGroup, $ChainsAbbrev, $AbbrevModifier) = ParseCLAbbrev($Abbrev);

  if ($HeadGroup !~ /\*/ ) {
    if (!(exists $CLSupportedHeadGroupMap{$HeadGroup})) {
      warn "Warning: Ignored CL compound abbreviation $Abbrev : Unknown head group $HeadGroup.\n";
      return 0;
    }
  }

  my(@ChainsAbbrevList) = ();
  @ChainsAbbrevList = split /\//, $ChainsAbbrev;
  if (@ChainsAbbrevList != 4) {
    warn "Warning: Ignored CL compound abbreviation $Abbrev due to incorrect format: Didn't specify chain abbreviations for sn1 and sn2 chain positions at 1' and 3' glycerol positions.\n";
    return 0;
  }

  my($Sn1Abbrev, $Sn2Abbrev, $Sn3Abbrev, $Sn4Abbrev) = @ChainsAbbrevList;

  if ($Sn1Abbrev eq "0:0") {
    warn "Warning: Ignoring CL compound abbreviation $Abbrev : sn1 abbreviation value of 0:0 is not allowed.\n";
    return 0;
  }

  if ($Sn3Abbrev eq "0:0") {
    warn "Warning: Ignoring CL compound abbreviation $Abbrev : sn3 abbreviation value of 0:0 is not allowed.\n";
    return 0;
  }

  if ($AbbrevModifier) {
    if (!($AbbrevModifier =~ /\[/ && $AbbrevModifier =~ /\]/)) {
      warn "Warning: Ignored CL compound abbreviation $Abbrev due to incorrect stereochemistry abbreviation: It must be enclosed by []\n";
      return 0;
    }
    $AbbrevModifier =~ s/(\[|\])//g;
    if ($AbbrevModifier !~ /^(rac|U)$/) {
      warn "Warning: Ignored CL compound abbreviation $Abbrev due to incorrect stereochemistry abbreviation: Valid values: rac, U\n";
      return 0;
    }
  }
  return 1;
}

# Setup CL chains abbreviation using individual chain abbreviations...
sub SetupCLChainsAbbrev {
  my($Sn1Abbrev, $Sn2Abbrev, $Sn3Abbrev, $Sn4Abbrev) = @_;
  my($CLChainsAbbrev);

  $CLChainsAbbrev = "1'-[${Sn1Abbrev}/${Sn2Abbrev}],3'-[${Sn3Abbrev}/${Sn4Abbrev}]";

  return $CLChainsAbbrev;
}

# Set up template data for a specific cmpd abbreviation and return a
# reference to the hash...
sub SetupCLCmpdAbbrevTemplateDataMap {
  my($AbbrevHeadGroup, $Abbrev, $ChainsAbbrev, $AbbrevModifier, $Sn1Abbrev, $Sn2Abbrev, $Sn3Abbrev, $Sn4Abbrev, $TemplateID, %AbbrevTemplateDataMap);

  ($Abbrev) = @_;
  %AbbrevTemplateDataMap = ();

  ($AbbrevHeadGroup, $ChainsAbbrev, $AbbrevModifier) = ParseCLAbbrev($Abbrev);

  ($Sn1Abbrev, $Sn2Abbrev, $Sn3Abbrev, $Sn4Abbrev) = split /\//, $ChainsAbbrev;

  my(@SnChainAdd) = ();
  $TemplateID = GetCLTemplateID($AbbrevHeadGroup, $ChainsAbbrev);

  if ($Sn2Abbrev eq "0:0") {
    if ($Sn4Abbrev eq "0:0") {
      push @SnChainAdd, (1,0,1,0);
    }
    else {
      push @SnChainAdd, (1,0,1,1);
    }
  }
  else {
    if ($Sn4Abbrev eq "0:0") {
      push @SnChainAdd, (1,1,1,0);
    }
    else {
      push @SnChainAdd, (1,1,1,1);
    }
  }
  @{$AbbrevTemplateDataMap{SnChainAdd}} = ();
  push @{$AbbrevTemplateDataMap{SnChainAdd}}, @SnChainAdd;

  @{$AbbrevTemplateDataMap{SnAbbrev}} = ();
  push @{$AbbrevTemplateDataMap{SnAbbrev}}, ($Sn1Abbrev, $Sn2Abbrev, $Sn3Abbrev, $Sn4Abbrev);

  $AbbrevModifier =~ s/(\[|\])//g;
  $AbbrevTemplateDataMap{AbbrevModifier} = $AbbrevModifier;

  ChainStr::SetupTemplateDataMap('CL', \%AbbrevTemplateDataMap, $CLTemplatesDataMap{$TemplateID});

  return \%AbbrevTemplateDataMap;
}

# Initialize CL data...
sub _InitializeData {
  _InitializeStrTemplates();
  _InitializeSupportedHeadGroupsData();
}

# Initialize structure template data for these supported templates:
#
sub _InitializeStrTemplates {
  %CLTemplatesDataMap = ();

  # CL templates...

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

  my($CLSn1AlkylSn2Sn3Sn4TemplateString)=<<ENDTEMPLATE;
CL sn1 alkyl, sn2, sn3 and sn4 acyl template structure
  LipdMAPS02060609152D

 36 35  0  0  0  0  0  0  0  0999 V2000
    2.5695    2.3676    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    2.5695    3.1925    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.3567    1.5706    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.3666    2.1547    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.9377    0.9881    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.7266    0.1926    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.1390   -0.5209    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8619    1.4575    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1371    0.2249    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1371    1.0450    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4170    1.4575    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4365    2.1709    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0093    2.1709    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7213    2.5833    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4112    1.4575    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.5237    0.4055    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.9296    0.4055    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.6441   -1.4460    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.8682   -1.0904    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    1.8682   -0.2820    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.6938    2.5364    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5773   -1.8341    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.3927    2.1328    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.7630   -1.2406    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.0486   -0.8281    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6659   -1.2406    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3804   -0.8281    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0948   -1.2406    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.8092   -0.8281    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.5237   -1.2406    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.8092   -0.0031    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0784   -1.9550    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7928   -2.3675    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5073   -1.9550    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7928   -3.1925    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2534   -1.9550    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0      
  1  3  1  0      
  1  4  1  0      
  3  5  1  0      
  5  6  1  0      
  6  7  1  0      
  9 10  2  0      
  8 10  1  0      
 10 11  1  0      
 12 14  1  0      
 13 14  1  0      
  6 16  1  1      
  6 17  1  1      
  7 18  1  0      
 18 19  1  0      
 19 20  2  0      
 19 22  1  0      
 21 13  1  0      
 21 23  1  0      
 23  1  1  0      
 24 19  1  0      
 13 11  1  6      
 13 15  1  1      
 24 25  1  0      
 26 25  1  0      
 27 26  1  0      
 28 27  1  0      
 29 28  1  0      
 30 29  1  0      
 29 31  2  0      
 26 32  1  6      
 33 32  1  0      
 34 33  1  0      
 33 35  2  0      
 26 36  1  1      
M  END

ENDTEMPLATE

  my($CLSn1AlkylSn3Sn4TemplateString)=<<ENDTEMPLATE;
CL sn1 alkyl, sn3 and sn4 acyl template structure
  LipdMAPS02060609152D

 33 32  0  0  0  0  0  0  0  0999 V2000
    2.5697    2.3678    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    2.5697    3.1927    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.3569    1.5707    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.3668    2.1549    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.9379    0.9882    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.7268    0.1926    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.1392   -0.5209    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4170    1.4576    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4366    2.1711    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0093    2.1711    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7214    2.5835    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4112    1.4576    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.5240    0.4055    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.9297    0.4055    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.6443   -1.4461    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.8683   -1.0905    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    1.8683   -0.2820    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.6939    2.5366    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5774   -1.8342    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.3928    2.1330    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.7631   -1.2407    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.0486   -0.8282    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6659   -1.2407    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3805   -0.8282    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0950   -1.2407    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.8094   -0.8282    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.5240   -1.2407    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.8094   -0.0031    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0785   -1.9551    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7929   -2.3677    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5075   -1.9551    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7929   -3.1927    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2534   -1.9551    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0      
  1  3  1  0      
  1  4  1  0      
  3  5  1  0      
  5  6  1  0      
  6  7  1  0      
  9 11  1  0      
 10 11  1  0      
  6 13  1  1      
  6 14  1  1      
  7 15  1  0      
 15 16  1  0      
 16 17  2  0      
 16 19  1  0      
 18 10  1  0      
 18 20  1  0      
 20  1  1  0      
 21 16  1  0      
 10  8  1  6      
 10 12  1  1      
 21 22  1  0      
 23 22  1  0      
 24 23  1  0      
 25 24  1  0      
 26 25  1  0      
 27 26  1  0      
 26 28  2  0      
 23 29  1  6      
 30 29  1  0      
 31 30  1  0      
 30 32  2  0      
 23 33  1  1      
M  END

ENDTEMPLATE

  my($CLSn1Sn2Sn3TemplateString)=<<ENDTEMPLATE;
CL sn1, sn2 and sn3 acyl template structure
  LipdMAPS02060609152D

 36 35  0  0  0  0  0  0  0  0999 V2000
    2.5695    1.6418    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    2.5695    2.4667    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.3567    0.8447    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.3666    1.4289    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.9377    0.2622    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.7266   -0.5332    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.1390   -1.2467    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8619    0.7316    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1371   -0.5009    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1371    0.3191    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4170    0.7316    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.8866    1.4451    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1663    2.6809    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1663    1.8575    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4365    1.4451    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0093    1.4451    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7213    1.8575    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4112    0.7316    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.5237   -0.3203    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.9296   -0.3203    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.6441   -2.1718    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.8682   -1.8162    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    1.8682   -1.0078    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.6938    1.8106    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5773   -2.5600    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.3927    1.4070    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.7630   -1.9664    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.0486   -1.5539    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6659   -1.9664    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3804   -1.5539    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0948   -1.9664    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.8092   -1.5539    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.5237   -1.9664    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.8092   -0.7289    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0784   -2.6809    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2534   -2.6809    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
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
 29 36  1  1      
M  END

ENDTEMPLATE

  my($CLSn1AlkylSn2Sn3TemplateString)=<<ENDTEMPLATE;
CL sn1 alkyl, sn2 and sn3 acyl template structure
  LipdMAPS02060609152D

 33 32  0  0  0  0  0  0  0  0999 V2000
    2.5694    1.7489    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    2.5694    2.5737    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.3567    0.9518    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.3665    1.5360    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.9376    0.3693    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.7265   -0.4261    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.1389   -1.1396    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8619    0.8387    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1371   -0.3938    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1371    0.4262    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4170    0.8387    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4365    1.5522    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0093    1.5522    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7213    1.9646    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4112    0.8387    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.5236   -0.2132    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.9296   -0.2132    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.6440   -2.0647    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.8682   -1.7091    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    1.8682   -0.9007    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.6938    1.9177    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5773   -2.4528    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.3927    1.5141    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.7630   -1.8593    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.0486   -1.4468    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6659   -1.8593    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3804   -1.4468    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0948   -1.8593    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.8091   -1.4468    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.5236   -1.8593    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.8091   -0.6218    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0784   -2.5737    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2534   -2.5737    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0      
  1  3  1  0      
  1  4  1  0      
  3  5  1  0      
  5  6  1  0      
  6  7  1  0      
  9 10  2  0      
  8 10  1  0      
 10 11  1  0      
 12 14  1  0      
 13 14  1  0      
  6 16  1  1      
  6 17  1  1      
  7 18  1  0      
 18 19  1  0      
 19 20  2  0      
 19 22  1  0      
 21 13  1  0      
 21 23  1  0      
 23  1  1  0      
 24 19  1  0      
 13 11  1  6      
 13 15  1  1      
 24 25  1  0      
 26 25  1  0      
 27 26  1  0      
 28 27  1  0      
 29 28  1  0      
 30 29  1  0      
 29 31  2  0      
 26 32  1  6      
 26 33  1  1      
M  END

ENDTEMPLATE

  my($CLSn1Sn3TemplateString)=<<ENDTEMPLATE;
CL sn1 and sn3 acyl template structure
  LipdMAPS02060609152D

 33 32  0  0  0  0  0  0  0  0999 V2000
    2.5694    1.6418    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    2.5694    2.4666    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.3567    0.8447    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.3665    1.4289    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.9376    0.2622    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.7265   -0.5332    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.1389   -1.2467    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4170    0.7316    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.8865    1.4451    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1663    2.6808    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1663    1.8575    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4365    1.4451    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0093    1.4451    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7213    1.8575    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4112    0.7316    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.5236   -0.3203    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.9296   -0.3203    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.6440   -2.1718    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.8682   -1.8162    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    1.8682   -1.0078    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.6938    1.8106    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5773   -2.5599    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.3927    1.4070    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.7630   -1.9664    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.0486   -1.5539    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6659   -1.9664    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3804   -1.5539    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0948   -1.9664    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.8091   -1.5539    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.5236   -1.9664    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.8091   -0.7289    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0784   -2.6808    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2534   -2.6808    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0      
  1  3  1  0      
  1  4  1  0      
  3  5  1  0      
  5  6  1  0      
  6  7  1  0      
 10 11  2  0      
  9 11  1  0      
 11 12  1  0      
 12 14  1  0      
 13 14  1  0      
  6 16  1  1      
  6 17  1  1      
  7 18  1  0      
 18 19  1  0      
 19 20  2  0      
 19 22  1  0      
 21 13  1  0      
 21 23  1  0      
 23  1  1  0      
 24 19  1  0      
 13  8  1  6      
 13 15  1  1      
 24 25  1  0      
 26 25  1  0      
 27 26  1  0      
 28 27  1  0      
 29 28  1  0      
 30 29  1  0      
 29 31  2  0      
 26 32  1  6      
 26 33  1  1      
M  END

ENDTEMPLATE

  my($CLSn1AlkylSn3TemplateString)=<<ENDTEMPLATE;
CL sn1 alkyl and sn3 acyl template structure
  LipdMAPS02060609152D

 30 29  0  0  0  0  0  0  0  0999 V2000
    2.5696    1.7490    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    2.5696    2.5739    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.3569    0.9519    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.3667    1.5361    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.9378    0.3693    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.7267   -0.4261    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.1391   -1.1397    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4170    0.8388    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4366    1.5523    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0093    1.5523    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7214    1.9647    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4112    0.8388    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.5239   -0.2132    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.9297   -0.2132    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.6442   -2.0649    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.8683   -1.7092    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    1.8683   -0.9008    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.6939    1.9178    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5774   -2.4530    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.3928    1.5142    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.7631   -1.8594    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.0486   -1.4469    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6659   -1.8594    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3805   -1.4469    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0950   -1.8594    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.8093   -1.4469    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.5239   -1.8594    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.8093   -0.6218    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0785   -2.5739    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2534   -2.5739    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0      
  1  3  1  0      
  1  4  1  0      
  3  5  1  0      
  5  6  1  0      
  6  7  1  0      
  9 11  1  0      
 10 11  1  0      
  6 13  1  1      
  6 14  1  1      
  7 15  1  0      
 15 16  1  0      
 16 17  2  0      
 16 19  1  0      
 18 10  1  0      
 18 20  1  0      
 20  1  1  0      
 21 16  1  0      
 10  8  1  6      
 10 12  1  1      
 21 22  1  0      
 23 22  1  0      
 24 23  1  0      
 25 24  1  0      
 26 25  1  0      
 27 26  1  0      
 26 28  2  0      
 23 29  1  6      
 23 30  1  1      
M  END

ENDTEMPLATE


# Format: ID => AbbrevID|HeadGroupName|Y1Sn1|Y2Sn1|Y1Sn2|Y2Sn2|Y1Sn3|Y2Sn3|Y1Sn4|Y2Sn4|Sn1ChainAtomNum|Sn2ChainAtomNum|Sn3ChainAtomNum|Sn4ChainAtomNum|Sn1ChainCarbons|Sn2ChainCarbons|Sn3ChainCarbons|Sn4ChainCarbons|Sn2CAtomNum|Sn2OAtomNum|Sn2HAtomNum|LMCategory|LMMainClass|LMSubClass|TemplateString
  # Note: Template structure for Both Alkyl and Alkenyl substituted CL is the same; Different entries are used to facilitate assigning LM identity
  %CLTemplatesDataMap = (
			       "CLSn1Sn2Sn3Sn4" => "CL|glycerol|2.0639|2.4763|0.9379|1.3504|-1.3477|-0.9352|-2.4747|-2.0622|12|8|33|37|2|2|2|2|6|20|19|GP|12|01|$CLSn1Sn2Sn3Sn4TemplateString",

			       "CLSn1AlkylSn2Sn3Sn4" => "CL|glycerol|2.1709|2.5833|1.0450|1.4575|-1.2406|-0.8281|-2.3675|-1.9550|12|8|30|34|0|2|2|2|6|17|16|GP|12|03|$CLSn1AlkylSn2Sn3Sn4TemplateString",
			       "CLSn1AlkenylSn2Sn3Sn4" => "CL|glycerol|2.1709|2.5833|1.0450|1.4575|-1.2406|-0.8281|-2.3675|-1.9550|12|8|30|34|0|2|2|2|6|17|16|GP|12|05|$CLSn1AlkylSn2Sn3Sn4TemplateString",

			       "CLSn1AlkylSn3Sn4" => "CL|glycerol|2.1711|2.5835|0|0|-1.2407|-0.8282|-2.3677|-1.9551|9|0|27|31|0|0|2|2|6|14|13|GP|12|08|$CLSn1AlkylSn3Sn4TemplateString",
			       "CLSn1AlkenylSn3Sn4" => "CL|glycerol|2.1711|2.5835|0|0|-1.2407|-0.8282|-2.3677|-1.9551|9|0|27|31|0|0|2|2|6|14|13|GP|12|10|$CLSn1AlkylSn3Sn4TemplateString",

			       "CLSn1Sn2Sn3" => "CL|glycerol|1.4451|1.8575|0.3191|0.7316|-1.9664|-1.5539|0|0|12|8|33|0|2|2|2|0|6|20|19|GP|12|02|$CLSn1Sn2Sn3TemplateString",

			       "CLSn1AlkylSn2Sn3" => "CL|glycerol|1.5522|1.9646|0.4262|0.8387|-1.8593|-1.4468|0|0|12|8|30|0|0|2|2|0|6|17|16|GP|12|04|$CLSn1AlkylSn2Sn3TemplateString",
			       "CLSn1AlkenylSn2Sn3" => "CL|glycerol|1.5522|1.9646|0.4262|0.8387|-1.8593|-1.4468|0|0|12|8|30|0|0|2|2|0|6|17|16|GP|12|06|$CLSn1AlkylSn2Sn3TemplateString",

			       "CLSn1Sn3" => "CL|glycerol|1.4451|1.8575|0|0|-1.9664|-1.5539|0|0|9|0|30|0|2|0|2|0|6|17|16|GP|12|07|$CLSn1Sn3TemplateString",
			       "CLSn1AlkylSn3" => "CL|glycerol|1.5523|1.9647|0|0|-1.8594|-1.4469|0|0|9|0|27|0|0|0|2|0|6|14|13|GP|12|09|$CLSn1AlkylSn3TemplateString",
			       "CLSn1AlkenylSn3" => "CL|glycerol|1.5523|1.9647|0|0|-1.8594|-1.4469|0|0|9|0|27|0|0|0|2|0|6|14|13|GP|12|11|$CLSn1AlkylSn3TemplateString",
			       );
}

# Initialize supported head groups...
sub _InitializeSupportedHeadGroupsData {
  my($GPType, $GPHeadGroup);
  %CLSupportedHeadGroupMap = ();

  for $GPType (keys %CLTemplatesDataMap) {
    ($GPHeadGroup) = split /\|/, $CLTemplatesDataMap{$GPType};
    if (!(exists $CLSupportedHeadGroupMap{$GPHeadGroup})) {
      $CLSupportedHeadGroupMap{$GPHeadGroup} = $GPHeadGroup;
    }
  }
}


1;

__END__

=head1 NAME

CLStr - Cardiolipins (CL) structure generation methods

=head1 SYNOPSIS

use CLStr;

use CLStr qw(:all);

=head1 DESCRIPTION

CLStr module provides these methods:

    GenerateCmpdOntologyData - Generate ontology data
    GenerateCmpdOntologySDDataLines - Generate ontology data lines for
                                      SD file
    GenerateCLChainStrData - Generate chain structure data
    GenerateSDFile - Generate SD file
    GetCLTemplatesData - Get templates data
    GetCLSupportedHeadGroupMap - Get supported headgroups data
    GetCLTemplateID - Get templates ID
    IsCLChainsAbbrevSupported - Is it a supported CL abbreviation
    ParseCLAbbrev - Parse CL abbreviation
    ProcessCLCmpdAbbrevs - Process CL abbreviation
    SetupCLCmpdAbbrevTemplateDataMap - Setup template structure data map
    ValidateCLAbbrev - Validate CL abbreviation

=head1 METHODS

=over 4

=item B<GenerateCmpdOntologyData>

    $DataHashRef = GenerateCmpdOntologyData($CmpDataRef);

Return a reference to a hash containing ontology data with hash keys and values
corresponding to property names and values.

=item B<GenerateCmpdOntologySDDataLines>

    $DataLinesArrayRef =
       GenerateCmpdOntologySDDataLines($CmpdDataRef);

Return a reference to an array containing ontology data lines suitable for
generate SD file data block.

=item B<GenerateCLChainStrData>

    ($AtomLinesArrayRef, $BondLinesArrayRef) =
       GenerateCLChainStrData($ChainType, $CmpdDataRef);

Return array references containing atom and bond data lines for SD file. Appropriate atom
and bond data lines are generated using chain type and abbreviation template data.

=item B<GenerateSDFile>

    GenerateSDFile($SDFileName, $CmdAbbrevsRef);

Generate a SD file for compound abbreviations. Structure data for specified abbreviation
is generated sequentially and written to SD file.

=item B<GetCLTemplatesData>

    $TemplatesDataRef = GetCLTemplatesData();

Return a reference to a hash containing CL templates data

=item B<GetCLSupportedHeadGroupMap>

    $SupportedHeadGroupDataRef = GetCLSupportedHeadGroupMap();

Return a reference to a hash containing supported head groups data.

=item B<GetCLTemplateID>

    $HeadGroupID = GetCLTemplateID($HeadGroupAbbrev, $ChainsAbbrev);

Return a supported template ID for compound abbreviation.

=item B<IsCLChainsAbbrevSupported>

    $Status = IsCLChainsAbbrevSupported($Abbrev, $PrintWarning);

Return 1 or 0 based on whether CL abbreviated is supported. For unsupported CL abbreviations,
a warning is printed unless PrintWarning flag is set.

=item B<ParseCLAbbrev>

    ($HeadGroup, $ChainsAbbrev, $AbbrevModifier) =
       ParseCLAbbrev($Abbrev);

Parse CL abbreviation and return these values: HeadGroup, ChainsAbbrev,
AbbrevModifier.

=item B<ProcessCLCmpdAbbrevs>

    ProcessCLCmpdAbbrevs($CmpdAbbrevsRef, $AllowArbitraryChainAbbrev,
                         $WriteSDFile, $SDFileName);

Process specified CL abbreviations to generate structures and write them out either
a SD file or simply report number of valid abbreviations.

=item B<SetupCLCmpdAbbrevTemplateDataMap>

    $AbbrevTemplateDataMapRef =
       SetupCLCmpdAbbrevTemplateDataMap($Abbrev);

Return a reference to a hash containing template data for compound abbreviation. The
template data is used to generate SD file for compound abbreviation.

=item B<ValidateCLAbbrev>

    $Status = ValidateCLAbbrev($Abbrev);

Return 1 or 0 based on whether a CL abbreviation is valid.

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
