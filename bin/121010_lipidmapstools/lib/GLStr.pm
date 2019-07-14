package GLStr;
#
# File: GLStr.pm
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
use LMAPSStr;
use ChainAbbrev;
use ChainStr;

our ($VERSION, @ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);

$VERSION = '2.00';
@ISA = qw(Exporter);
@EXPORT = qw();
@EXPORT_OK = qw(GenerateCmpdOntologyData GenerateCmpdOntologyDataLines GenerateGLChainStrData GenerateSDFile IsGLAbbrevSupported ParseGLAbrev ProcessGLCmpdAbbrevs SetupGLCmpdAbbrevTemplateDataMap ValidateGLAbbrev);
%EXPORT_TAGS = (
		all  => [@EXPORT, @EXPORT_OK],
		);

# Initialize GL data...
my(%GLTemplatesDataMap);
_InitializeData();

# Process specified compound abbreviations containing any wild cards
# and count the number of valid abbreviations along with optional generation
# of a SD file.
#
sub ProcessGLCmpdAbbrevs {
  my($CmpdAbbrevsRef, $AllowArbitraryChainAbbrev, $WriteSDFile, $SDFileName) = @_;

  $AllowArbitraryChainAbbrev = defined $AllowArbitraryChainAbbrev ? $AllowArbitraryChainAbbrev : 0;

  $WriteSDFile = defined $WriteSDFile ? $WriteSDFile : 0;
  if ($WriteSDFile &&  IsEmpty($SDFileName)) {
    warn "Warning: GLStr::ProcessGLCmpdAbbrevs: No SD file name specified. Ingoring structure generation.\n";
    return;
  }

  if ($WriteSDFile) {
    print "Generating SD file $SDFileName...\n";

    open SDFILE, ">$SDFileName" or die "Error: Couldn't open $SDFileName: $! \n";

    _ProcessGLCmpdAbbrevs($CmpdAbbrevsRef, $AllowArbitraryChainAbbrev, $WriteSDFile, \*SDFILE);

    close SDFILE;
  }
  else {
    _ProcessGLCmpdAbbrevs($CmpdAbbrevsRef, $AllowArbitraryChainAbbrev, $WriteSDFile);
  }
}

# Process specified compound abbreviations containing any wild cards
# and count the number of valid abbreviations along with optional generation
# of a SD file.
#
sub _ProcessGLCmpdAbbrevs {
  my($CmpdAbbrevsRef, $AllowArbitraryChainAbbrev, $WriteSDFile, $SDFileRef) = @_;
  my($Abbrev, $TemplateID, $AbbrevCount, $GLIsomersAbbrev, $Sn1Abbrev, $Sn2Abbrev, $Sn3Abbrev, $GLType, $ChainsAbbrev, $AbbrevModifier, $AllowSubstituents, $AllowRings);

  $AbbrevCount = 0;

  ABBREV: for $Abbrev (@{$CmpdAbbrevsRef}) {
    $Sn1Abbrev = "0:0"; $Sn2Abbrev = "0:0"; $Sn3Abbrev = "0:0";

    if (!ValidateGLAbbrev($Abbrev)) {
      next ABBREV;
    }

    ($GLType, $ChainsAbbrev, $AbbrevModifier) = ParseGLAbbrev($Abbrev);
    ($Sn1Abbrev, $Sn2Abbrev, $Sn3Abbrev) = split /\//, $ChainsAbbrev;

    ($AllowSubstituents, $AllowRings) = (1, 0);
    if (!(ChainAbbrev::IsChainAbbrevOkay($Sn1Abbrev, $AllowSubstituents, $AllowRings, $AllowArbitraryChainAbbrev) && ChainAbbrev::IsChainAbbrevOkay($Sn2Abbrev, $AllowSubstituents, $AllowRings, $AllowArbitraryChainAbbrev) && ChainAbbrev::IsChainAbbrevOkay($Sn3Abbrev, $AllowSubstituents, $AllowRings, $AllowArbitraryChainAbbrev))) {
      warn "Warning: Ignoring GL compound abbreviation $Abbrev\n";
      next ABBREV;
    }
    if (!IsGLAbbrevSupported($Abbrev, 1)) {
      warn "Warning: Ignoring compound abbreviation $Abbrev\n";
      next ABBREV;
    }

    my($WildCardInGLType) = ($GLType =~ /\*/ ) ? 1 : 0;;

    if (!($WildCardInGLType || ChainAbbrev::IsWildCardInChainAbbrev($Sn1Abbrev) || ChainAbbrev::IsWildCardInChainAbbrev($Sn2Abbrev) || ChainAbbrev::IsWildCardInChainAbbrev($Sn3Abbrev))) {
      if ($AbbrevModifier =~ /iso/) {
	my($IsomerGLType, $IsomerChainsAbbrev, $IsomerAbbrevModifier);

	for $GLIsomersAbbrev (ExpandGLIsomersAbbrev($Abbrev)) {
	  ($IsomerGLType, $IsomerChainsAbbrev, $IsomerAbbrevModifier) = ParseGLAbbrev($Abbrev);
	  $TemplateID = _GetGLTemplateID($IsomerChainsAbbrev);

	  if (!exists $GLTemplatesDataMap{$TemplateID}) {
	    warn "Warning: Ignored GL isomer compound abbreviation $GLIsomersAbbrev : Abbreviation doesn't match any template\n";
	    next ABBREV;
	  }

	  $AbbrevCount++;
	  if ($WriteSDFile) {
	    _GenerateAndWriteCmdDataString($SDFileRef, $GLIsomersAbbrev);
	  }
	}
      }
      else {
	$TemplateID = _GetGLTemplateID($ChainsAbbrev);
	if (!exists $GLTemplatesDataMap{$TemplateID}) {
	  warn "Warning: Ignored GL compound abbreviation $Abbrev : Abbreviation doesn't match any template\n";
	  next ABBREV;
	}
	$AbbrevCount++;
	if ($WriteSDFile) {
	  _GenerateAndWriteCmdDataString($SDFileRef, $Abbrev);
	}
      }
      next ABBREV;
    }

    # Abbreviation modifier is not supported for with wild cards...
    if ($AbbrevModifier =~ /iso/i) {
      warn "Warning: Ignoring GL compound abbreviation $Abbrev : isomer abbreviation, $AbbrevModifier, is not supported with wild cards in any part of the abbreviation\n";
      next ABBREV;
    }

    # Arbitrary acyl chain abbreviation is not supported with wild cards...
    if ($AllowArbitraryChainAbbrev) {
      warn "Warning: Ignoring GL compound abbreviation $Abbrev : Allow arbitrary chain abbreviation option is not supported with wild cards in any part of the abbreviation\n";
      next ABBREV;
    }

    # Substituents are not supported with wild cards...
    if (ChainAbbrev::IsSubstituentInChainAbbrev($Sn1Abbrev) || ChainAbbrev::IsSubstituentInChainAbbrev($Sn2Abbrev) || ChainAbbrev::IsSubstituentInChainAbbrev($Sn3Abbrev)) {
      warn "Warning: Ignoring GL compound abbreviation $Abbrev : Substituent specifications are not supported with wild cards in any part of the abbreviation\n";
      next ABBREV;
    }

    # Get expanded abbreviation for each position...
    my($Sn1ExpandedAbbrevRef) = ChainAbbrev::ExpandChainAbbrev($Sn1Abbrev);
    my($Sn2ExpandedAbbrevRef) = ChainAbbrev::ExpandChainAbbrev($Sn2Abbrev);
    my($Sn3ExpandedAbbrevRef) = ChainAbbrev::ExpandChainAbbrev($Sn3Abbrev);

    # Enumerate various possibilities...
    my($GLTypeAbbrev, $ExpandedAbbrev, $ExpandedSn1Abbrev, $ExpandedSn2Abbrev, $ExpandedSn3Abbrev);
    $GLTypeAbbrev = $GLType;

    if ($WildCardInGLType) {
      # Enumerate MG, DG and TG...
      my(@GLTypeList) = ();
      push @GLTypeList, ('MG', 'DG', 'TG');
      for $GLTypeAbbrev (@GLTypeList) {
	SN1ABBREV: for $ExpandedSn1Abbrev (@$Sn1ExpandedAbbrevRef) {
	  SN2ABBREV: for $ExpandedSn2Abbrev (@$Sn2ExpandedAbbrevRef) {
	    if ($ExpandedSn2Abbrev =~ /^(O-|P-)/) {
	      next SN2ABBREV;
	    }
	    SN3ABBREV: for $ExpandedSn3Abbrev (@$Sn3ExpandedAbbrevRef) {
	      if ($ExpandedSn3Abbrev =~ /^(O-|P-)/) {
		next SN3ABBREV;
	      }
	      if ($GLTypeAbbrev eq 'MG') {
		$ExpandedAbbrev = $GLTypeAbbrev . '(' . $ExpandedSn1Abbrev . '/0:0/0:0)';
	      }
	      elsif ($GLTypeAbbrev eq 'DG') {
		$ExpandedAbbrev = $GLTypeAbbrev . '(' . $ExpandedSn1Abbrev . '/' . $ExpandedSn2Abbrev . '/0:0)' . $AbbrevModifier;
	      }
	      else {
		$ExpandedAbbrev = $GLTypeAbbrev . '(' . $ExpandedSn1Abbrev . '/' . $ExpandedSn2Abbrev . '/' . $ExpandedSn3Abbrev . ')' . $AbbrevModifier;
	      }

	      $AbbrevCount++;
	      if ($WriteSDFile) {
		_GenerateAndWriteCmdDataString($SDFileRef, $ExpandedAbbrev);
	      }

	      if ($GLTypeAbbrev =~ /^MG$/i) {
		next SN1ABBREV;
	      }
	      elsif ($GLTypeAbbrev =~ /^DG$/i) {
		next SN2ABBREV;
	      }
	    }
	  }
	}
      }
    }
    else {
      for $ExpandedSn1Abbrev (@$Sn1ExpandedAbbrevRef) {
	SN2ABBREV: for $ExpandedSn2Abbrev (@$Sn2ExpandedAbbrevRef) {
	  if ($ExpandedSn2Abbrev =~ /^(O-|P-)/) {
	    next SN2ABBREV;
	  }
	  SN3ABBREV: for $ExpandedSn3Abbrev (@$Sn3ExpandedAbbrevRef) {
	    if ($ExpandedSn3Abbrev =~ /^(O-|P-)/) {
	      next SN3ABBREV;
	    }
	    $ExpandedAbbrev = $GLTypeAbbrev . '(' . $ExpandedSn1Abbrev . '/' . $ExpandedSn2Abbrev . '/' . $ExpandedSn3Abbrev . ')' . $AbbrevModifier;

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

# Exapand isomers abbreviation...
sub ExpandGLIsomersAbbrev {
  my($Abbrev) = @_;
  my(@ExpandedIsomersAbbrevs) = ();

  my($GLType, $ChainsAbbrev, $AbbrevModifier) = ParseGLAbbrev($Abbrev);
  my($Sn1Abbrev, $Sn2Abbrev, $Sn3Abbrev) = split /\//, $ChainsAbbrev;
  my($IsomerAbbrev) = $AbbrevModifier;
  $IsomerAbbrev =~ s/(\[|\])//g;

  if ($IsomerAbbrev =~ /^iso2$/) {
    push @ExpandedIsomersAbbrevs, "${GLType}(${Sn1Abbrev}/${Sn2Abbrev}/${Sn3Abbrev})";
    push @ExpandedIsomersAbbrevs, "${GLType}(${Sn2Abbrev}/${Sn1Abbrev}/${Sn3Abbrev})";
  }
  elsif ($IsomerAbbrev =~ /^iso3$/) {
    my(@ChainsAbbrevList) = ();
      # Setup similar chains at first two positions...
      if ($Sn1Abbrev eq $Sn2Abbrev) {
	push @ChainsAbbrevList, ($Sn1Abbrev, $Sn1Abbrev, $Sn3Abbrev);
      }
      elsif ($Sn1Abbrev eq $Sn3Abbrev) {
	push @ChainsAbbrevList, ($Sn1Abbrev, $Sn3Abbrev, $Sn2Abbrev);
      }
      elsif ($Sn2Abbrev eq $Sn3Abbrev) {
	push @ChainsAbbrevList, ($Sn2Abbrev, $Sn2Abbrev, $Sn1Abbrev);
      }
    push @ExpandedIsomersAbbrevs, "${GLType}($ChainsAbbrevList[0]/$ChainsAbbrevList[1]/$ChainsAbbrevList[2])";
    push @ExpandedIsomersAbbrevs, "${GLType}($ChainsAbbrevList[0]/$ChainsAbbrevList[2]/$ChainsAbbrevList[1])";
    push @ExpandedIsomersAbbrevs, "${GLType}($ChainsAbbrevList[2]/$ChainsAbbrevList[0]/$ChainsAbbrevList[1])";
  }
  elsif ($IsomerAbbrev =~ /^iso6$/) {
    my(@ChainsAbbrevList) = ();
    push @ChainsAbbrevList, ($Sn1Abbrev, $Sn2Abbrev, $Sn3Abbrev);
    my($Chain1Abbrev, $Chain2Abbrev, $Chain3Abbrev);
    for $Chain1Abbrev (@ChainsAbbrevList) {
      CHAIN2: for $Chain2Abbrev (@ChainsAbbrevList) {
	if ($Chain1Abbrev eq $Chain2Abbrev) {
	  next CHAIN2;
	}
	CHAIN3: for $Chain3Abbrev (@ChainsAbbrevList) {
	  if (($Chain3Abbrev eq $Chain1Abbrev) || ($Chain3Abbrev eq $Chain2Abbrev) ) {
	    next CHAIN3;
	  }
	  push @ExpandedIsomersAbbrevs, "${GLType}(${Chain1Abbrev}/${Chain2Abbrev}/${Chain3Abbrev})";
	}
      }
    }
  }

  return @ExpandedIsomersAbbrevs;
}

# Generate a SD file for valid abbreviations...
sub GenerateSDFile {
  my($SDFileName, $CmdAbbrevsRef) = @_;
  my($Abbrev, $CmpdDataString);

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
  $CmpdAbbrevTemplateDataMapRef = SetupGLCmpdAbbrevTemplateDataMap($Abbrev);

  # Generate structure data for chains...
  my($Sn1AtomBlockLinesRef, $Sn1BondBlockLinesRef) = GLStr::GenerateGLChainStrData('Sn1', $CmpdAbbrevTemplateDataMapRef);
  my($Sn2AtomBlockLinesRef, $Sn2BondBlockLinesRef) = GLStr::GenerateGLChainStrData('Sn2', $CmpdAbbrevTemplateDataMapRef);
  my($Sn3AtomBlockLinesRef, $Sn3BondBlockLinesRef) = GLStr::GenerateGLChainStrData('Sn3', $CmpdAbbrevTemplateDataMapRef);

  # Generate data block lines including various desriptors...
  my($OntologyDataLinesRef) = GenerateCmpdOntologySDDataLines($CmpdAbbrevTemplateDataMapRef);

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
sub GenerateGLChainStrData {
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

# Multiple bond counts...
sub _SetupChainLengthAndMultipleBondsCountInfo {
  my($CmpdAbbrevTemplateDataMapRef, $OntologyDataMapRef) = @_;
  my($ChainIndex, $ChainAbbrev, $ChainCount, $ChainLengthAbbrev, $DoubleBondCountAbbrev, $ChainLength, $DoubleBonds, @ChainIndices);

  @ChainIndices = (0, 1, 2);

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

# Abbreviation and systematic name...
sub _SetupAbbrevAndSystematicName {
  my($CmpdAbbrevTemplateDataMapRef, $OntologyDataMapRef) = @_;
  my($ChainsSystematicName, $GLTypeAbbrev, $AbbrevModifier, $ChainsAbbrev, @ChainIndices);

  @ChainIndices = (0, 1, 2);

  $ChainsSystematicName = ChainAbbrev::SetupChainsSystematicName($CmpdAbbrevTemplateDataMapRef, \@ChainIndices);
  $OntologyDataMapRef->{'Systematic Name'} = "${ChainsSystematicName}-sn-glycerol";

  $GLTypeAbbrev = $CmpdAbbrevTemplateDataMapRef->{AbbrevID};
  $AbbrevModifier = $CmpdAbbrevTemplateDataMapRef->{AbbrevModifier};
  if ($AbbrevModifier) {
    $AbbrevModifier = "[$AbbrevModifier]";
  }
  $ChainsAbbrev = ChainAbbrev::SetupChainsAbbrev($CmpdAbbrevTemplateDataMapRef, \@ChainIndices);

  $OntologyDataMapRef->{Abbrev} = "${GLTypeAbbrev}(${ChainsAbbrev})${AbbrevModifier}";
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

# Parse abbreviation...
sub ParseGLAbbrev {
  my($Abbrev) = @_;

  my($GLType, $ChainsAbbrev, $AbbrevModifier);
  if ($Abbrev =~ /\]$/) {
    ($GLType, $ChainsAbbrev, $AbbrevModifier) = $Abbrev =~ /^(.+?)\((.+?)\)\[(.+?)\]$/;
    $AbbrevModifier = '[' . $AbbrevModifier . ']';
  }
  else {
    ($GLType, $ChainsAbbrev) = $Abbrev =~ /^(.+?)\((.+?)\)$/;
    $AbbrevModifier = '';
  }
  return ($GLType, $ChainsAbbrev, $AbbrevModifier);
}

# Does template exist to handle this abbreviation?
#
# At this point in time, alkyl/alkenyl chains at sn2/sn3 are not supported for GL: Templates need to be created.
#
# For GL: alkenyl at sn2/sn3; alkyl at sn2 without alkyl at sn1
#
sub IsGLAbbrevSupported {
  my($Abbrev, $PrintWarning) = @_;
  my($GLType, $ChainsAbbrev, $AbbrevModifier, $Sn1Abbrev, $Sn2Abbrev, $Sn3Abbrev, @AbbrevList);

  ($GLType, $ChainsAbbrev, $AbbrevModifier) = ParseGLAbbrev($Abbrev);

  $Sn1Abbrev = "0:0"; $Sn2Abbrev = "0:0"; $Sn3Abbrev = "0:0";
  @AbbrevList = split /\//, $ChainsAbbrev;
  if (@AbbrevList == 1) {
    ($Sn1Abbrev) = @AbbrevList;
  }
  if (@AbbrevList == 2) {
    ($Sn1Abbrev, $Sn2Abbrev) = @AbbrevList;
  }
  if (@AbbrevList == 3) {
    ($Sn1Abbrev, $Sn2Abbrev, $Sn3Abbrev) = @AbbrevList;
  }

  if (ChainAbbrev::IsAlkenylChainAbbrev($Sn2Abbrev)) {
    if ($PrintWarning) {
      warn "Warning: Ignoring GL compound abbreviation $Abbrev : alkenyl chain at sn2 position.\n";
    }
    return 0;
  }
  if (ChainAbbrev::IsAlkylChainAbbrev($Sn2Abbrev) && !ChainAbbrev::IsAlkylChainAbbrev($Sn1Abbrev)) {
    if ($PrintWarning) {
      warn "Warning: Ignoring GL compound abbreviation $Abbrev : alkyl chain at sn2 position.\n";
    }
    return 0;
  }
  if (ChainAbbrev::IsAlkylChainAbbrev($Sn3Abbrev) || ChainAbbrev::IsAlkenylChainAbbrev($Sn3Abbrev)) {
    if ($PrintWarning) {
      warn "Warning: Ignoring GL compound abbreviation $Abbrev : alkyl/alkenyl chain at sn3 position.\n";
    }
    return 0;
  }

  return 1;
}

#
# Check out the validity of GL abbreviation...
sub ValidateGLAbbrev {
  my($Abbrev) = @_;

  if (!($Abbrev =~ /\(/ && $Abbrev =~ /\)/)) {
    warn "Warning: Ignored GL compound abbreviation $Abbrev due to incorrect format. Valid format: GLTypeAbbrev(Chain1Abbrev/Chain2Abbrev/Chain3Abbrev).\n";
    return 0;
  }

  # Make sure GL type is these...
  if ($Abbrev =~ /^\(/) {
    warn "Warning: Ignored GL compound abbreviation $Abbrev due to incorrect format. Valid format: GLTypeAbbrev(Chain1Abbrev/Chain2Abbrev/Chain3Abbrev).\n";
    return 0;
  }

  # Number of '(' and ')' parenthesis must match...
  my(@LeftParenthesis) = $Abbrev =~ /\(/g;
  my(@RightParenthesis) = $Abbrev =~ /\)/g;
  if ($#LeftParenthesis != $#RightParenthesis) {
    warn "Warning: Ignored GL compound abbreviation $Abbrev due to incorrect format: Number of ( and ) doesn't match. Valid format: GLTypeAbbrev(Chain1Abbrev/Chain2Abbrev/Chain3Abbrev).\n";
    return 0;
  }

  @LeftParenthesis = $Abbrev =~ /\[/g;
  @RightParenthesis = $Abbrev =~ /\]/g;
  if ($#LeftParenthesis != $#RightParenthesis) {
    warn "Warning: Ignored GL compound abbreviation $Abbrev due to incorrect format: Number of [ and ] doesn't match. Valid format: GLTypeAbbrev(Chain1Abbrev/Chain2Abbrev/Chain3Abbrev).\n";
    return 0;
  }

  my($GLType, $ChainsAbbrev, $AbbrevModifier) = ParseGLAbbrev($Abbrev);

  if ($GLType !~ /^(MG|DG|TG|\*)$/) {
    warn "Warning: Ignored GL compound abbreviation $Abbrev due to incorrect GL type abbreviation. Valid GL type abbreviations: MG, DG or TG. \n";
    return 0;
  }

  my(@ChainsAbbrevList) = ();
  @ChainsAbbrevList = split /\//, $ChainsAbbrev;
  if (@ChainsAbbrevList != 3) {
    warn "Warning: Ignored GL compound abbreviation $Abbrev due to incorrect format: Didn't specify abbreviations for all three sn1, sn2, and sn3 chains.\n";
    return 0;
  }

  my($Sn1Abbrev, $Sn2Abbrev, $Sn3Abbrev) = @ChainsAbbrevList;

  if ($Sn1Abbrev eq "0:0") {
    warn "Warning: Ignoring GL compound abbreviation $Abbrev : sn1 abbreviation value of 0:0 is not allowed.\n";
    return 0;
  }

  if ($GLType =~ /^MG$/) {
    if ($Sn2Abbrev !~ /^0:0$/ || $Sn3Abbrev !~ /^0:0$/) {
      warn "Warning: Ignoring GL compound abbreviation $Abbrev : sn2 and sn3 abbreviations must be 0:0 for MG\n";
      return 0;
    }
  }
  elsif ($GLType =~ /^DG$/) {
    if ($Sn2Abbrev =~ /^0:0$/ && $Sn3Abbrev =~ /^0:0$/) {
      warn "Warning: Ignoring GL compound abbreviation $Abbrev : one of sn2 and sn3 abbreviation must not be 0:0 for DG\n";
      return 0;
    }
  }
  elsif ($GLType =~ /^TG$/) {
    if ($Sn2Abbrev =~ /^0:0$/ || $Sn3Abbrev =~ /^0:0$/) {
      warn "Warning: Ignoring GL compound abbreviation $Abbrev : sn2 and sn3 abbreviations must not be 0:0 for TG\n";
      return 0;
    }
  }

  if ($AbbrevModifier) {
    if (!($AbbrevModifier =~ /\[/ && $AbbrevModifier =~ /\]/)) {
      warn "Warning: Ignored GL compound abbreviation $Abbrev due to incorrect isomer/sn2 stereochemistry abbreviation: It must be enclosed by []\n";
      return 0;
    }
    $AbbrevModifier =~ s/(\[|\])//g;
    if ($AbbrevModifier !~ /^(rac|R|S|U|iso2|iso3|iso6)$/) {
      warn "Warning: Ignored GL compound abbreviation $Abbrev due to incorrect isomer or sn2 stereochemistry abbreviation: Valid values: rac, R, S, U, iso2, iso3, or iso6 \n";
      return 0;
    }
    # Makes sure for isomers specification, appropriate chains abbrev are specified...
    if ($GLType =~ /^MG$/ && $AbbrevModifier =~ /iso/) {
      warn "Warning: Ignoring GL compound abbreviation $Abbrev : isomer abbreviation [$AbbrevModifier] is not allowed for MG\n";
      return 0;
    }
    if ($GLType =~ /^DG$/ && $AbbrevModifier =~ /^(iso3|iso6)$/) {
      warn "Warning: Ignoring GL compound abbreviation $Abbrev : isomer abbreviation [$AbbrevModifier] is not valid for DG\n";
      return 0;
    }
    if ($GLType =~ /^TG$/ && $AbbrevModifier =~ /^iso2$/) {
      warn "Warning: Ignoring GL compound abbreviation $Abbrev : isomer abbreviation [$AbbrevModifier] is not valid for TG\n";
      return 0;
    }

    if ($AbbrevModifier =~ /^iso2$/) {
      if ($Sn1Abbrev eq $Sn2Abbrev) {
	warn "Warning: Ignoring GL compound abbreviation $Abbrev : isomer abbreviation [iso2] is not valid for same chain abbrev at sn1 and sn2 positions\n";
	return 0;
      }
    }
    elsif ($AbbrevModifier =~ /^iso3$/) {
      if (!(($Sn1Abbrev eq $Sn2Abbrev && $Sn1Abbrev ne $Sn3Abbrev) || ($Sn1Abbrev ne $Sn2Abbrev && $Sn1Abbrev eq $Sn3Abbrev) || ($Sn1Abbrev ne $Sn2Abbrev && $Sn2Abbrev eq $Sn3Abbrev) )) {
	warn "Warning: Ignoring GL compound abbreviation $Abbrev : isomer abbreviation [iso3] is not valid for similar or different chain abbrev at all sn1, sn2, and sn3 positions\n";
	return 0;
      }
    }
    elsif ($AbbrevModifier =~ /^iso6$/) {
      if (!($Sn1Abbrev ne $Sn2Abbrev && $Sn1Abbrev ne $Sn3Abbrev && $Sn2Abbrev ne $Sn3Abbrev)) {
	warn "Warning: Ignoring GL compound abbreviation $Abbrev : isomer abbreviation [iso6] is not valid for similar chain abbrev at any sn1, sn2, and sn3 positions\n";
	return 0;
      }
    }
  }
  return 1;
}

# Set up template data for a specific cmpd abbreviation and return a
# reference to the hash...
sub SetupGLCmpdAbbrevTemplateDataMap {
  my($Abbrev, $GLType, $ChainsAbbrev, $AbbrevModifier, $SnAbbrev, $TemplateID, @SnAbbrevs, @SnChainAdd, %AbbrevTemplateDataMap);

  ($Abbrev) = @_;
  %AbbrevTemplateDataMap = ();

  ($GLType, $ChainsAbbrev, $AbbrevModifier) = ParseGLAbbrev($Abbrev);

  $TemplateID = _GetGLTemplateID($ChainsAbbrev);

  @SnAbbrevs = split /\//, $ChainsAbbrev;
  @SnChainAdd = ();

  for $SnAbbrev (@SnAbbrevs) {
    push @SnChainAdd, ($SnAbbrev eq "0:0") ? 0 : 1;
  }

  @{$AbbrevTemplateDataMap{SnChainAdd}} = ();
  push @{$AbbrevTemplateDataMap{SnChainAdd}}, @SnChainAdd;

  @{$AbbrevTemplateDataMap{SnAbbrev}} = ();
  push @{$AbbrevTemplateDataMap{SnAbbrev}}, @SnAbbrevs;

  $AbbrevModifier =~ s/(\[|\])//g;
  $AbbrevTemplateDataMap{AbbrevModifier} = $AbbrevModifier;

  ChainStr::SetupTemplateDataMap('GL', \%AbbrevTemplateDataMap, $GLTemplatesDataMap{$TemplateID});
  return \%AbbrevTemplateDataMap;
}

# Get GL template ID...
sub _GetGLTemplateID {
  my($ChainsAbbrev) = @_;
  my($TemplateID, $SnAbbrev, $SnChainType, $ChainNum, @SnAbbrevs, @SnChainTypes);

  @SnAbbrevs = split /\//, $ChainsAbbrev;
  @SnChainTypes = ();

  $ChainNum = 0;
  for $SnAbbrev (@SnAbbrevs) {
    $ChainNum++;
    if ($SnAbbrev eq "0:0") {
      $SnChainType = "";
    }
    else {
      $SnChainType = (ChainAbbrev::IsAlkylChainAbbrev($SnAbbrev)) ? "Sn${ChainNum}Alkyl" : ((ChainAbbrev::IsAlkenylChainAbbrev($SnAbbrev)) ? "Sn${ChainNum}Alkenyl" : "Sn${ChainNum}");
    }
    push @SnChainTypes, $SnChainType;
  }
  $TemplateID = join "", @SnChainTypes;

  return $TemplateID;
}

# Initialize GL data...
sub _InitializeData {
  _InitializeStrTemplates();
}

# Initialize structure template data for these supported templates:
#
# . MG:
#       Sn1/0.0/0.0 - Sn1 Acyl group - 16:0/0.0/0.0- GL0101
#       Sn1Alkyl/0.0/0.0 - Sn1 Alkyl group - 16:0e/0.0/0.0 - GL0102
#       Sn1Alkenyl/0.0/0.0 - Sn1 Alkenyl group - 16:0p/0.0/0.0 - GL0103
#
# . DG:
#       Sn1/Sn2/0.0 - Sn1 and Sn2 Acyl group - 16:0/18:0/0.0 - GL0201
#       Sn1Alkyl/Sn2/0.0 - Sn1 Alkyl and Sn2 Acyl group - 16:0e/18:0/0.0 - GL0202
#       Sn1Alkenyl/Sn2/0.0 - Sn1 Alkenyl and Sn2 Acyl group - 16:p/18:0/0:0 - GL0204
#
#       Sn1/0:0/Sn3 - Sn1 and Sn3 acyl group - 16:0/0:0/18.0 - GL0201
#       Sn1Alkyl/0.0/Sn3 - Sn1 Alkyl and Sn3 acyl group - 16:0e/0:0/0.0 - GL0202
#       Sn1Alkenyl/0.0/Sn3 - Sn1 Alenkyl and Sn3 acyl group - 16:p/0:0/18:0 - GL0204
#
#       Sn1/Sn2/Sn3 - Sn1, Sn2 and Sn3 acyl group - 16:0/18:0/18:0 - GL0301
#       Sn1Alkyl/Sn2/Sn3 - Sn1Alkyl, Sn2 and Sn3 acyl group - 16:0e/18:0/18:0 - GL0302
#       Sn1Alkenyl/Sn2/Sn3 - Sn1Alenkyl, Sn2 and Sn3 acyl group - 16:0e/18:0/18:0 - GL0304
#
sub _InitializeStrTemplates {
  %GLTemplatesDataMap = ();

  # MG - sn1/0.0/0.0
  my($MGTemplateString)=<<ENDMGTEMPLATE;
MG sn1 structure template
  LipdMAPS02060609152D

 10  9  0  0  0  0  0  0  0  0999 V2000
    2.1432   -0.2612    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.4289    0.1500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7143   -0.2612    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.1500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7145   -0.2612    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4288    0.1500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4288    0.9757    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.1272   -0.9757    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.3014   -0.9757    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1432   -0.2612    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  3  9  1  6      
  3  8  1  1      
  4  5  1  0      
  5  6  1  0      
  6  7  2  0      
  6 10  1  0      
M  END

ENDMGTEMPLATE

  # MG - sn1/0.0/0.0
  my($MGAlkylTemplateString)=<<ENDMGALKYLTEMPLATE;
MG sn1 alkyl structure template
  LipdMAPS02060609152D

  7  6  0  0  0  0  0  0  0  0999 V2000
    1.4299    0.1518    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.7151    0.5633    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0001    0.1518    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7149    0.5633    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4299    0.1518    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.4132   -0.5633    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4133   -0.5633    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  3  7  1  6      
  3  6  1  1      
  4  5  1  0      
M  END

ENDMGALKYLTEMPLATE

  # DG - sn1/sn2/0.0
  my($DGSn1Sn2TemplateString)=<<ENDDGSN1SN2TEMPLATE;
DG sn1sn2 structure template
  LipdMAPS02060609152D

 13 12  0  0  0  0  0  0  0  0999 V2000
    2.1432    0.3581    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.4289    0.7694    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7143    0.3581    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.7694    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7145    0.3581    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4288    0.7694    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4288    1.5951    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.1272   -0.3564    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.3014   -0.3564    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4129   -0.7692    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4129   -1.5951    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1272   -0.3564    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1432    0.3581    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  3  9  1  6      
  3  8  1  1      
  4  5  1  0      
  5  6  1  0      
  6  7  2  0      
  6 13  1  0      
  9 10  1  0      
 10 11  2  0      
 10 12  1  0      
M  END

ENDDGSN1SN2TEMPLATE

  # DG - sn1/sn2/0.0
  my($DGSn1AlkylSn2TemplateString)=<<ENDDGSN1ALKYLSN2TEMPLATE;
DG sn1sn2 alkyl structure template
  LipdMAPS02060609152D

 10  9  0  0  0  0  0  0  0  0999 V2000
    1.6362    0.7714    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.9214    1.1829    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.2064    0.7714    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5083    1.1829    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2232    0.7714    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.6196    0.0564    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2067    0.0564    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9214   -0.3565    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9214   -1.1829    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6362    0.0564    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  3  7  1  6      
  3  6  1  1      
  4  5  1  0      
  7  8  1  0      
  8  9  2  0      
  8 10  1  0      
M  END

ENDDGSN1ALKYLSN2TEMPLATE

  # DG - sn1/sn2/0.0
  my($DGSn1AlkylSn2AlkylTemplateString)=<<ENDDGSN1ALKYLSN2ALKYLTEMPLATE;
DG sn1sn2 alkyl structure template
  LipdMAPS02060609152D

  8  7  0  0  0  0  0  0  0  0999 V2000
    1.4290    0.3580    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.7146    0.7693    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0001    0.3580    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7145    0.7693    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4290    0.3580    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.4129   -0.3566    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4130   -0.3566    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1274   -0.7693    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  3  7  1  6      
  3  6  1  1      
  4  5  1  0      
  7  8  1  0      
M  END

ENDDGSN1ALKYLSN2ALKYLTEMPLATE

  # DG - sn1/0.0/sn3
  my($DGSn1Sn3TemplateString)=<<ENDDGSN1SN3TEMPLATE;
DG sn1sn3 structure template
  LipdMAPS02060609152D

 13 12  0  0  0  0  0  0  0  0999 V2000
    1.1372   -0.5538    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4226   -0.9652    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2917   -0.5538    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0063   -0.9652    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7206   -0.5538    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7206    0.2719    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.8355   -1.6796    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.0097   -1.6796    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4350   -0.9652    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.1372    0.2712    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.7205    0.8546    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.7205    1.6796    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.4350    0.4420    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  2  8  1  6      
  2  7  1  1      
  3  4  1  0      
  4  5  1  0      
  5  6  2  0      
  5  9  1  0      
  1 10  1  0      
 10 11  1  0      
 11 12  1  0      
 11 13  2  0      
M  END

ENDDGSN1SN3TEMPLATE

  # DG - sn1/0.0/sn3
  my($DGSn1AlkylSn3TemplateString)=<<ENDDGSN1ALKYLSN3TEMPLATE;
DG sn1sn3 alkyl structure template
  LipdMAPS02060609152D

 10  9  0  0  0  0  0  0  0  0999 V2000
    0.4229   -0.5538    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2918   -0.9652    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0061   -0.5538    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7207   -0.9652    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.1212   -1.6796    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7047   -1.6796    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.4229    0.2712    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.0062    0.8446    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0062    1.6796    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.7207    0.4420    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  2  6  1  6      
  2  5  1  1      
  3  4  1  0      
  1  7  1  0      
  7  8  1  0      
  8  9  1  0      
  8 10  2  0      
M  END

ENDDGSN1ALKYLSN3TEMPLATE

  # TG - sn1/sn2/sn3
  my($TGTemplateString)=<<ENDTGTEMPLATE;
TG structure template
  LipdMAPS02060609152D

 16 15  0  0  0  0  0  0  0  0999 V2000
    1.1372    0.0656    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4226   -0.3458    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2917    0.0656    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0063   -0.3458    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7206    0.0656    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7206    0.8913    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.8355   -1.0602    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.0097   -1.0602    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7047   -1.4731    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7047   -2.2990    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4190   -1.0602    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4350   -0.3458    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.1372    0.8906    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.7205    1.4850    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.7205    2.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.4350    1.0614    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  2  8  1  6      
  2  7  1  1      
  3  4  1  0      
  4  5  1  0      
  5  6  2  0      
  5 12  1  0      
  8  9  1  0      
  9 10  2  0      
  9 11  1  0      
  1 13  1  0      
 13 14  1  0      
 14 15  1  0      
 14 16  2  0      
M  END

ENDTGTEMPLATE

  # TG - sn1/sn2/sn3
  my($TGSn1AlkylTemplateString)=<<ENDTGSN1ALKYLTEMPLATE;
TG sn1 alkyl structure template
  LipdMAPS02060609152D

 13 12  0  0  0  0  0  0  0  0999 V2000
    0.6292    0.0656    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0854   -0.3458    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7997    0.0656    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5143   -0.3458    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.3275   -1.0602    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4983   -1.0602    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2127   -1.4731    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2127   -2.2990    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9270   -1.0602    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6292    0.8906    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.2125    1.4950    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2125    2.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.9270    1.0614    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  2  6  1  6      
  2  5  1  1      
  3  4  1  0      
  6  7  1  0      
  7  8  2  0      
  7  9  1  0      
  1 10  1  0      
 10 11  1  0      
 11 12  1  0      
 11 13  2  0      
M  END

ENDTGSN1ALKYLTEMPLATE

  my($TGSn1AlkylSn2AlkylTemplateString)=<<ENDTGSN1ALKYLSN2ALKYLTEMPLATE;
TG sn1/sn2 alkyl structure template
  LipdMAPS02060609152D

 11 10  0  0  0  0  0  0  0  0999 V2000
    0.4228   -0.3473    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2917   -0.7586    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0059   -0.3473    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7204   -0.7586    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.1211   -1.4729    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7045   -1.4729    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4188   -1.8857    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4228    0.4776    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.0060    1.0512    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0060    1.8857    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.7204    0.6483    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  2  6  1  6      
  2  5  1  1      
  3  4  1  0      
  6  7  1  0      
  1  8  1  0      
  8  9  1  0      
  9 10  1  0      
  9 11  2  0      
M  END

ENDTGSN1ALKYLSN2ALKYLTEMPLATE

  # Format: ID => AbbrevID|Y1Sn1|Y2Sn1|Y1Sn2|Y2Sn2|Y1Sn3|Y2Sn3|Sn1AtomNum|Sn2AtomNum|Sn3AtomNum|Sn1CarbonCount|Sn2CarbonCount|Sn3CarbonCount|Sn2CAtomNum|Sn2OAtomNum|Sn2HAtomNum|LMCategory|LMMainClass|LMSubClass|TemplateString
  # Note: Template structure for Both Alkyl and Alkenyl substituted GL is the same; Different entries are used to facilitate assigning LM identity
  % GLTemplatesDataMap = (
		       "Sn1" => "MG|-0.2612|0.1500|0|0|0|0|10|0|0|2|0|0|3|9|8|GL|01|01|$MGTemplateString",
		       "Sn1Alkyl" => "MG|0.1518|0.5633|0|0|0|0|5|0|0|0|0|0|3|7|6|GL|01|02|$MGAlkylTemplateString",
		       "Sn1Alkenyl" => "MG|0.1518|0.5633|0|0|0|0|5|0|0|0|0|0|3|7|6|GL|01|03|$MGAlkylTemplateString",

		       "Sn1Sn2" => "DG|0.3581|0.7694|-0.7692|-0.3564|0|0|13|12|0|2|2|0|3|9|8|GL|02|01|$DGSn1Sn2TemplateString",
		       "Sn1AlkylSn2" => "DG|0.7714|1.1829|-0.3565|0.0564|0|0|5|10|0|0|2|0|3|7|6|GL|02|02|$DGSn1AlkylSn2TemplateString",
		       "Sn1AlkylSn2Alkyl" => "DG|0.3580|0.7693|-0.7693|-0.3566|0|0|5|8|0|0|1|0|3|7|6|GL|02|03|$DGSn1AlkylSn2AlkylTemplateString",
		       "Sn1AlkenylSn2" => "DG|0.7714|1.1829|-0.3565|0.0564|0|0|5|10|0|0|2|0|3|7|6|GL|02|04|$DGSn1AlkylSn2TemplateString",

		       "Sn1Sn3" => "DG|-0.9652|-0.5538|0|0|1.6796|2.0911|9|0|12|2|0|2|2|8|7|GL|02|01|$DGSn1Sn3TemplateString",
		       "Sn1AlkylSn3" => "DG|-0.9652|-0.5538|0|0|1.6796|2.0911|4|0|9|0|0|2|2|6|5|GL|02|02|$DGSn1AlkylSn3TemplateString",
		       "Sn1AlkenylSn3" => "DG|-0.9652|-0.5538|0|0|1.6796|2.0911|4|0|9|0|0|2|2|6|5|GL|02|04|$DGSn1AlkylSn3TemplateString",

		       "Sn1Sn2Sn3" => "TG|-0.3458|0.0656|-1.4731|-1.0602|2.2991|2.7115|12|11|15|2|2|2|2|8|7|GL|03|01|$TGTemplateString",
		       "Sn1AlkylSn2Sn3" => "TG|-0.3458|0.0656|-1.4731|-1.0602|2.2987|2.7115|4|9|12|0|2|2|2|6|5|GL|03|02|$TGSn1AlkylTemplateString",
		       "Sn1AlkylSn2AlkylSn3" => "TG|-0.7586|-0.3473|-1.8857|-1.4729|1.8857|2.2973|4|7|10|0|1|2|2|6|5|GL|03|03|$TGSn1AlkylSn2AlkylTemplateString",
		       "Sn1AlkenylSn2Sn3" => "TG|-0.3458|0.0656|-1.4731|-1.0602|2.2987|2.7115|4|9|12|0|2|2|2|6|5|GL|03|04|$TGSn1AlkylTemplateString",
		      );

}

1;

__END__

=head1 NAME

GLStr - Glycerolipids (GL) structure generation methods

=head1 SYNOPSIS

use GLStr;

use GLStr qw(:all);

=head1 DESCRIPTION

GLStr module provides these methods:

    GenerateCmpdOntologyData - Generate ontology data
    GenerateCmpdOntologySDDataLines - Generate ontology data lines for SD file
    GenerateGLChainStrData - Generate chain structure data
    GenerateSDFile - Generate SD file
    IsGLAbbrevSupported - Is it a supported GL abbreviation
    ParseGLAbrev - Parse GL abbreviation
    ProcessGLCmpdAbbrevs - Process GL abbreviation
    SetupGLCmpdAbbrevTemplateDataMap - Setup template structure data map
    ValidateGLAbbrev - Validate GL abbreviation

=head1 METHODS

=over 4

=item B<ExpandGLCmpdAbbrevs>

    $ExpandedAbbrevArrayRef = ExpandGLCmpdAbbrevs($CmpdAbbrev);

Return a reference to an array containing complete GL abbreviations. Wild card
characters in GL abbreviation name are expanded to generate fully qualified
GL abbreviations.

=item B<GenerateCmpdOntologyData>

    $DataHashRef = GenerateCmpdOntologyData($CmpdDataRef);

Return a reference to a hash containing ontology data with hash keys and values
corresponding to property names and values.

=item B<GenerateCmpdOntologySDDataLines>

    $DataLinesArrayRef = GenerateCmpdOntologySDDataLines($CmpDataRef);

Return a reference to an array containing ontology data lines suitable for
generate SD file data block.


=item B<GenerateGLChainStrData>

    ($AtomLinesArrayRef, $BondLinesArrayRef) =
       GenerateGLChainStrData($ChainType, $CmpdDataRef);

Return array references containing atom and bond data lines for SD file. Appropriate atom
and bond data lines are generated using chain type and abbreviation template data.

=item B<GenerateSDFile>

    GenerateSDFile($SDFileName, $CmdAbbrevsRef);

Generate a SD file for compound abbreviations. Structure data for specified abbreviation
is generated sequentially and written to SD file.

=item B<IsGLAbbrevSupported>

    $Status = IsGLAbbrevSupported($Abbrev, $PrintWarning);

Return 1 or 0 based on whether GL abbreviated is supported. For unsupported GL abbreviations,
a warning is printed unless PrintWarning flag is set.

=item B<ParseGLAbrev>

    ($GLType, $ChainsAbbrev, $AbbrevModifier) = ParseGLAbrev($Abbrev);

Parse GL abbreviation and return these values: GLType, ChainsAbbrev,
AbbrevModifier.

=item B<ProcessGLCmpdAbbrevs>

    ProcessGLCmpdAbbrevs($CmpdAbbrevsRef, $AllowArbitraryChainAbbrev,
                         $WriteSDFile, $SDFileName);

Process specified GL abbreviations to generate structures and write them out either
a SD file or simply report number of valid abbreviations.

=item B<SetupGLCmpdAbbrevTemplateDataMap>

    $AbbrevTemplateDataMapRef =
       SetupGLCmpdAbbrevTemplateDataMap($Abbrev);

Return a reference to a hash containing template data for compound abbreviation. The
template data is used to generate SD file for compound abbreviation.

=item B<ValidateGLAbbrev>

    $Status = ValidateGLAbbrev($Abbrev);

Return 1 or 0 based on whether a GL abbreviation is valid.

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
