package FAStr;
#
# File: FAStr.pm
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
use ChainStr;

our ($VERSION, @ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);

$VERSION = '2.00';
@ISA = qw(Exporter);
@EXPORT = qw();
@EXPORT_OK = qw(GenerateCmpdOntologyData GenerateCmpdOntologySDDataLines GenerateFAChainStrData GenerateSDFile IsFAAbbrevSupported ParseFAAbrev ProcessFACmpdAbbrevs SetupFACmpdAbbrevTemplateDataMap ValidateFAAbbrev);
%EXPORT_TAGS = (
		all  => [@EXPORT, @EXPORT_OK],
		);

# Initialize FA data...
my(%FATemplatesDataMap);
_InitializeData();

#
# Process specified compound abbreviations containing any wild cards
# and count the number of valid abbreviations along with optional generation
# of a SD file.
#
sub ProcessFACmpdAbbrevs {
  my($CmpdAbbrevsRef, $WriteSDFile, $SDFileName) = @_;

  $WriteSDFile = defined $WriteSDFile ? $WriteSDFile : 0;
  if ($WriteSDFile &&  IsEmpty($SDFileName)) {
    warn "Warning: FAStr::ProcessFACmpdAbbrevs: No SD file name specified. Ingoring structure generation.\n";
    return;
  }

  if ($WriteSDFile) {
    print "Generating SD file $SDFileName...\n";

    open SDFILE, ">$SDFileName" or die "Error: Couldn't open $SDFileName: $! \n";

    _ProcessFACmpdAbbrevs($CmpdAbbrevsRef, $WriteSDFile, \*SDFILE);

    close SDFILE;
  }
  else {
    _ProcessFACmpdAbbrevs($CmpdAbbrevsRef, $WriteSDFile);
  }
}

# Process specified compound abbreviations containing any wild cards
# and count the number of valid abbreviations along with optional generation
# of a SD file.
#
sub _ProcessFACmpdAbbrevs {
  my($CmpdAbbrevsRef, $WriteSDFile, $SDFileRef) = @_;
  my($Abbrev, $AbbrevCount, $ChainAbbrev, $ExpandedChainAbbrev, $ExpandedAbbrevRef, @ExpandedCmpdAbbrevs);

  $AbbrevCount = 0;

 ABBREV: for $Abbrev (@{$CmpdAbbrevsRef}) {
    if (!ValidateFAAbbrev($Abbrev)) {
      next ABBREV;
    }

    ($ChainAbbrev) = ParseFAAbbrev($Abbrev);

    #
    # Allow any reasonable combination of chain lengths and bonds...
    #
    if (!IsFAAbbrevSupported($ChainAbbrev)) {
      warn "Warning: Ignoring compound abbreviation $Abbrev\n";
      next ABBREV;
    }

    if (!(ChainAbbrev::IsWildCardInChainAbbrev($ChainAbbrev))) {
      $AbbrevCount++;
      if ($WriteSDFile) {
	_GenerateAndWriteCmdDataString($SDFileRef, $Abbrev);
      }
      next ABBREV;
    }

    # Expand the chain abbreviation...
    $ExpandedAbbrevRef = ChainAbbrev::ExpandChainAbbrev($ChainAbbrev);
    EXPANDEDABBREV: for $ExpandedChainAbbrev (@{$ExpandedAbbrevRef}) {
      if (ChainAbbrev::IsAlkylChainAbbrev($ExpandedChainAbbrev) || ChainAbbrev::IsAlkenylChainAbbrev($ExpandedChainAbbrev)) {
	# Skip these...
	next EXPANDEDABBREV;
      }
      $AbbrevCount++;
      if ($WriteSDFile) {
	_GenerateAndWriteCmdDataString($SDFileRef, $ExpandedChainAbbrev);
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

# Generate appropriate compound data string containing structure data and ontology...
sub _GenerateCmdDataString {
  my($Abbrev) = @_;
  my($CmpdDataString, $AtomBlockLinesRef, $BondBlockLinesRef, $CmpdAbbrevTemplateDataMapRef, $OntologyDataLinesRef, @CmpdDataLines);

  $CmpdDataString = '';
  @CmpdDataLines = ();

  # Setup template data for a specific compound abbreviation...
  $CmpdAbbrevTemplateDataMapRef = SetupFACmpdAbbrevTemplateDataMap($Abbrev);

  # Generate structure data for for the chain...
  ($AtomBlockLinesRef, $BondBlockLinesRef) = GenerateFAChainStrData('Sn1', $CmpdAbbrevTemplateDataMapRef);

  # Generate data block lines including various desriptors...
  ($OntologyDataLinesRef) = GenerateCmpdOntologySDDataLines($CmpdAbbrevTemplateDataMapRef);

  # Write out first four SD file lines: Name, MiscInfo, Comments, Count
  $CmpdDataString .= "$Abbrev\n";
  $CmpdDataString .= LMAPSStr::GenerateCmpdMiscInfoLine(). "\n";
  $CmpdDataString .= LMAPSStr::GenerateCmpdCommentsLine() . "\n";
  $CmpdDataString .= ChainStr::GenerateCmpdCountsLine($CmpdAbbrevTemplateDataMapRef) . "\n";

  # Atom lines for template and chains...
  $CmpdDataString .= ChainStr::GenerateAtomBlockLines($CmpdAbbrevTemplateDataMapRef, $AtomBlockLinesRef) . "\n";

  # Bond lines for template and chains...
  $CmpdDataString .= ChainStr::GenerateBondBlockLines($CmpdAbbrevTemplateDataMapRef, $BondBlockLinesRef) . "\n";

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
  my($ChainIndex, $ChainAbbrev, $ChainLength, $DoubleBondCount, $TripleBondCount);

  # Sn1 chain length and double bond count...
  $ChainIndex = 0;
  $ChainAbbrev = $CmpdAbbrevTemplateDataMapRef->{SnChainAdd}[$ChainIndex] ? $CmpdAbbrevTemplateDataMapRef->{SnAbbrev}[$ChainIndex] : '0:0';

  ($ChainLength, $DoubleBondCount, $TripleBondCount) = ChainAbbrev::GetChainLengthAndMultipleBondCount($ChainAbbrev);

  $OntologyDataMapRef->{'Chain Length'} = $ChainLength;
  $OntologyDataMapRef->{'Double Bonds'} = $DoubleBondCount;
  $OntologyDataMapRef->{'Triple Bonds'} = $TripleBondCount;

}

# Abbreviation and systematic name...
#
#   Setup:
#  . Chain name prefix using chain length.
#  . Sorted double bond position and geomety
#  . Sorted substituents position and name along with handling of multiple substituents at chain carbon position
#  . Based on substitutent at position 1, figure out name ending: NH2 - amine; CHO - aldehyde
#  . Based on number of double bonds modify chain name prefix and combine it with name ending
#
#  And generate a complete name.
#
#  Examples:
#
#  tetradecanoic acid
#  5Z,8Z,11Z,14Z-eicosatetraenoic acid
#  4-hydroxy-heptanoic acid
#  9,12-dimethoxy-13-hydroxy-10-octadecenoic acid
#
#  10E,12Z-hexadecadien-1-ol
#  2R,9R-dihydroxy-3S,4S,7S,8S-diepoxy-5E,10-undecadien-1-ol
#  6-hydroxy-2,4-hexadienal
#
#  octadecanamide
#
# 6-oxo-9S,11R,15S-trihydroxy-13E-prostenoic acid - Abbrev: 20:1(13E)(6Ke,9OH[S],11OH[R],15OH[S]){8a,12b}
# 9S,11R,15S-trihydroxy-5Z,13E-prostadienoic acid - Abbrev: 20:2(5Z,13E)(9OH[S],11OH[R],15OH[S]){8a,12b}
#
sub _SetupAbbrevAndSystematicName {
  my($CmpdAbbrevTemplateDataMapRef, $OntologyDataMapRef) = @_;
  my($ChainIndex, $ChainAbbrev, $ChainName,  $ChainNamePrefix, $CountPrefix, $ChainLenToNamePrefixMap, $CountToNamePrefixMap, $ChainLength, $DoubleBondCount, $DoubleBondGeometry, $Substituents, $Rings, $Substituent, $SubstituentsName, $SystematicName, $DoubleBondNames, $ProstaglandinAbbrev, $TemplateID);

  # Parse chain abbreviation...
  $ChainIndex = 0;
  $ChainAbbrev = $CmpdAbbrevTemplateDataMapRef->{SnChainAdd}[$ChainIndex] ? $CmpdAbbrevTemplateDataMapRef->{SnAbbrev}[$ChainIndex] : '0:0';
  ($ChainLength, $DoubleBondCount, $DoubleBondGeometry, $Substituents, $Rings) = ChainAbbrev::ParseChainAbbrev($ChainAbbrev);

  # Initialize systematic name...
  $SystematicName = '';
  $OntologyDataMapRef->{Abbrev} = $ChainAbbrev;
  $OntologyDataMapRef->{'Systematic Name'} = $SystematicName;

  # Get chain length name and count prefix...
  $ChainLenToNamePrefixMap = ChainAbbrev::GetChainLenToNamePrefixMap();
  $CountToNamePrefixMap = ChainAbbrev::GetCountToNamePrefixMap();
  if (!exists $ChainLenToNamePrefixMap->{$ChainLength}) {
    return;
  }
  $ChainNamePrefix = $ChainLenToNamePrefixMap->{$ChainLength};

  # Is this abbreviation corresponds to prostaglandin subclass?
  $ProstaglandinAbbrev = 0;
  if ($Rings) {
    $ProstaglandinAbbrev = _IsProstaglandinAbbrev($CmpdAbbrevTemplateDataMapRef, $ChainAbbrev);
    if ($ProstaglandinAbbrev) {
      $ChainNamePrefix = 'prost';
    }
  }

  # Setup double bond names by sorting double bond specification by position.
  #
  $DoubleBondNames = '';
  if ($DoubleBondCount) {
    my($DoubleBondPos, $DoubleBondGeometry);
    for $DoubleBondPos (sort {$a <=> $b} keys %{$CmpdAbbrevTemplateDataMapRef->{SnDbleBondPosToGeometry}[$ChainIndex]}) {
      $DoubleBondGeometry = $CmpdAbbrevTemplateDataMapRef->{SnDbleBondPosToGeometry}[$ChainIndex]{$DoubleBondPos};
      if ($DoubleBondNames) {
	$DoubleBondNames .= ",${DoubleBondPos}${DoubleBondGeometry}";
      }
      else {
	$DoubleBondNames = "${DoubleBondPos}${DoubleBondGeometry}";
      }
    }
  }

  # Setup substituent name..
  $SubstituentsName = '';
  if (!IsEmpty($Substituents)) {
    $SubstituentsName = ChainAbbrev::SetupChainSubstituentsName($CmpdAbbrevTemplateDataMapRef, $ChainIndex);
  }

  $ChainName = '';
  # Modify chain name based on number of double bonds...
  if ($DoubleBondCount) {
    $CountPrefix = '';
    if ($DoubleBondCount > 1) {
      $CountPrefix = (exists $CountToNamePrefixMap->{$DoubleBondCount}) ? $CountToNamePrefixMap->{$DoubleBondCount} : '';
      $ChainName = "${ChainNamePrefix}a${CountPrefix}en";
    }
    else {
      $ChainName = "${ChainNamePrefix}en";
    }
  }
  else {
    $ChainName = "${ChainNamePrefix}an";
  }

  # Setup chain name...
  $TemplateID = $CmpdAbbrevTemplateDataMapRef->{TemplateID};
  if ($TemplateID =~ /^OH$/i) {
    if ($SubstituentsName || $DoubleBondNames) {
      $ChainName .= "-1-ol";
    }
    else {
      $ChainName = "1-${ChainName}ol";
    }
  }
  elsif ($TemplateID =~ /^NH2$/i) {
    # Primary amides
    $ChainName .= "amide";
  }
  elsif ($TemplateID =~ /^CHO$/i) {
    # Fatty aldehydes
    $ChainName .= "al";
  }
  else {
    # It's an acid...
    $ChainName .= "oic acid";
  }

  # Setup systematic name...
  $SystematicName = '';
  if ($SubstituentsName) {
    $SystematicName .= "${SubstituentsName}-";
  }
  if ($DoubleBondNames) {
    $SystematicName .= "${DoubleBondNames}-";
  }
  $SystematicName .= "${ChainName}";

  $OntologyDataMapRef->{Abbrev} = $ChainAbbrev;
  $OntologyDataMapRef->{'Systematic Name'} = $SystematicName;
}

# LM classification info...
sub _SetupLMClassificationInfo {
  my($CmpdAbbrevTemplateDataMapRef, $OntologyDataMapRef) = @_;
  my($LMCategory, $LMMainClass, $LMSubClass, $ChainIndex, $ChainAbbrev, $TemplateID,);

  if (IsEmpty($CmpdAbbrevTemplateDataMapRef->{LMCategory})) {
    return;
  }

  $LMCategory = $CmpdAbbrevTemplateDataMapRef->{LMCategory};
  $LMMainClass = '';
  $LMSubClass = '';

  # Setup chain abbreviation...
  $ChainIndex = 0;
  $ChainAbbrev = $CmpdAbbrevTemplateDataMapRef->{SnChainAdd}[$ChainIndex] ? $CmpdAbbrevTemplateDataMapRef->{SnAbbrev}[$ChainIndex] : '0:0';
  $TemplateID = $CmpdAbbrevTemplateDataMapRef->{TemplateID};

  ($LMMainClass, $LMSubClass) = _GetLMMainAndSubClasses($CmpdAbbrevTemplateDataMapRef, $TemplateID, $ChainAbbrev);

  $OntologyDataMapRef->{'LM Category'} = $LMCategory;
  $OntologyDataMapRef->{'LM Main Class'} = $LMCategory . $LMMainClass;
  $OntologyDataMapRef->{'LM Sub Class'} = $LMCategory . $LMMainClass . $LMSubClass;
}

# Generate atom and bond block lines for a chain
sub GenerateFAChainStrData {
  my($ChainType, $CmpdAbbrevTemplateDataMapRef) = @_;

  return ChainStr::GenerateChainStrData($ChainType, $CmpdAbbrevTemplateDataMapRef);
}

# Parse abbreviation...
sub ParseFAAbbrev {
  my($Abbrev) = @_;
  my($ChainAbbrev, $AbbrevModifier);

  if ($Abbrev =~ /\]$/) {
    ($ChainAbbrev, $AbbrevModifier) = $Abbrev =~ /^(.+?)\[(.+?)\]$/;
    $AbbrevModifier = '[' . $AbbrevModifier . ']';
  }
  else {
    $ChainAbbrev = $Abbrev;
    $AbbrevModifier = '';
  }
  return ($ChainAbbrev, $AbbrevModifier);
}

# Does template exist to handle this abbreviation?
#
sub IsFAAbbrevSupported {
  my($ChainAbbrev) = @_;
  my($ChainLength, $DoubleBondCount, $DoubleBondGeometry, $Substituents, $Rings, $AllowRings);

  $AllowRings = 1;

  ($ChainLength, $DoubleBondCount, $DoubleBondGeometry, $Substituents, $Rings) = ChainAbbrev::ParseChainAbbrev($ChainAbbrev);

  if (!_IsFAChainLengthOkay($ChainAbbrev, $ChainLength)) {
    return 0;
  }

  if ($DoubleBondCount && $DoubleBondGeometry) {
    if (!ChainAbbrev::IsDoubleBondsAbbrevOkay($ChainAbbrev, $ChainLength, $DoubleBondCount, $DoubleBondGeometry)) {
      return 0;
    }
  }

  if ($Substituents) {
    if (!ChainAbbrev::IsSubstituentsAbbrevOkay($ChainAbbrev, $ChainLength, $DoubleBondCount, $DoubleBondGeometry, $Substituents)) {
      return 0;
    }
  }

  if ($Rings) {
    if (!$AllowRings) {
      warn "Warning: Ignoring chain abbreviation $ChainAbbrev : Rings are not allowed\n";
      return 0;
    }
    if (!ChainAbbrev::IsRingsAbbrevOkay($ChainAbbrev, $ChainLength, $Rings)) {
      return 0;
    }
    # Only five membered ring is allowed to support prostaglandin structure generation...
    my($RingPos1, $RingPos2, $StereoChemistry1, $StereoChemistry2, $RingSize, @RingWords);
    @RingWords = quotewords(',', 0, $Rings);
    ($RingPos1, $StereoChemistry1) = ChainAbbrev::ParseRingAbbrev($RingWords[0]);
    ($RingPos2, $StereoChemistry2) = ChainAbbrev::ParseRingAbbrev($RingWords[1]);
    $RingSize = $RingPos2 - $RingPos1 + 1;
    if ($RingSize != 5) {
      warn "Warning: Ignoring chain abbreviation $ChainAbbrev : Unknown ring size, $RingSize: Supported ring size: 5...\n";
      return 0;
    }
    # Don't allow triple bond with ring specifications...
    if ($DoubleBondCount && $DoubleBondGeometry) {
      if ($DoubleBondGeometry =~ /Y/i) {
	# Triple bond specification...
	warn "Warning: Ignoring chain abbreviation $ChainAbbrev : Triple bond on the chain is not supported with ring specification $Rings..\n";
	return 0;
      }
    }

    # Ring specification is only allowed for chain with COOH at position 1: Prostaglandin supporty...
    if ($Substituents) {
      my($RingPos1SubstituentCount, $RingPos2SubstituentCount, $Substituent, $SubstituentPos, $SubstituentAbbrev, $SubstituentStereoChemistry, @SubstituentsList);
      ($RingPos1SubstituentCount, $RingPos2SubstituentCount) = (0, 0);
      (@SubstituentsList) = split /\,/, $Substituents;
      for $Substituent (@SubstituentsList) {
	($SubstituentPos, $SubstituentAbbrev, $SubstituentStereoChemistry) = ChainAbbrev::ParseSubstituentAbbrev($Substituent);
	if ($SubstituentPos <= 1 && $SubstituentAbbrev !~ /^(COOH)$/i) {
	  warn "Warning: Ignoring chain abbreviation $ChainAbbrev : Substituent specification, $Substituent, in not valid: Ring specification, $Rings, is only allowed with COOH substituent.\n";
	  return 0;
	}
	if ($SubstituentPos == $RingPos1) {
	  $RingPos1SubstituentCount++;
	}
	elsif ($SubstituentPos == $RingPos2) {
	  $RingPos2SubstituentCount++;
	}
      }
      if ($RingPos1SubstituentCount > 1) {
	warn "Warning: Ignoring chain abbreviation $ChainAbbrev : Number of substituents, $RingPos1SubstituentCount, at ring position, $RingPos1, must be equal to 1...\n";
	return 0;
      }
      if ($RingPos2SubstituentCount > 1) {
	warn "Warning: Ignoring chain abbreviation $ChainAbbrev : Number of substituents, $RingPos2SubstituentCount, at ring position, $RingPos2, must be equal to 1...\n";
	return 0;
      }
    }
  }
  return 1;
}

# Is it a valid chain length abbreviation...
sub _IsFAChainLengthOkay {
  my($ChainAbbrev, $ChainLength) = @_;

  if ($ChainLength =~ /\*/) {
    return 1;
  }

  if ($ChainLength !~ /[0-9]/) {
    warn "Warning: Ignoring chain abbreviation $ChainAbbrev : Chain length, $ChainLength, must be an integar...\n";
    return 0;
  }
  my($ChainLenToNamePrefixMap);
  $ChainLenToNamePrefixMap = ChainAbbrev::GetChainLenToNamePrefixMap();
  if (!exists $ChainLenToNamePrefixMap->{$ChainLength}) {
    warn "Warning: Ignoring chain abbreviation $ChainAbbrev : Systematic name generation for chain length, $ChainLength, is not supported...\n";
    return 0;
  }

  return 1;
}

#
# Check out the validity of FA abbreviation...
sub ValidateFAAbbrev {
  my($Abbrev) = @_;

  if ($Abbrev =~ /\//) {
    warn "Warning: Ignored FA compound abbreviation $Abbrev due to incorrect format: \"/\" is not allowed. Valid format: ChainLenSpec:DoubleBondCount(DoubleBondsSpec(SubstituentsSpec){RingSpec}.\n";
    return 0;
  }

  # Number of '(' and ')' parenthesis must match...
  my(@LeftParenthesis) = $Abbrev =~ /\(/g;
  my(@RightParenthesis) = $Abbrev =~ /\)/g;
  if ($#LeftParenthesis != $#RightParenthesis) {
    warn "Warning: Ignored FA compound abbreviation $Abbrev due to incorrect format: Number of ( and ) doesn't match. Valid format: ChainLenSpec:DoubleBondCount(DoubleBondsSpec(SubstituentsSpec){RingSpec}.\n";
    return 0;
  }

  @LeftParenthesis = $Abbrev =~ /\[/g;
  @RightParenthesis = $Abbrev =~ /\]/g;
  if ($#LeftParenthesis != $#RightParenthesis) {
    warn "Warning: Ignored FA compound abbreviation $Abbrev due to incorrect format: Number of [ and ] doesn't match. Valid format: ChainLenSpec:DoubleBondCount(DoubleBondsSpec(SubstituentsSpec){RingSpec}.\n";
    return 0;
  }

  my($ChainsAbbrev, $AbbrevModifier) = ParseFAAbbrev($Abbrev);
  if ($AbbrevModifier) {
    warn "Warning: Ignored FA compound abbreviation $Abbrev due to incorrect format: Abbreviation modifier is not allowed.  Valid format: ChainLenSpec:DoubleBondCount(DoubleBondsSpec(SubstituentsSpec){RingSpec}.\n";
    return 0;
  }

  my(@ChainsAbbrevList) = ();
  @ChainsAbbrevList = split /\//, $ChainsAbbrev;
  if (@ChainsAbbrevList != 1) {
    warn "Warning: Ignored FA compound abbreviation $Abbrev due to incorrect format: Only one chain is allowed.\n";
    return 0;
  }

  my($Sn1Abbrev) = $ChainsAbbrev;
  if ($Sn1Abbrev eq "0:0") {
    warn "Warning: Ignoring FA compound abbreviation $Abbrev : abbreviation value of 0:0 is not allowed.\n";
    return 0;
  }

  # Alkyl and Alkenyl chain abbbreviatios not supported...
  if (ChainAbbrev::IsAlkylChainAbbrev($Sn1Abbrev)) {
    warn "Warning: Ignoring FA compound abbreviation $Abbrev : Alkyl chain abbreviation is not allowed.\n";
    return 0;
  }
  if (ChainAbbrev::IsAlkenylChainAbbrev($Sn1Abbrev)) {
    warn "Warning: Ignoring FA compound abbreviation $Abbrev : Alkenyl chain abbreviation is not allowed.\n";
    return 0;
  }

  return 1;
}

# Set up template data for a specific cmpd abbreviation and return a
# reference to the hash...
sub SetupFACmpdAbbrevTemplateDataMap {
  my($Abbrev, $AbbrevModifier, $Sn1Abbrev, $Sn2Abbrev, $Sn3Abbrev, $TemplateID, @SnChainAdd, %AbbrevTemplateDataMap);

  ($Abbrev) = @_;

  ($Sn1Abbrev, $Sn2Abbrev, $Sn3Abbrev) = ('0:0') x 3;
  @SnChainAdd = ();
  %AbbrevTemplateDataMap = ();

  ($Sn1Abbrev) = ParseFAAbbrev($Abbrev);
  $TemplateID = _GetTemplateID($Sn1Abbrev);

  push @SnChainAdd, (1,0,0);

  @{$AbbrevTemplateDataMap{SnChainAdd}} = ();
  push @{$AbbrevTemplateDataMap{SnChainAdd}}, @SnChainAdd;

  @{$AbbrevTemplateDataMap{SnAbbrev}} = ();
  push @{$AbbrevTemplateDataMap{SnAbbrev}}, ($Sn1Abbrev, $Sn2Abbrev, $Sn3Abbrev);

  $AbbrevTemplateDataMap{TemplateID} = $TemplateID;
  $AbbrevTemplateDataMap{AbbrevModifier} = '';

  ChainStr::SetupTemplateDataMap('FA', \%AbbrevTemplateDataMap, $FATemplatesDataMap{$TemplateID});
  return \%AbbrevTemplateDataMap;
}

# Assign LM MainClass and SubClass...
sub _GetLMMainAndSubClasses {
  my($CmpdAbbrevTemplateDataMapRef, $TemplateID, $ChainAbbrev) = @_;
  my($LMMainClass, $LMSubClass) = ('') x 2;

  if ($TemplateID =~ /^OH$/i) {
    # Fatty alcohols
    $LMMainClass = '05'; $LMSubClass = '00';
  }
  elsif ($TemplateID =~ /^CHO$/i) {
    # Fatty aldehydes
    $LMMainClass = '06'; $LMSubClass = '00';
  }
  elsif ($TemplateID =~ /^NH2$/i) {
    # Primary amides
    $LMMainClass = '09'; $LMSubClass = '01';
  }
  elsif (_IsProstaglandinAbbrev($CmpdAbbrevTemplateDataMapRef, $ChainAbbrev)) {
    $LMMainClass = '03'; $LMSubClass = '01';
  }
  else {
    # It's an acid..
    $LMMainClass = '01';
    $LMSubClass = _GetFattyAcidsAndConjugatesSubClass($CmpdAbbrevTemplateDataMapRef, $ChainAbbrev);
  }

  return ($LMMainClass, $LMSubClass);
}

# Get sub class...
sub _GetFattyAcidsAndConjugatesSubClass {
  my($CmpdAbbrevTemplateDataMapRef, $ChainAbbrev) = @_;
  my($LMSubClass, $ChainLength, $DoubleBondCount, $DoubleBondGeometry, $Substituents, $Rings);

  $LMSubClass = '';
  ($ChainLength, $DoubleBondCount, $DoubleBondGeometry, $Substituents, $Rings) = ChainAbbrev::ParseChainAbbrev($ChainAbbrev);

  # Sub class priorities for supported classess:
  # Halogenated fatty acids [FA0109]
  # Thia fatty acids [FA0113]
  # Dicarboxylic acids [FA0117]
  # Oxo fatty acids [FA0106]
  # Hydroperoxy fatty acids [FA0104]
  # Epoxy fatty acids [FA0107]
  # Methoxy fatty acids [FA0108]
  # Hydroxy fatty acids [FA0105]
  # Nitro fatty acids [FA0112]
  # Amino fatty acids [FA0110]
  # Unsaturated fatty acids [FA0103]
  # Cyano fatty acids [FA0111]
  # Methyl branched fatty acids [FA0102]
  # Straight chain fatty acids [FA0101]
  #
  # *Carbocyclic fatty acids [FA0114]
  # *Heterocyclic fatty acids [FA0115]
  # *Mycolic acids [FA0116]
  #
  # Br,Cl,F,SH,COOH,KETONES,OOH,EPOXIDES,OME,OH,NH2,TRIPLEBONDS,DOUBLEBONDS,PR,ET,METHYLENES,METHYLS
  #
  SUBCLASS: {
    if ($Substituents =~ /(F|Cl|Br|I)/i) { $LMSubClass = '09'; last SUBCLASS;}
    if ($Substituents =~ /SH/i) { $LMSubClass = '13'; last SUBCLASS;}
    if ($Substituents =~ /COOH/i) { $LMSubClass = '17'; last SUBCLASS;}
    if ($Substituents =~ /Ke/i) { $LMSubClass = '06'; last SUBCLASS;}
    if ($Substituents =~ /OOH/i) { $LMSubClass = '04'; last SUBCLASS;}
    if ($Substituents =~ /Ep/i) { $LMSubClass = '07'; last SUBCLASS;}
    if ($Substituents =~ /OMe/i) { $LMSubClass = '08'; last SUBCLASS;}
    if ($Substituents =~ /OH/i) { $LMSubClass = '05'; last SUBCLASS;}
    if ($Substituents =~ /NO2/i) { $LMSubClass = '12'; last SUBCLASS;}
    if ($Substituents =~ /NH2/i) { $LMSubClass = '10'; last SUBCLASS;}
    if ($DoubleBondCount) { $LMSubClass = '03'; last SUBCLASS;}
    if ($Substituents =~ /CN/i) { $LMSubClass = '11'; last SUBCLASS;}
    if ($Substituents =~ /(Me|Et|Pr|My)/i) { $LMSubClass = '02'; last SUBCLASS;}
    $LMSubClass = '01';
  }

  return $LMSubClass;
}

# Get template ID...
sub _GetTemplateID {
  my($ChainAbbrev, $TemplateID, $ChainLength, $DoubleBondCount, $DoubleBondGeometry, $Substituents, $Rings, $SubstituentPos, $SubstituentAbbrev, $SubstituentStereoChemistry, $Substituent, @SubstituentsList);

  ($ChainAbbrev) = @_;

  ($ChainLength, $DoubleBondCount, $DoubleBondGeometry, $Substituents, $Rings) = ChainAbbrev::ParseChainAbbrev($ChainAbbrev);

  $TemplateID = '';
  (@SubstituentsList) = split /\,/, $Substituents;
  SUBSTITUENT: for $Substituent (@SubstituentsList) {
    ($SubstituentPos, $SubstituentAbbrev, $SubstituentStereoChemistry) = ChainAbbrev::ParseSubstituentAbbrev($Substituent);
    if ($SubstituentPos == 1) {
      if ($SubstituentAbbrev =~ /^OH$/i) {
	$TemplateID = 'OH';
	last SUBSTITUENT;
      }
      elsif ($SubstituentAbbrev =~ /^NH2$/i) {
	$TemplateID = 'NH2';
	last SUBSTITUENT;
      }
      elsif ($SubstituentAbbrev =~ /^CHO$/i) {
	$TemplateID = 'CHO';
	last SUBSTITUENT;
      }
    }
  }
  if ($Rings) {
    my($RingPosRequiresReversal, $BondGeometryRequiresReversal);

    $RingPosRequiresReversal = 0;
    $BondGeometryRequiresReversal = 0;

    # Is the ring starting position odd numbered and required template reversal?
    my($RingPos1, $RingPos2, $StereoChemistry1, $StereoChemistry2, @RingWords);
    @RingWords = quotewords(',', 0, $Rings);
    ($RingPos1, $StereoChemistry1) = ChainAbbrev::ParseRingAbbrev($RingWords[0]);
    ($RingPos2, $StereoChemistry2) = ChainAbbrev::ParseRingAbbrev($RingWords[1]);

    $RingPosRequiresReversal = ($RingPos1 % 2) ? 1 : 0;

    if ($DoubleBondGeometry) {
      # Does bond geometry contain cis double bonds at postions requiring reverse template?
      my($BondGeometry, @DoubleBondGeometryList);
      @DoubleBondGeometryList = split /\,/, $DoubleBondGeometry;
      GEOMERTY: for $BondGeometry (@DoubleBondGeometryList) {
	if ($BondGeometry =~ /^(3Z|4Z|5Z|6Z)$/i) {
	  $BondGeometryRequiresReversal = 1;
	  last GEOMERTY;
	}
      }
    }
    if ($RingPosRequiresReversal && $BondGeometryRequiresReversal) {
      # No reversal...
      $TemplateID = 'COOH';
    }
    elsif ($RingPosRequiresReversal || $BondGeometryRequiresReversal) {
      $TemplateID = 'COOHReversed';
    }

  }
  if (!$TemplateID) {
    $TemplateID = 'COOH';
  }
  return $TemplateID;
}

# Is it a prostaglandin abbreviation?
sub _IsProstaglandinAbbrev {
  my($CmpdAbbrevTemplateDataMapRef, $ChainAbbrev) = @_;
  my($ChainIndex, $RingSize, $ChainLength, $StartRingPos, $EndRingPos, $Status);

  $ChainIndex = 0;
  if (!$CmpdAbbrevTemplateDataMapRef->{SnRing}[$ChainIndex]) {
    return 0;
  }
  $ChainLength = $CmpdAbbrevTemplateDataMapRef->{SnChainLength}[$ChainIndex];
  $RingSize = $CmpdAbbrevTemplateDataMapRef->{SnRingPosInfo}[$ChainIndex]{RingSize};
  $StartRingPos = $CmpdAbbrevTemplateDataMapRef->{SnRingPosInfo}[$ChainIndex]{StartRingPos};
  $EndRingPos = $CmpdAbbrevTemplateDataMapRef->{SnRingPosInfo}[$ChainIndex]{EndRingPos};

  $Status = ($ChainLength == 20 && $RingSize == 5 && $StartRingPos == 8 && $EndRingPos == 12) ? 1 : 0;

  return $Status;
}


# Initialize FA data...
sub _InitializeData {
  _InitializeStrTemplates();
}

# Initialize structure template data for these supported templates:
#
sub _InitializeStrTemplates {
  %FATemplatesDataMap = ();

  # COOH
  my($COOHTemplateString)=<<ENDTEMPLATE;
COOH structure template
  LipdMAPS02280709152D

  3  2  0  0  0  0  0  0  0  0999 V2000
   -0.3572   -0.2063    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3572   -0.6187    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3572    0.6187    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  2  1  1  0      
  1  3  2  0      
M  END

ENDTEMPLATE

  # COOHReversed
  my($COOHReversedTemplateString)=<<ENDTEMPLATE;
COOHReversed structure template
  LipdMAPS02280709152D

  3  2  0  0  0  0  0  0  0  0999 V2000
   -0.3572    0.1031    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3572    0.5156    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3572   -0.5156    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  2  1  1  0      
  1  3  2  0      
M  END

ENDTEMPLATE

  # OH
  my($OHTemplateString)=<<ENDTEMPLATE;
OH structure template
  LipdMAPS02280709152D

  2  1  0  0  0  0  0  0  0  0999 V2000
   -0.3572   -0.2063    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3572   -0.6187    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  2  1  1  0      
M  END

ENDTEMPLATE

  # Amindes: NH2
  my($NH2TemplateString)=<<ENDTEMPLATE;
NH2 structure template
  LipdMAPS02280709152D

  3  2  0  0  0  0  0  0  0  0999 V2000
   -0.3572   -0.2063    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3572   -0.6187    0.0000 NH2 0  0  0  0  0  0  0  0  0  0  0  0
   -0.3572    0.6187    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  2  1  1  0      
  1  3  2  0      
M  END

ENDTEMPLATE

  # CHO
  my($CHOTemplateString)=<<ENDTEMPLATE;
CHO structure template
  LipdMAPS02280709152D

  3  2  0  0  0  0  0  0  0  0999 V2000
   -0.3572   -0.2063    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3572   -0.6187    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3572    0.6187    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  2  1  1  0      
  1  3  2  0      
M  END

ENDTEMPLATE

  my($COOMeTemplateString)=<<ENDTEMPLATE;
COOMe structure template
  LipdMAPS02280709152D

  4  3  0  0  0  0  0  0  0  0999 V2000
   -0.7144   -0.1240    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0000   -0.5366    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7144    0.5366    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.7144   -0.1241    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  1  3  2  0      
  2  4  1  0      
M  END

ENDTEMPLATE

  # Format: ID => AbbrevID|Y1Sn1|Y2Sn1|Y1Sn2|Y2Sn2|Y1Sn3|Y2Sn3|Sn1AtomNum|Sn2AtomNum|Sn3AtomNum|Sn1CarbonCount|Sn2CarbonCount|Sn3CarbonCount|Sn2CAtomNum|Sn2OAtomNum|Sn2HAtomNum|LMCategory|LMMainClass|LMSubClass|TemplateString
  % FATemplatesDataMap = (
		       "COOH" => "FA|-0.6187|-0.2063|0|0|0|0|1|0|0|1|0|0|0|0|0|FA|||$COOHTemplateString",
		       "COOHReversed" => "FA|0.1031|0.5156|0|0|0|0|1|0|0|1|0|0|0|0|0|FA|||$COOHReversedTemplateString",

		       "OH" => "FA|-0.6187|-0.2063|0|0|0|0|1|0|0|1|0|0|0|0|0|FA|||$OHTemplateString",
		       "NH2" => "FA|-0.6187|-0.2063|0|0|0|0|1|0|0|1|0|0|0|0|0|FA|||$NH2TemplateString",
		       "CHO" => "FA|-0.6187|-0.2063|0|0|0|0|1|0|0|1|0|0|0|0|0|FA|||$CHOTemplateString",

		       "COOMe" => "FA|-0.5366|-0.1240|0|0|0|0|1|0|0|1|0|0|0|0|0|FA|||$COOMeTemplateString",
		      );

}

1;

__END__

=head1 NAME

FAStr - Fatty Acyls (FA) structure generation methods

=head1 SYNOPSIS

use FAStr;

use FAStr qw(:all);

=head1 DESCRIPTION

FAStr module provides these methods:

    GenerateCmpdOntologyData - Generate ontology data
    GenerateCmpdOntologySDDataLines - Generate ontology data lines for SD file
    GenerateFAChainStrData - Generate chain structure data
    GenerateSDFile - Generate SD file
    IsFAAbbrevSupported - Is it a supported FA abbreviation
    ParseFAAbrev - Parse FA abbreviation
    SetupFACmpdAbbrevTemplateDataMap - Setup template structure data map
    ValidateFAAbbrev - Validate FA abbreviation

=head1 METHODS

=over 4

=item B<GenerateCmpdOntologyData>

    $DataHashRef = GenerateCmpdOntologyData($CmpdDataRef);

Return a reference to a hash containing ontology data with hash keys and values
corresponding to property names and values.

=item B<GenerateCmpdOntologySDDataLines>

    $DataLinesArrayRef = GenerateCmpdOntologySDDataLines($CmpDataRef);

Return a reference to an array containing ontology data lines suitable for
generate SD file data block.

=item B<GenerateFAChainStrData>

    ($AtomLinesArrayRef, $BondLinesArrayRef) =
       GenerateFAChainStrData($ChainType, $CmpdDataRef);

Return array references containing atom and bond data lines for SD file. Appropriate atom
and bond data lines are generated using chain type and abbreviation template data.

=item B<GenerateSDFile>

    GenerateSDFile($SDFileName, $CmdAbbrevsRef);

Generate a SD file for compound abbreviations. Structure data for specified abbreviation
is generated sequentially and written to SD file.

=item B<IsFAAbbrevSupported>

    $Status = IsFAAbbrevSupported($Abbrev, $PrintWarning);

Return 1 or 0 based on whether FA abbreviated is supported. For unsupported FA abbreviations,
a warning is printed unless PrintWarning flag is set.

=item B<ParseFAAbrev>

    ($ChainAbbrev, $AbbrevModifier) = ParseFAAbrev($Abbrev);

Parse FA abbreviation and return these values: ChainsAbbrev and AbbrevModifier.

=item B<ProcessFACmpdAbbrevs>

    ProcessFACmpdAbbrevs($CmpdAbbrevsRef, $WriteSDFile, $SDFileName);

Process specified FA abbreviations to generate structures and write them out either
a SD file or simply report number of valid abbreviations.

=item B<SetupFACmpdAbbrevTemplateDataMap>

    $AbbrevTemplateDataMapRef =
       SetupFACmpdAbbrevTemplateDataMap($Abbrev);

Return a reference to a hash containing template data for compound abbreviation. The
template data is used to generate SD file for compound abbreviation.

=item B<ValidateFAAbbrev>

    $Status = ValidateFAAbbrev($Abbrev);

Return 1 or 0 based on whether a FA abbreviation is valid.

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
