package TextUtil;
#
# $RCSfile: TextUtil.pm,v $
# $Date: 2007/09/05 15:55:59 $
# $Revision: 1.1 $
#
# Author: Manish Sud <msud@san.rr.com>
#
# Copyright (C) 2004-2005 Manish Sud. All rights reserved.
#
# This library is free software; you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the Free
# Software Foundation; either version 2.1 of the License, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
# or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General
# Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program; if not, write to the Free Software Foundation Inc.,
# 59 Temple Place, Suite 330, Boston, MA  02111-1307, USA.
#
use 5.006;
use strict;
use Exporter;

our ($VERSION, @ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);

$VERSION = '1.00';
@ISA = qw(Exporter);
@EXPORT = qw(AddNumberSuffix ContainsWhiteSpaces GetTextLine IsInteger IsPositiveInteger IsFloat IsEmpty IsNotEmpty IsNumerical JoinWords QuoteAWord RemoveLeadingWhiteSpaces RemoveTrailingWhiteSpaces RemoveLeadingAndTrailingWhiteSpaces WrapText);
@EXPORT_OK = qw();
%EXPORT_TAGS = (
		all  => [@EXPORT, @EXPORT_OK],
		);

# Add number suffix...
sub AddNumberSuffix {
  my($Value) = @_;
  my($ValueWithSuffix, $Suffix);

  $ValueWithSuffix = $Value;
  if (!IsPositiveInteger($Value)) {
    return $ValueWithSuffix;
  }
  $Suffix = "th";
  if ($Value < 10 || $Value > 20) {
    my $Remainder = $Value % 10;
    $Suffix = ($Remainder == 1) ? "st" : (($Remainder == 2) ? "nd" : (($Remainder == 3) ? "rd" : "th"));
  }
  $ValueWithSuffix = "${ValueWithSuffix}${Suffix}";
  return $ValueWithSuffix;
}

# Check out the string: Doen it contain any white space characters?
sub ContainsWhiteSpaces {
  my($TheString) = @_;
  my($Status) = 0;

  if (defined($TheString) && length($TheString)) {
    $Status = ($TheString =~ /[ \t\r\n\f]/ ) ? 1 : 0;
  }
  return $Status;
}

# Read the line, change to UNIX new line char, and chop off new line char as well...
sub GetTextLine {
  my($TextFile) = @_;
  my($Line) = "";

  # Get the next non empty line...
  LINE: while (<$TextFile>) {
      # Change Windows and Mac new line char to UNIX...
      s/(\r\n)|(\r)/\n/g;
      chomp;
      $Line = $_;
      if ($Line) {
	last LINE;
      }
      else {
	next LINE;
      }
    }
  return $Line;
}

# Check out the string: Is it an integer?
sub IsInteger {
  my($TheString) = @_;
  my($Status) = 0;

  if (defined($TheString) && length($TheString)) {
    $TheString = RemoveLeadingAndTrailingWhiteSpaces($TheString);
    $TheString =~ s/^[+-]//;
    $Status = ($TheString =~ /[^0-9]/) ? 0 : 1;
  }
  return $Status;
}

# Check out the string: Is it an integer with value > 0?
sub IsPositiveInteger {
  my($TheString) = @_;
  my($Status) = 0;

  $Status = IsInteger($TheString) ? ($TheString > 0 ? 1 : 0) : 0;

  return $Status;
}


# Check out the string: Is it a float?
sub IsFloat {
  my($TheString) = @_;
  my($Status) = 0;

  if (defined($TheString) && length($TheString)) {
    $TheString = RemoveLeadingAndTrailingWhiteSpaces($TheString);
    $TheString =~ s/^[+-]//;
    $Status = ($TheString =~ /[^0-9.eE]/) ? 0 : 1;
  }
  return $Status;
}

# Check out the string: Is it defined and has a non zero length?
sub IsEmpty {
  my($TheString) = @_;
  my($Status) = 0;

  $Status = (defined($TheString) && length($TheString)) ? 0 : 1;

  return $Status;
}

# Check out the string: Is it defined and has a non zero length?
sub IsNotEmpty {
  my($TheString) = @_;
  my($Status) = 0;

  $Status = (defined($TheString) && length($TheString)) ? 1 : 0;

  return $Status;
}

# Check out the string: Does it only contain numerical data?
sub IsNumerical {
  my($TheString) = @_;
  my($Status) = 0;

  if (defined($TheString) && length($TheString)) {
    $TheString = RemoveLeadingAndTrailingWhiteSpaces($TheString);
    $TheString =~ s/^[+-]//;
    $Status = ($TheString =~ /[^0-9.eE]/) ? 0 : 1;
  }
  return $Status;
}

# Join different words using delimiter and quote parameters. And return as
# a string value.
sub JoinWords {
  my($Words, $Delim, $Quote) = @_;
  my($JoinedWords) = "";

  if (@$Words) {
    my(@NewWords) = map { (defined($_) && length($_)) ? $_ : "" } @$Words;
    if ($Quote) {
      @NewWords = map { "\"$_\"" } @NewWords;
    }
    $JoinedWords = join $Delim, @NewWords;
  }
  return ($JoinedWords);
}

# Based on quote parameter, figure out what to do
sub QuoteAWord {
  my($Word, $Quote) = @_;
  my($QuotedWord);

  $QuotedWord = "";
  if ($Word) {
    $QuotedWord = $Word;
  }
  if ($Quote) {
    $QuotedWord = "\"$QuotedWord\"";
  }
  return ($QuotedWord);
}

# Remove leading white space characters from the string...
sub RemoveLeadingWhiteSpaces {
  my($InString) = @_;
  my($OutString, $TrailingString, $LeadingWhiteSpace);

  $OutString = $InString;
  if (length($InString) && ContainsWhiteSpaces($InString)) {
    $OutString =~ s/^([ \t\r\n\f]*)(.*?)$/$2/;
  }
  return $OutString;
}

# Remove Trailing white space characters from the string...
sub RemoveTrailingWhiteSpaces {
  my($InString) = @_;
  my($OutString, $LeadingString, $TrailingWhiteSpace);

  $OutString = $InString;
  if (length($InString) && ContainsWhiteSpaces($InString)) {
    $OutString =~ s/^(.*?)([ \t\r\n\f]*)$/$1/;
  }
  return $OutString;
}

# Remove both leading and trailing white space characters from the string...
sub RemoveLeadingAndTrailingWhiteSpaces {
  my($InString) = @_;
  my($OutString);

  $OutString = $InString;
  if (length($InString) &&ContainsWhiteSpaces($InString)) {
    $OutString =~ s/^([ \t\r\n\f]*)(.*?)([ \t\r\n\f]*)$/$2/;
  }
  return $OutString;
}

# Wrap text string...
sub WrapText {
  my($InString, $WrapLength, $WrapDelimiter);
  my($OutString);

  $WrapLength = 40;
  $WrapDelimiter = "\n";
  if (@_ == 3) {
    ($InString, $WrapLength, $WrapDelimiter) = @_;
  }
  elsif (@_ == 2) {
    ($InString, $WrapLength) = @_;
  }
  else {
    ($InString, $WrapLength) = @_;
  }
  $OutString = $InString;
  if ($InString && (length($InString) > $WrapLength)) {
    $OutString = "";
    my($Index, $Length, $FirstPiece, $StringPiece);
    $Index = 0; $Length = length($InString);
    $FirstPiece = 1;
    for ($Index = 0; $Index <= $Length; $Index += $WrapLength) {
      if (($Index + $WrapLength) < $Length) {
	$StringPiece = substr($InString, $Index, $WrapLength);
      }
      else {
	# Last piece of the string...
	$StringPiece = substr($InString, $Index, $WrapLength);
      }
      if ($FirstPiece) {
	$FirstPiece = 0;
	$OutString = $StringPiece;
      }
      else {
	$OutString .= "${WrapDelimiter}$ {StringPiece}";
      }
    }
  }
  return $OutString;
}

1;

__END__

=head1 NAME

TextUtil - Text processing utilities

=head1 SYNOPSIS

use TextUtil;

use TextUtil qw(:all);

=head1 DESCRIPTION

TextUtil module provides method to process text lines and words:

    AddNumberSuffix - Add number suffix
    ContainsWhiteSpaces - Does string contains white spaces
    GetTextLine - Get next line of text
    IsInteger - Is string an integer
    IsPositiveInteger - Is string an positve integer
    IsFloat - Is string a float
    IsNotEmpty - Is sring not empty
    JoinWords - Join words using specified delimiter
    QuoteAWord - Quote a word
    RemoveLeadingWhiteSpaces - Remove leading white spaces
    RemoveTrailingWhiteSpaces - Remove trailing white spaces
    RemoveLeadingAndTrailingWhiteSpaces - Remove white spaces

=head1 METHODS

=over 4

=item B<AddNumberSuffix>

    $NumberWithSuffix = AddNumberSuffix($IntegerValue);

Returns number with appropriate suffix: 0, 1st, 2nd, 3rd, 4th, and so on.

=item B<ContainsWhiteSpaces>

    $Status = ContainsWhiteSpaces($TheString);

=item B<GetTextLine>

    $Line = GetTextLine(\*TEXTFILE);

Read next line from an already opened text file and take out any carriage return.
And return it as a string. NULL is returned for EOF.

=item B<IsInteger>

    $Status = IsInteger($TheString);

=item B<IsPositiveInteger>

    $Status = IsPositiveInteger($TheString);

=item B<IsFloat>

    $Status = IsFloat($TheString);

=item B<IsNotEmpty>

    $Status = IsNotEmpty($TheString);

=item B<IsNumerical>

    $Status = IsNumerical($TheString);

=item B<JoinWords>

    $JoinedWords = JoinWords($Words, $Delim, $Quote);

Join different words using delimiter and quote parameters. And return it
as a string.

=item B<QuoteAWord>

    $QuotedWord = QuoteAWord($Word, $Quote);

Return a quoted string based on $Quote flag.

=item B<RemoveLeadingWhiteSpaces>

    $OutString = RemoveLeadingWhiteSpaces($InString);

=item B<RemoveTrailingWhiteSpaces>

    $OutString = RemoveTrailingWhiteSpaces($InString);

=item B<RemoveLeadingAndTrailingWhiteSpaces>

    $OutString = RemoveLeadingAndTrailingWhiteSpaces($InString);

=item B<WrapText>

    $OutString = WrapText($InString, [$WrapLength, $WrapDelimiter]);

=back

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 COPYRIGHT

Copyright (C) 2004-2005 Manish Sud. All rights reserved.

This library is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 2.1 of the License, or (at your option)
any later version.

=head1 SEE ALSO

=cut
