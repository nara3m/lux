package FileUtil;
#
# $RCSfile: FileUtil.pm,v $
# $Date: 2007/09/05 15:55:58 $
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
@EXPORT = qw(CheckFileType ConvertCygwinPath ExpandFileNames GetUsageFromPod ParseFileName);
@EXPORT_OK = qw();
%EXPORT_TAGS = (
		all  => [@EXPORT, @EXPORT_OK],
		);

# Check to see path contains cygdrive and convert it into windows path...
sub ConvertCygwinPath {
  my($Path) = @_;
  my($NewPath, $OSName);

  $NewPath = $Path; $OSName = $^O;
  if ($OSName =~ /cygwin/i && $Path =~ /cygdrive/i ) {
    my(@PathParts) = split "\/", $Path;
    my($Drive) = $PathParts[2];
    shift @PathParts; shift @PathParts; shift @PathParts;
    $NewPath = join "\/", @PathParts;
    $NewPath = $Drive. ":\/" . $NewPath;
  }
  return $NewPath;
}

# Based on the file name extension, figure out its type.
sub CheckFileType {
  my($FileName, $FileExts) = @_;
  my($Status, @FileExtsList, $Index, $Ext);

  $Status = 0;
  @FileExtsList = split " ", $FileExts;
  for $Index (0 .. $#FileExtsList) {
    $Ext = $FileExtsList[$Index];
    if ($FileName =~ /(\.$Ext)$/i) {
      $Status = 1;
    }
  }
  return ($Status);
}

# For any directory name specified in Files array reference parameter, generate
# all files names in that directory which correspond to the extensions specified.
# And combine them with rest of the names. Return the expanded list.
sub ExpandFileNames {
  my($Files, $FileExts) = @_;
  my($FileName, @FilesList, $Index, $Delimiter);

  for ($Index = 0; $Index < @$Files; $Index++) {
    $FileName = @$Files[$Index];
    $Delimiter = "\/";
    if ($FileName =~ /\\/ ) {
      $Delimiter = "\\";
    }
    if (-d $FileName && $FileExts) {
      my($DirName) = $FileName;
      # glob doesn't appear to work during command line invocation from Windows...
      # So, use opendir to make it work...
      #
      # push @FilesList,  map {glob("$DirName/*.$_")} split " ", $FileExts;
      if (opendir DIRNAME, $DirName) {
	# Collect files with specified file extensions...
	my($FileExtsPattern) = join "|", map { "\.$_" } split " ", $FileExts;
	push @FilesList, grep /$FileExtsPattern/, map { "$DirName$Delimiter$_"  } readdir DIRNAME;
	closedir DIRNAME;
      }
      else {
	warn "Warning: Ignoring directory $DirName: Couldn't open it: $!";
      }
    }
    elsif ($FileName =~ /\*/) {
      # Filenames are not expanded during command line invocation from Windows...
      my($FileDir, $Name, $FileExt) = ParseFileName($FileName);
      if (opendir FILEDIR, $FileDir) {
	if (length($Name) == 1) {
	  push @FilesList, grep /\.$FileExt/, map { "$FileDir$Delimiter$_"  } readdir FILEDIR;
	}
	else {
	  push @FilesList, grep /\.$FileExt/, grep /$Name/, map { "$FileDir$Delimiter$_"  } readdir FILEDIR;
	}
	closedir FILEDIR;
      }
      else {
	warn "Warning: Ignoring files $FileName: Couldn't open directory $FileDir: $!";
      }
    }
    else {
      push @FilesList, $FileName;
    }
  }
  return @FilesList;
}

# Get Usage from Pod...
sub GetUsageFromPod {
  my($Usage, $ScriptPath);

  ($ScriptPath) = @_;
  $Usage = "Script usage not available: pod2text/pod2text.bat failure in your Perl installtion\n";

  # Get pod documentation: try pod2text first followed by perdoc.bat in case it fails to
  # to handle ActiveState Perl; otherwise, give up and just return an empty string...
  my($PodStatus) = (open CMD, "pod2text $ScriptPath|") ? 1 : ((open CMD, "pod2text.bat $ScriptPath|") ? 1 : 0);
  if (!$PodStatus) {
    return $Usage;
  }
  my(@LineWords);
  PODLINE: while (<CMD>) {
    chomp;
    if (/^SYNOPSIS/) {
      $_ = <CMD>;
      chomp;
      s/^ +//g;
      (@LineWords) = split / /;
      $Usage = qq(Usage: $LineWords[0] [-options]... );
      shift @LineWords;
      $Usage .= join " ", @LineWords;
      $Usage .= qq(\n);
    }
    elsif (/^DESCRIPTION/) {
      $Usage .= qq(Description:\n);
      while (<CMD>) {
	if (/^OPTIONS/) {
	  $Usage .= qq(Options:\n);
	  while (<CMD>) {
	    # Collect options lines...
	    if (/^[ ]+\-/) {
	      # Option specifier line: Put back in <> which pod2text replaced
	      # to **
	      my($OptionLine) = qq($_);
	    OPTIONLINE: while (<CMD>) {
		if (/^(    )/) {
		  $OptionLine .= qq($_);
		}
		else {
		  $OptionLine =~ s/\*(([a-zA-Z0-9])|(\[)|(\#)|(\"))/"\<" . substr($&, -1, 1)/e;
		  $OptionLine =~ s/(([a-zA-Z0-9])|(\])|(\#)|(\"))\*/substr($&, 0, 1) . "\>"/e;
		  $Usage .= qq($OptionLine$_);
		  last OPTIONLINE;
		}
	      }
	    }
	    else {
	      $Usage .= qq($_);
	    }
	  }
	}
	else {
	  # Collect description lines...
	  $Usage .= qq($_);
	}
      }
    }
  }
  close CMD;

  # Take out * which pod2text puts in for <>
  $Usage =~ s/\*(([a-zA-Z0-9;#-])|(\")|(\()|(\[)|(\.))/substr($&, -1, 1)/eg;
  $Usage =~ s/(([a-zA-Z0-9;#-])|(\")|(\))|(\])|(\.))\*/substr($&, 0, 1)/eg;

  return $Usage;
}

# Split full file name into directory path, file name, and the extension.
sub ParseFileName {
  my($FullName) = @_;
  my($FileDir, $FileName, $FileExt, @FullFileNameParts, @FileNameParts, $Delimiter);

  $Delimiter = "\/";
  if ($FullName =~ /\\/ ) {
    $Delimiter = "\\";
    $FullName =~ s/\\/\//g;
  }
  $FileDir = ""; $FileName = ""; $FileExt = "";
  @FullFileNameParts = (); @FileNameParts = ();

  @FullFileNameParts = split "\/", $FullName;
  @FileNameParts = split /\./, $FullFileNameParts[$#FullFileNameParts];

  # Setup file dir...
  if (@FileNameParts == 1) {
    $FileDir = $FullName;
  }
  else {
    if (@FullFileNameParts == 1) {
      $FileDir = "\.$Delimiter";
    }
    else {
      pop @FullFileNameParts;
      $FileDir = join $Delimiter, @FullFileNameParts;
    }
  }

  # Setup file name and ext...
  if (@FileNameParts == 2) {
    $FileName = $FileNameParts[0];
    $FileExt = $FileNameParts[1];
  }
  return ($FileDir, $FileName, $FileExt);
}

1;

__END__

=head1 NAME

FileUtil - File processing utilities

=head1 SYNOPSIS

use FileUtil;

use FileUtil qw(:all);

=head1 DESCRIPTION

FileUtil module provides methods for processing files.

    CheckFileType - Figure out file type using its extension
    ConvertCygwinPath - Convert CygWin path into Windows
    ExpandFileNames - Retrieve full file names
    GetUsageFromPod - Retrieve script usage from pod
    ParseFileName - Split file name into dir, name and ext

=head1 METHODS

=over 4

=item B<CheckFileType>

    $Status = CheckFileType($FileName, $FileExts);

Based on $FileExts, decide type of $FileName and return its value as a flag.

=item B<ConvertCygwinPath>

    $NewPath = ConvertCygwinPath($Path);

Check to see $Path contains cygdrive and convert it into Windows path...

=item B<ExpandFileNames>

    (@FilesList) = ExpandFileNames(\@Files, $FileExts);

For a directory name in @Files, generate all file names which correspond to extensions
in $FileExts and return these values as a list.

=item B<GetUsageFromPod>

    $ScriptUsage = GetUsageFromPod($AbsoluteScriptPath);

Get pod documentation from script using pod2text first followed by perdoc.bat
in case it fails to to handle ActiveState Perl and then convert into a usage
description; otherwise, give up and just return an emty string.

=item B<ParseFileName>

    ($FileDir, $FileName, $FileExt) = ParseFileName($FullFileName);

Splits $FullFileName into directory name, file name, and extension. $FileDir is
NULL for absence of directory name in $FullFileName.

=back

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 COPYRIGHT

Copyright (C) 2004-2005 Manish Sud. All rights reserved.

This library is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 2.1 of the License, or (at your option)
any later version.

=cut
