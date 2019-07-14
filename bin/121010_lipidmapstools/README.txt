INTRODUCTION

LIPID MAPS Tools package provides various command line Perl scripts to generate
SD files containing structure and ontological data for various lipid categories.

INSTALLATION

1. Add <YOUR LIPID MAPS TOOLS DIR>/bin to your PATH environment variable.

      C shell example:

          set path = "$path <YOUR LIPID MAPS TOOLS DIR>/bin"

      Bourne shell example:

          PATH = "$PATH:<YOUR LIPID MAPS TOOLS DIR>/bin"

      Windows:

          o Right click on "My Computer" icon and select properties. And
             Select "Environment Variables" button under Advanced tab.

          o Select path under "System Variables" and click on Edit button. And 
             system variable "path" to include LIPID MAPS TOOLS bin directory.

2. And you're all set to try out the command line scripts.

CAVEATS

On some systems, command line scripts may need to be invoked using
"perl -S <ScriptName.pl>" instead of "<ScriptName.pl>". Example:

      % perl -S GLStrGen.pl

FEEDBACK

webmaster@lipidmaps.org

AUTHOR

Manish Sud

CONTRIBUTOR(S)

Eoin Fahy

COPYRIGHT

Copyright (C) 2006-2012. The Regents of the University of California. All Rights Reserved.

LICENSE

Modified BSD License
