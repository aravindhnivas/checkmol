checkmol/matchmol, (C) Norbert Haider, University of Vienna, 2003-2013
norbert.haider@univie.ac.at

A utility program for functional group analysis of chemical structures
and for structure/substructure comparison of two or more molecules.


===============================================================================

This software is published under the terms of the GNU General Public 
License (GPL). For a detailed description of this license, see 
http://www.gnu.org/copyleft/gpl.html

If invoked as "checkmol", this program reads 2D and 3D molecular structure
files in various formats (Tripos Alchemy "mol", Tripos SYBYL "mol2", MDL "mol")
and describes the molecule either by its functional groups or by a set of 
descriptors useful for database pre-screening (number of rings, sp2-hybridized
carbons, aromatic bonds, etc.).

If invoked as "matchmol", the program reads two individual 2D or 3D molecular
structure files in various formats (see above) and checks if the first molecule
(the "needle") is a substructure of the second one (the "haystack").
"Haystack" can also be a MDL SD-file (containing multiple MOL files);
if invoked with "-" as file argument, both "needle" and "haystack" are
read as an SD-file from standard input, assuming the first entry in
the SDF to be the "needle"; output: entry number + ":F" (false) or ":T" (true)

===============================================================================


Quick installation guide


Compile with fpc (Free Pascal, see http://www.freepascal.org), using
the -Sd or -S2 option (Delphi mode; IMPORTANT!)

example for compilation (with optimization) and installation:

fpc checkmol.pas -S2 -O3 -Op3 

as "root", do the following:

cp checkmol /usr/local/bin    (or any other directory in your path)
cd /usr/local/bin
ln checkmol matchmol          (ATTENTION: a symbolic link does not work!)


Note that checkmol and matchmol are the same executable, but the program
behaves differently depending on the name it was invoked with. Of course,
you can also copy "checkmol" to "matchmol", but then it takes twice as
much disk space (under Windows, this is the only possibility, as there are 
no hard links available there).

