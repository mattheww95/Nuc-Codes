nuc-codes

A simple help you convert those pesky gene-coordinates to genome coordinates.


Quickstart:
Install globally or in an environment, I don't care what matters is that you use 
python >3.7

>pip install .

Let it install
Then follow the help messages by entering nuc-codes help.

Then to initialize a *project*
>nuc-codes init reference_fasta.fa reference_gff.gff project_name

You will then not have to specify a fasta and gff each time you run nuc-codes.

Then to list attributes you can query your sequence for, which come from the gff file
enter:
>nuc-codes *project* attributes

This will display all of the attribute keys in your gff file for you to query, it is best
to already be familiar with what values those keys may hold in your gff.

Then to convert positions enter:
>nuc-codes *project* gff_attribute gff_value position_1 etc.

You can enter multiple positions or just one, the full help message is as follows:

>nuc-code help options
>A program to help Cody and hopefully make him help me more.
>
>Initialize Program to save a gff and reference sequence.
>nuc-code init <reference_fasta> <reference_gff> <project_name>
>
>Check possible attributes to use for gff extraction.
>nuc-code <project_name> attributes
>e.g. nuc-code <project_name> attributes
>
>Query a gene positions, you can query one or multiple conditions
>nuc-code <project_name> <gff attribute key> <gff attribute key value> <Amino acid positions>
>e.g. nuc-code SARS Name orf1ab 4 5 6
>
>Help note
>If no attribute is displayed check that the values you are using exist in your gff.
>Note of Importance:
>The Phase attribute of GFF file is not included in the tabulation of nucleotides.
>It is also assumed that the GFF being used is base 1 indexed.

Important note:
Phase of the codon within the GFF file is not used in tabulation of codons


A little bit about the program:
nuc-codes actually stands for Nuclear Cody, not nucleotide coordinates. Just as how in
my mind gff does not stand for General Feature Format but instead Good Friends File.

Why nuclear cody you might ask:
Well during this whole COVID-19 debacle we were tasked with tallying various mutations from
variants of concern (This was awful). There multiple different reference sequences and
gff files for the reference sequence, and most literature would only cite gene coordinates
which were not always the easiest to convert.

Needless to say this made Cody go Nuclear! Everytime he saw gene coordinates, it was like watching two beryllium
spheres over a plutonium-gallium core. Which was understandable, as it is an incredibly frustrating excercise
and Cody is incredibly talented and would go through a rigorous process to make sure we always identified the right SNP.

Unfortunately I could not rely on Cody forever, as his attention can be fleeting and without his aid I was useless
. So I did the only thing a Python developer could do.... dedicate a simple utility to them and beg them for aid.

This is the story of nuc-codes, let us hope it is enough to win his future efforts.