###
### COSEG
###

Background
==========

  A program to identify repeat subfamilies using significant
 co-segregating ( 2-3 bp ) mutations.

   This program is derived from three C programs and several perl
 scripts written by Alkes Price as part of an analysis of Alu 
 elements in the human genome ( Whole-genome analysis of Alu repeat
 elements reveals complex evolutionary history, Alkes L. Price,
 Eleazar Eskin, and Pavel Pevzner, 2004 Genome Research ). The program
 was first adapted for use with other repeat families and then
 extended to support consideration of three co-segregating mutations
 using Alkes statistical model.  In 2008 with the help of Andy Siegel
 an alternative statistical model was developed and the codebase
 repackaged into the single source file.

 There are a few caveats:  

  * The input sequences must be full length alignments to a *single*
    reference sequence.

  * The longer the sequence the harder this problem is to solve.  Shorter
    families or short subregions are recommended for this process.

 Current Authors: Robert Hubley <rhubley@systemsbiology.org>
                  Andy Siegel
                  Arian Smit

 Original Code Authors: Alkes Price

 The original code can be found here:
            http://www.cse.ucsd.edu/~ppevzner/download/alucode.tar.gz


Installation
============

  1. Build COSEG c program
         make

  2. Edit preprocessAlignments.pl
      
       Change line 114:

              use lib "/usr/local/RepeatMasker";
  
       To point to your RepeatMasker installation.

  3. OPTIONAL:  If you plan to run the experimental 
                script refineConsSeqs.pl you will need
                to edit the script.

        Change line ~110:

              my $RepeatModelerDir = "/usr/local/RepeatModeler";

        To point to your RepeatModeler location.


 
ALU Example Run:

  1. Given the sample ALU dataset provided by Alkes Price in the
     original distribution of his code ( files have been renamed ):

            ALU.seqs: Human sequence data from alignment to AluSx consensus )
            ALU.ins : Insertion sequences from the alignments 
            ALU.cons: AluSx consensus 
                      Note: that positions 116-135 (inclusive) are 
                      lowercase while all remaining positions are uppercase.
                      This is how we encode positions within the consensus
                      that have lower quality and should not be considered
                      in this analysis.

  2. Run analysis:

         runcoseg.pl -k -d -m 50 -c ALU.cons -s ALU.seqs -i ALU.ins

         NOTE: In this run we must tell coseg to treat lower case
               characters in the consensus as a blacklist designation
               by using the "-d" flag.  Also we are using the 
               original pValue calculation as defined by Alkes et al.
               by specifying the "-k" flag to closely approximate the
               functioning of Alkes original program.
 
  3. Create png/svg files
 
       Coseg outputs a visualization of the subfamily tree in 
       both GraphViz and SVG formats.  The SVG format can be viewed
       directly in any web browser.  The GraphViz format needs to 
       be converted to a graphics format ( ie. PNG ) using
       the GraphViz command line tools ( http://www.graphviz.org/ )
       or loaded using the graphViz webapp ( https://mdaines.github.io/viz.js/ ).

       The SVG output may be resized by changing the first line
       in the SVG file to reflect the output size you would like.

       For example the default is to display the entire graph in
       the current size of the web browser window:

          <svg height="100%" width="100%" viewbox="0 0 90.625 146.25">

       Changing height/width to:

          <svg height="147" width="91" viewbox="0 0 90.625 146.25">
 
       Would display the graph using the full scale.  The viewbox
       values give you the bounding box of the graph at full scale.


Output Files  
============

  The following files are named after the aligned sequence file
  name as a prefix.

  *.log -- A log of mutation sites found and the order of the
           clustering.
  *.subfamililies.seq -- name, count, P-value and consensus sequence of
         each subfamily found by our algorithm.  (For subfamilies not in the
         original scaffold, we also include in parentheses the P-value of
         the scaffold subfamily from which it is derived).
  *.assign -- for each of the elements, lists the
         subfamily to which the algorithm has assigned it.
  *.tree.svg -- evolutionary tree of the subfamilies, in SVG format.
  *.tree.viz -- evolutionary tree of the subfamilies, in GRAPHVIZ format.


Input Files 
===========

  *.cons -- Consensus File:
    The consensus used to generate the aligned seqeunce file.  This is
    in fasta format.

  *.seqs -- Aligned Sequences File:
    Near full length aligned sequences to the consensus file ( one per line ).
    Deletions are represented by 1 or more "-" characters.  Insertions of 
    any length are represented by a single "+" character.

  *.ins -- Inserts File:
    Inserts pulled from the alignments are stored in this file.  The 
    aligned sequence file ordering is used in this file.  Each insert
    is encoded as:

            [<position>:<sequence>][<position>:<sequence>]..

    ie.
             3:AA 10:AGTTTA

  A script is provided to build these files automatically given either
  wu-blast or cross_match alignment files -- see below.



Running using your own data
===========================

  1. Cross_match a reference sequence against a genome or database. 

       cross_match line1copies consensus -M 25p41g.matrix 
                   -gap_init -25 -gap_ext -5 -minscore 200 
                   -minmatch 6 -alignments -bandwidth 50 -word_raw > LINE1

       The example file LINE1 included in this distribution was created
       using the command line above and can be used directly in the following
       steps.

  2. Determine consensus range to use for analysis ( ie. 298 - 797 bp )

  3. Create input files to alkes programs:

       preprocessAlignments.pl -maxEdgeGap 10
                                -minConsRange 298
                                -maxConsRange 797
                                -alignments LINE1

       This will create 3 new files: LINE1.seqs
                                     LINE1.ins
                                     LINE1.cons
      
       NOTE: Use the -w flag to preprocessAlignments if you use WUBlast
             to perform the alignments.

  4. Run analysis:

       runcoseg.pl -t -m 50 -c LINE1.cons -s LINE1.seqs 
                   -i LINE1.ins

       NOTE: In this examplewe use the "-t" flag to indicate we want to use
             3 bp co-segregating mutations as well as 2bp co-segregating 
             mutations when developing subfamilies.  Also we use the
             new pValue calculation ( default ) by leaving out the "-k"
             flag.

  5. Open up a web browser and point it at the file LINE1.seqs.tree.svg.
     Most browsers support zooming in on svg files.  If you want to render 
     the SVG file larger by default simply edit the *.svg file and change the
     line:

          <svg height="100%" width="100%" viewbox="0 0 90.625 146.25">

     to reflect a fixed size for the graph.  The "viewbox" values give
     the absolute size of the drawing so a 1:1 scale would be:
    
          <svg height="145.25" width="90.625" viewbox="0 0 90.625 146.25">

  6. If you would like to produce a SVG file without node label or 
     using either the divergence, P-value or subfamily size as the 
     label, simply rerun ./postprocess.pl without the "-l" flag or
     with "-l div", "-l pv", "-l c" accordingly.


Utility Programs
================

refineConsSeqs.pl  

Coseg uses a rather crude method of building consensus sequences
for each subfamily it finds.  This script makes use of the refiner
script in the RepeatModeler package to build and refine consensus
sequences based on subfamily members assigned by coseg.

Usage:

  1. From the directory where your coseg results can be
     found run:

          refineConsSeqs.pl -subConsFile mycosegrun.subfamilies.seq

     Where "mycosegrun" is the prefix of the coseg run.


extractSubSeqs.pl

A utility script to extract the sequences for a given coseg subfamily.
This uses the coseg input sequence file along with the *.assign file
to determine which sequences belong to the queried subfamily number.

Usage:




Version History
---------------
  0.2.3:
  * IUB codes in input sequences caused the code to segfault.
    Coseg will now randomly choose a nucleotide to substitute
    each time it encounters one in in the input sequence. It
    will also inform the user when doing so.  Thanks to David
    Ray for reporting this and suggesting the fix.
  * Fixed the svg tag so that the files
    will directly load in HTML5 web browsers.
  * Calculation of divergence has been improved.
    We now use kimura substition distance with CpG 
    site accounting modifications instead of the
    mixed substition and indel calculation.

  0.2.2:  
  * Create a *.svg graph file without the need
    to download/use GraphViz.  The layout is
    handled by an adaptation of algorithm 
    developed by Atze van der Ploeg. The SVG
    file produced supports various labeling 
    options and subfamily details displayed
    when a node is hovered over.
  * Changed the default colormap for the graph
    output.  Now warm colors denote more diverged
    subfamilies in the tree while cooler colors
    represent younger subfamilies.  To restore
    the original color scheme use the new "-o" 
    flag to coseg.
  * Added parameter to control the minimum distance
    between diagnostic sites. Now the user can override
    the historic value of 10 using the -u flag.
  * Improved error reporting when there is a mismatch
    between an individual sequence length and the 
    consensus length in the input files.
  * Fixed a bug that caused coseg to segfault.
  * Added experimental script refineConsSeqs.pl.  This
    script uses the RepeatModeler application to build
    and refine the consensus sequences for each 
    subfamily.

  0.2.1:  
  * Improved code documentation
  * Single mutation significance cutoff ( SIGMATHRESH ) was 
    pre-calculated for Alkes Alu analysis and hardcoded.  This
    version calculates the correct sigma cutoff using the length
    of the input sequence.
  * Fixed bug with implementation of Siegel's pValue
    calculation which caused a segfault -- found by Neal Platt.
  * Switched default pvalue method to Andy Siegel's method and
    provided a new "-k" switch to use Alkes Price's method.
  * Fixed bug where the program was exiting when calculations 
    fell below the precision of the machine ( epsilon ). Message
    given was "Below epsilon..." and the runcoseg.pl script
    moved on even though coseg failed.

