#!/usr/bin/perl
=head1 NAME

  postprocess.pl - Post processes the coseg results 

=head1 SYNOPSIS

  postprocess.pl -s <sequences file> -i <insertions file> 
                 -c <consensus file> [-l sub|div|c|pv]


  options:  
      -l : Place a label inside each node:
               sub = Subfamily Number
               div = Percent divergence
                 c = Subfamily size
                pv = P-Value for subfamily
  
=cut

use strict;
use Getopt::Long;
use Data::Dumper;


## TODO: Create cytoscape *.sif file:
##              Node#<tab>-<tab>Node#
##              ..
##       Create cytoscape *.noa files ( one for each parameter ):
##              pVal (java.lang.String)
##              Node# = Value
##              ..       
## Could also create an xgmml file -- which has all these things together:
##     http://www.cs.rpi.edu/~puninj/XGMML/
##<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
##<graph label="cyto.sif.1" xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:xlin
##k="http://www.w3.org/1999/xlink" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syn
##tax-ns#" xmlns:cy="http://www.cytoscape.org" xmlns="http://www.cs.rpi.edu/XGMML"
## >
##  <att name="documentVersion" value="1.1"/>
##  <att name="networkMetadata">
##    <rdf:RDF>
##      <rdf:Description rdf:about="http://www.cytoscape.org/">
##        <dc:type>Protein-Protein Interaction</dc:type>
##        <dc:description>N/A</dc:description>
##        <dc:identifier>N/A</dc:identifier>
##        <dc:date>2008-07-24 11:08:00</dc:date>
##        <dc:title>cyto.sif.1</dc:title>
##        <dc:source>http://www.cytoscape.org/</dc:source>
##        <dc:format>Cytoscape-XGMML</dc:format>
##      </rdf:Description>
##    </rdf:RDF>
##  </att>
##  <att type="string" name="backgroundColor" value="#ffffff"/>
##  <node label="237" id="-249">
##    <att type="real" name="divergence" value="0.14"/>
##    <att type="integer" name="p" value="2"/>
##    <att type="string" name="subfamily" value="sub237"/>
##    <att type="real" name="pVal" value="0.1"/>
##    <att type="integer" name="count" value="94"/>
##    <att type="string" name="mutations" value="235:t-&amp;amp;gt;g"/>
##    <att type="string" name="canonicalName" value="237"/>
##  </node>
##  <edge label="237 (-) 236" source="-249" target="-246">
##    <att type="string" name="XGMML Edge Label" value="237 (-) 236"/>
##    <att type="string" name="interaction" value="-"/>
##    <att type="string" name="canonicalName" value="237 (-) 236"/>
##    <graphics width="1" fill="#0000ff" cy:sourceArrow="0" cy:targetArrow="0" cy:sourceArrowColor="#000000" cy:targetArrowColor="#000000" cy:edgeLabelFont="Default-0-10" cy:edgeLineType="SOLID" cy:curved="STRAIGHT_LINES"/>
##  </edge>
##</graph>


#
# Option processing
#  e.g.
#   -t: Single letter binary option
#   -t=s: String parameters
#   -t=i: Number paramters
#
my @getopt_args = (
    '-s=s', 
    '-c=s',
    '-l=s',
    '-i=s'
);

my %options = ();
Getopt::Long::config("noignorecase", "bundling_override");
unless (GetOptions(\%options, @getopt_args)) {
    usage();
}

sub usage {
  exec "pod2text $0";
  exit;
}

usage() if ( ! $options{'s'} || ! $options{'i'} || ! $options{'c'} );


# input file listing insertion values of each element
my $insertSourceFile = $options{'i'};  

# input file of assignments to all subfamilies 
my $assignFile = $options{'s'} . ".assign"; 

# Consensus variation output from program.
my $subFile = $options{'s'} . ".subfamilies";   

# Original consensus file
my $consFile = $options{'c'};   

# Main log file
my $logFile = $options{'s'} . ".log";

# Graph file
my $vizFile = $options{'s'} . ".tree.viz";

# Graph SVG file
my $svgFile = $options{'s'} . ".tree.svg";

# resolved insertions
my $tmpInsertConsFile = $options{'s'} . ".subfamilies.ins";   

# Final sequence file for each subfamily with insertions resolved.
my $seqFile = $options{'s'} . ".subfamilies.seq";   


print "Postprocessing coseg results\n";

##
##  Read in mutation information
##
my @parent = ();
open MUT, "<$logFile" ||
   die "Error: Could not open file $logFile\n\n";
my @mutations = ();
while ( <MUT> )
{
  if ( /Building subfamily\s+(\d+)\s+\(parent\s+(\d+), logpvalue.*:\s+(.*)/ )
  {
    my $subfam = $1;
    my $muts = $3;
    $parent[$1] = $2;
    my $mutString = "parent=$2 muts=";
    while ( $muts =~ /pos\s+(\d+)\s+(\S)\s+to\s+(\S)/ig )
    {
      $mutString .= " $1:$2/$3";
    } 
    $mutations[$subfam] = $mutString;
  }
}
close MUT;
 
 
## 
## Read in pvalues
##
my @pvalues = ();
# First the scaffold pvalues
open SUB, "<$subFile" ||
   die "Error: Could not open file $subFile\n\n";
while ( <SUB> )
{
 if ( /Subfamily\s+(\d+):.*(mstLogPValue|logpvalue)\s+([-\d\.]+)/ ) 
 {
   my $aa = abs($3*.434294481);
   my $b = int($aa);
   my $c = $aa - $b;
   my $d = exp(-$c/.434294481);
   my $x; 
   my $y;
   if($d >= 0.95)
   {
     $x = 1;
     $y = -$b;
   }
   else
   {
     $x = int(10*$d + 0.5);
     $y = -($b+1);
   }
   $pvalues[$1] = "${x}e$y";
 } 
}
close( LISTFILE );


##
## Build the insert file
##
open TMPINSERT, ">$tmpInsertConsFile" || 
   die("Error: Could not open file for writing: $tmpInsertConsFile!\n\n");
open(SUB,"$subFile") || 
   die "Error: Could not open file $subFile\n\n";
my $s = 0;
my $consensus;
my $consensuscount;
my $n;
my @count = ();
while( <SUB> )
{
  #Subfamily 4: count 739 mutrate 0.027/0.023 mstLogPValue -413.400874
  #Subfamily 58: count 129 mutrate 0.206/0.200 sigma 17.839535, logpvalue -149.335610
  #Subfamily 26: count -1614591685 mutrate 0.114/0.087 sigma 10.109653, logpvalue -42.332037
  #if ( /^Subfamily\s+(\d+):\s+count\s+\d+\s+mutrate\s+([\d\.]+)\/([\d\.]+)/ )
  if ( /^Subfamily\s+(\d+):/ )
  { 
    print TMPINSERT $_;
    next;
  }

  #next unless ( /^\d+\:/ );
  if ( ! /^\d+\:/ )
  {
    print "Hmmm: $_\n";
  }

  chomp;
  my @array = split / +/;
  for( my $l=0; $l<=$#array; $l++ )
  {
    print TMPINSERT "$array[$l]";
    # If "+" appears in the a subfamily entry
    if( $array[$l] =~ /(\d+):\+/ )
    {
      # NOW, print insertion consensus for subfamily $s position $x

      # The insertion position
      my $insPos = $1;
 
      print "  Found an insertion site: $array[$l]\n";

      # The assignment file and the insertion sequences files
      #print "Opening $insertSourceFile\n";
      open INSERTS,"<$insertSourceFile" || 
          die "Error: Could not open file $insertSourceFile\n\n";
      open ASSIGNS,"<$assignFile" || 
          die "Error: Could not open file $assignFile\n\n";

      # Clear consensus variables
      $consensus = "~~~";
      $consensuscount = 0;
      $n = 0;

      print "    Calculating consensus for insertion site:";
      while(<ASSIGNS>)
      {
        print "." if ( ( $n % 1000 ) ==  0 );
        if ( /(\d+)\s+(\d+)/ ) 
        {
          my $lineb = <INSERTS>;
          if ( $1 != $n ) 
            { print("Problem with assign file $1 != $n\n"); exit(1); }
          if ( $2 != $s ) 
            { $n++; next; }


          chomp $lineb;
          my @insArray = split(/ +/,$lineb);
          for(my $l4=0; $l4<=$#insArray; $l4++)
          {
            if ( $insArray[$l4] =~ /(\d+):(\S+)/ )
            {
              my $pos = $1; 
              my $seq = $2;
              if( $pos != $insPos ) { next; }
              $count[$s][$insPos]{$seq}++;
              if($count[$s][$insPos]{$seq} > $consensuscount)
              {
                $consensuscount = $count[$s][$insPos]{$seq};
                $consensus = $seq;
              }
            }
          }
          $n++;
        }
      }
      print "\n";
      print "    Consensus = $consensus\n";
      if($consensus =~ /\~/) 
         { die "OOPS s=$s x=$insPos consensus=$consensus\n"; }
      print TMPINSERT "$consensus";
      close( INSERTS );
      close( ASSIGNS );
    }
    print TMPINSERT " ";
  }
  print TMPINSERT "\n";
  $s++;
}
close( SUB );
close( TMPINSERT );

# next, build seqFile
open TMPINSERT,"<$tmpInsertConsFile" || 
    die "Error: Could not re-open file $tmpInsertConsFile\n\n";
open ASSIGNS,"$assignFile" || 
    die "Error: Could not re-open file $assignFile\n\n";
open(SEQ,">$seqFile") || die("COF");
open CONS,"<$consFile" ||
   die "Error: Could not open file $consFile\n\n";

my @ccount = ();
while(my $line = <ASSIGNS>)
{
  my @array = split(/ +/,$line);
  $ccount[$array[1]] += 1;
}
close(ASSIGNS);
my $maxCount = 0;
foreach my $count ( @ccount )
{
  $maxCount = $count if ( $maxCount < $count );
}

my $conSeq;
while(my $line = <CONS>)
{
  if($line =~ />/) { next; }
  chomp($line);
  $conSeq .= lc($line);
}
my $L = length( $conSeq );
my @array = split(//,$conSeq);
my @Sx = ();
for(my $x=0; $x<$L; $x++) { $Sx[$x] = $array[$x]; }
my $n=0;
my @this = ();
my @this2 = ();
my @array2 = ();
while(my $line = <TMPINSERT>)
{
  for(my $x=0; $x<$L; $x++) { $this[$x] = $Sx[$x]; $this2[$x] = $Sx[$x]; }
  $line = <TMPINSERT>;
  chomp($line);
  @array = split(/ +/,$line);
  my $length = @array;
  for(my $l=0; $l<$length; $l++)
  { 
    @array2 = split(/:/,$array[$l]);
    my $x = $array2[0];
    my $string = $array2[1];
    if($string =~ /\+/)
    {
      $string =~ s/\+//;
      $this[$x] .= $string;
      $this2[$x] .= "+";
    }
    else
    {
      $this2[$x] = $string;
      $string =~ s/-//;
      $this[$x] = $string;
    }
  }

  print SEQ ">subfamily$n count=$ccount[$n] pvalue=$pvalues[$n]";
  if ( $parent[$n] >= 0 )
  {
    print SEQ " ( parent_pvalue=$pvalues[$parent[$n]] )";
  }
  print SEQ " " . $mutations[$n];
  print SEQ "\n";

  for(my $x=0; $x<$L; $x++) { print SEQ "$this[$x]"; }
  printf SEQ ("\n");

  $n++;
}
close( TMPINSERT );
close( ASSIGNS );
close( SEQ );
close( CONS );
print("Done building $seqFile\n");
unlink( $tmpInsertConsFile );


## Convert the VIZ file to SVG
open VIZ, "<$vizFile" ||
   die "Could not open $vizFile\n";

open OUT ,">$svgFile" ||
   die "Could not create $svgFile\n";

#
#
#
my $minNodeHeight = 30;

# 
# Parse the Graphviz tree definition file *.viz
# and build a tree datastructure compatible with
# the layout algorithm.  Also determine minimum
# node size and graph height and scale nodes 
# accordingly.
#
my @nodes = ();
my $root = undef;
my $minHeight = 100000;
while ( <VIZ> )
{
  if ( /^\s*(\d+)\s+\[(.*)\];\s*$/ )
  {
    # Node data line
    my $nodeIdx = $1;
    my $nodeData = $2;
    my ( $label, $height, $color, $isScaffold );
    if ( $nodeData =~ /label="([^"]+)"/ )
    {
      # Convert "\n" used in Graphviz labels to a comma separated 
      # string for SVG's tooltip-like display. 
      $label = $1;
      $label =~ s/\\n/, /g;
    }
    if ( $nodeData =~ /height=([\d\.]+)/ ) 
    {
      # Coseg uses circles to represent nodes in both 
      # Graphviz and SVG.  Height = Width, so we only
      # need to record one.
      $height = $1;
      $minHeight = $height if ( $minHeight > $height );
    }
    if ( $nodeData =~ /color="([^"]+)"/ )
    {
      my $colorData = $1;
      if ( $colorData =~ /([\d\.]+),([\d\.]+),([\d\.]+)/ )
      {
        # Convert to HSL
        # GraphViz using HSV.  Families defined by double or 
        # tripple mutations are darker and more saturated than
        # nodes with only single point mutations:
        # 
        # Single:          H     S     V
        #                0-1.0  .15   .95
        #    In SVG HSL :  H     S     L
        #                0-216  58.7%  87.6%
        #
        # Double/Tripple:  H     S     V
        #                0-1.0  .9    .6
        #    In SVG HSL :  H     S     L
        #                0-216  82.1%  32.9%
        #
        if ( $2 == 0.15 )
        {
          $isScaffold = 0;
          $color = ($1*360) . ",58.7%,87.6%";
        }
        else
        {
          $isScaffold = 1;
          $color = ($1*360) . ",82.1%,32.9%";
        }
      }
    }
    $nodes[$nodeIdx] = { 'h' => $height, 
                         'w' => $height,
                         'label' => $label,
                         'isScaffold' => $isScaffold,
                         'color' => $color };
  }elsif ( /^\s*(\d+)\s+->\s+(\d+);\s*$/ )
  {
    # Relationship line
    push @{$nodes[$1]->{'children'}}, $nodes[$2];
    if ( ! defined $root )
    {
      $root = $1;
    }elsif ( $root == $2 )
    {
      die "Multiple definitions of root in *.viz file!\n";
    }
  }
}
close VIZ;
if ( $minHeight < $minNodeHeight  )
{
  my $scale = $minNodeHeight / $minHeight;
  for ( my $i = 0; $i <= $#nodes; $i++ )
  {
    $nodes[$i]->{'h'} *= $scale;
    $nodes[$i]->{'w'} *= $scale;
  }
}
&layout( $nodes[$root], 20, 20 );
my ( $boundingHeight, $boundingWidth ) = findExactHeightWidth( $nodes[$root], 0, 0);
#print OUT "<svg height=\"100%\" width=\"100%\" viewbox=\"0 0 " . ($boundingWidth+5) . " " . ($boundingHeight+5) . "\">\n";
print OUT "<svg xmlns=\"http://www.w3.org/2000/svg\" height=\"" . ($boundingHeight+5) . "px\" width=\"" . ($boundingWidth+5) . "px\" viewbox=\"0 0 " . ($boundingWidth+5) . " " . ($boundingHeight+5) . "\">\n";
printSVG( $nodes[$root], $options{'l'} );
print OUT "</svg>\n";
close OUT; 



############################################################################
# Support routines for use with non-layered tree layout algorithm
#

# Example input tree
# my $tree = {
#             'w'        => 70,
#             'h'        => 70,
#             'label'    => "sub3,c=100,pv=3e-18,div=0.311",
#             'color'    => "238,58.7%,87.6%",
#             'isScaffold' => 1,
#             'children' => [
#                             {
#                               'w'        => 40,
#                               'h'        => 40,
#                               'label'    => "sub1,c=83,pv=3e-18,div=0.311",
#                               'color'    => "19,58.7%,87.6%",
#                               'isScaffold' => 1,
#                               'children' => [
#                                               {
#                                                 'w'        => 10,
#                                                 'h'        => 10,
#                                                 'children' => []
#                                               }
#                               ]
#                             },
#             ]
#};
#

# NOTE: Depends on complete layout
sub findExactHeightWidth
{
  my $tree = shift;

  my $height = $tree->{'y'} + $tree->{'h'};
  my $width = $tree->{'x'} + $tree->{'w'};
  for ( my $i = 0 ; $i <= $#{ $tree->{'children'} } ; $i++ )
  {
    my ($retHeight,$retWidth) = &findExactHeightWidth( $tree->{'children'}->[ $i ] );
    $height = $retHeight if ( $retHeight > $height );
    $width = $retWidth if ( $retWidth > $width );
  }
  return ( $height, $width );
}

#
# Convert the post-layout tree into SVG for display
# 
sub printSVG
{
  my $tree = shift;
  my $showLabel = shift;

  print OUT "<circle cx=\""
      . ( $tree->{'x'} + ( $tree->{'w'} / 2 ) )
      . "\" cy=\""
      . ( $tree->{'y'} + ( $tree->{'h'} / 2 ) + 1 )
      . "\" r=\""
      . ( $tree->{'w'} / 2 )
      .  "\" stroke=\"black\" ";
  if ( $tree->{'isScaffold'} == 1 )
  {
    print OUT "stroke-width=\"5\" ";
  }else
  {
    print OUT "stroke-width=\"1\" "; 
  }
  print OUT "style=\"fill:hsl(" . $tree->{'color'} . ");\">\n";
  print OUT "<title>" . $tree->{'label'} . "</title></circle>\n";

  my $label = undef;
  if ( $showLabel eq "sub" ) {
    $label = $1 if ( $tree->{'label'} =~ /(sub\d+)/ );
  }elsif ( $showLabel eq "div" ) 
  {
    $label = $1 if ( $tree->{'label'} =~ /div=([\d\.]+)/ );
    $label = sprintf("%0.0f%", $label*100 );
  }elsif ( $showLabel eq "pv" )
  {
    $label = $1 if ( $tree->{'label'} =~ /pv=([\d\.Ee\-]+)/ );
  }elsif ( $showLabel eq "c" )
  {
    $label = $1 if ( $tree->{'label'} =~ /c=(\d+)/ );
  }
 
  if ( $label ne "" ) {
    print OUT "<text font-size=\"10\" alignment-baseline=\"middle\" " 
        . "text-anchor=\"middle\" x=\""
        . ( $tree->{'x'} + ( $tree->{'w'} / 2 ) )
        . "\" y=\""
        . ( $tree->{'y'} + ( $tree->{'h'} / 2 ) + 1 )
        . "\" fill=\"black\">"
        . $label
        . "</text>\n";
  }

  return if ( !exists $tree->{'children'} || !@{ $tree->{'children'} } );
  for ( my $i = 0 ; $i <= $#{ $tree->{'children'} } ; $i++ )
  {
    my $child          = $tree->{'children'}->[ $i ];
    my $parent_start_x = ( $tree->{'x'} + ( $tree->{'w'} / 2 ) );
    my $parent_start_y = ( $tree->{'y'} + $tree->{'h'} + 1 );
    my $child_end_x    = ( $child->{'x'} + ( $child->{'w'} / 2 ) );
    my $child_end_y    = $child->{'y'} + 1;
    my $height         = ( $child_end_y - $parent_start_y ) / 2;

    print OUT "<path d=\"M$parent_start_x $parent_start_y C $parent_start_x "
        . ( $parent_start_y + $height )
        . ", $child_end_x "
        . ( $child_end_y - $height )
        . ", $child_end_x $child_end_y\" stroke=\"black\" fill=\"transparent\"/>\n";
    &printSVG( $tree->{'children'}->[ $i ], $showLabel );
  }
}

#################################################################################
# This is a non-layerd, variable node size layout algorithm adapted 
# from:
#    The extended Reingold-Tilford algorithm as described in the paper
#    "Drawing Non-layered Tidy Trees in Linear Time" by Atze van der Ploeg
#
#       Ploeg, Atze. "Drawing non-layered tidy trees in linear time." 
#       Software: Practice and Experience 44.12 (2014): 1467-1484.
#       PDF: http://oai.cwi.nl/oai/asset/21856/21856B.pdf
#       Java Code: https://github.com/cwi-swat/non-layered-tidy-trees
#    
#    This code is in the public domain, use it any way you wish. 
#    A reference to the paper is appreciated!
# 
# A D3 adaptation by Chris Maloney may be found here:
#     https://github.com/Klortho/d3-flextree
#
#
# Tree Datastructure:
# NODE:
#    {  'w'   : Input width of the node
#       'h'   : Input height of the node
#       'x'   : Output X coordinate ( top-left of node bounding box )
#       'y'   : Output Y coordinate ( top-left of node bounding box ) 
#       'tl'  : Internal - left thread
#       'tr'  : Internal - right thread
#       'el'  : Internal - extreme left node
#       'er'  : Internal - extreme right node
#       'msel': Internal - sum of modifiers at the extreme left node
#       'mser': Internal - sum of modifiers at the extreme right node
#       'mod' : Internal - modifier
#       'prelim" : Internal - preliminary x position
#       'children' :  Input pointer to an array of children NODEs.
#     }
#

sub layout
{
  my $tree            = shift;
  my $neighborSpacing = shift;
  my $levelSpacing    = shift;
  &setLevelSpacing( $tree, $levelSpacing, 0 );
  &firstWalk( $tree, $neighborSpacing );
  my $minX = &secondWalk( $tree, 0, undef );
  &normalize( $tree, -$minX );
}

sub normalize
{
  my $tree = shift;
  my $adj  = shift;

  $tree->{'x'} += $adj;

  return if ( !exists $tree->{'children'} || !@{ $tree->{'children'} } );
  for ( my $i = 0 ; $i <= $#{ $tree->{'children'} } ; $i++ )
  {
    &normalize( $tree->{'children'}->[ $i ], $adj );
  }
}

# RMH: Set the level spacing
sub setLevelSpacing
{
  my $tree         = shift;
  my $levelSpacing = shift;
  my $currY        = shift;

  $tree->{'y'}   = $currY;
  $tree->{'mod'} = 0;
  return if ( !exists $tree->{'children'} || !@{ $tree->{'children'} } );
  for ( my $i = 0 ; $i <= $#{ $tree->{'children'} } ; $i++ )
  {
    &setLevelSpacing( $tree->{'children'}->[ $i ],
                      $levelSpacing, $currY + $tree->{'h'} + $levelSpacing );
  }
}

sub firstWalk
{
  my $tree            = shift;
  my $neighborSpacing = shift;

  if ( !exists $tree->{'children'} || !@{ $tree->{'children'} } )
  {
    &setExtremes( $tree );
    return;
  }
  &firstWalk( $tree->{'children'}->[ 0 ], $neighborSpacing );

  # Create siblings in contour minimal vertical coordinate and index list.
  my $ih = updateIYL( bottom( $tree->{'children'}->[ 0 ]->{'el'} ), 0, undef );

  for ( my $i = 1 ; $i <= $#{ $tree->{'children'} } ; $i++ )
  {
    &firstWalk( $tree->{'children'}->[ $i ], $neighborSpacing );

    # Store lowest vertical coordinate while extreme nodes still point 
    # in current subtree.
    my $minY = bottom( $tree->{'children'}->[ $i ]->{'er'} );
    &separate( $tree, $i, $ih, $neighborSpacing );
    $ih = updateIYL( $minY, $i, $ih );
  }
  &positionRoot( $tree );
  &setExtremes( $tree );
}

sub setExtremes
{
  my $tree = shift;
  if ( !exists $tree->{'children'} || !@{ $tree->{'children'} } )
  {
    $tree->{'el'}   = $tree;
    $tree->{'er'}   = $tree;
    $tree->{'msel'} = $tree->{'mser'} = 0;
  } else
  {
    $tree->{'el'}   = $tree->{'children'}->[ 0 ]->{'el'};
    $tree->{'msel'} = $tree->{'children'}->[ 0 ]->{'msel'};
    $tree->{'er'} = $tree->{'children'}->[ $#{ $tree->{'children'} } ]->{'er'};
    $tree->{'mser'} =
        $tree->{'children'}->[ $#{ $tree->{'children'} } ]->{'mser'};
  }
}

sub separate
{
  my $tree            = shift;
  my $i               = shift;
  my $ih              = shift;
  my $neighborSpacing = shift;

  # Right contour node of left siblings and its sum of modfiers.
  my $sr = undef;
  $sr = $tree->{'children'}->[ $i - 1 ] if ( @{ $tree->{'children'} } );
  my $mssr = $sr->{'mod'};

  # Left contour node of current subtree and its sum of modfiers.
  my $cl = undef;
  $cl = $tree->{'children'}->[ $i ] if ( @{ $tree->{'children'} } );
  my $mscl = $cl->{'mod'};

  while ( $sr && $cl )
  {
    $ih = $ih->{'next'} if ( &bottom( $sr ) > $ih->{'lowY'} );

    # How far to the left of the right side of sr is the left side of cl?
    my $dist =
        ( $mssr + $sr->{'prelim'} + $sr->{'w'} ) - ( $mscl + $cl->{'prelim'} );

    # RMH: Set additional neighbor spacing
    $dist += $neighborSpacing;

    if ( $dist > 0 )
    {
      $mscl += $dist;
      &moveSubtree( $tree, $i, $ih->{'index'}, $dist );
    }

    my $sy = &bottom( $sr );
    my $cy = &bottom( $cl );

    # Advance highest node(s) and sum(s) of modifiers.
    if ( $sy <= $cy )
    {
      $sr = &nextRightContour( $sr );
      $mssr += $sr->{'mod'} if ( $sr );
    }
    if ( $sy >= $cy )
    {
      $cl = &nextLeftContour( $cl );
      $mscl += $cl->{'mod'} if ( $cl );
    }
  }

  # Set threads and update extreme nodes.
  # In the first case, the current subtree must be taller than 
  # the left siblings.
  if ( !defined $sr && $cl )
  {
    &setLeftThread( $tree, $i, $cl, $mscl );
    # In this case, the left siblings must be taller than the current subtree.
  } elsif ( $sr && !defined $cl )
  {
    &setRightThread( $tree, $i, $sr, $mssr );
  }
}

sub moveSubtree
{
  my $tree = shift;
  my $i    = shift;
  my $si   = shift;
  my $dist = shift;

  # Move subtree by changing mod.
  $tree->{'children'}->[ $i ]->{'mod'}  += $dist;
  $tree->{'children'}->[ $i ]->{'msel'} += $dist;
  $tree->{'children'}->[ $i ]->{'mser'} += $dist;
  &distributeExtra( $tree, $i, $si, $dist );
}

sub nextLeftContour
{
  my $tree = shift;
  if ( !defined $tree->{'children'} || !@{ $tree->{'children'} } )
  {
    return $tree->{'tl'};
  } else
  {
    return $tree->{'children'}->[ 0 ];
  }
}

sub nextRightContour
{
  my $tree = shift;
  if ( !defined $tree->{'children'} || !@{ $tree->{'children'} } )
  {
    return $tree->{'tr'};
  } else
  {
    return $tree->{'children'}->[ $#{ $tree->{'children'} } ];
  }
}

sub bottom
{
  my $tree = shift;
  return ( $tree->{'y'} + $tree->{'h'} );
}

sub setLeftThread
{
  my $tree     = shift;
  my $i        = shift;
  my $cl       = shift;
  my $modsumcl = shift;

  my $li = $tree->{'children'}->[ 0 ]->{'el'};
  $li->{'tl'} = $cl;

  # Change mod so that the sum of modifier after following thread is correct.
  my $diff =
      ( $modsumcl - $cl->{'mod'} ) - $tree->{'children'}->[ 0 ]->{'msel'};
  $li->{'mod'} += $diff;

  # Change preliminary x coordinate so that the node does not move.
  $li->{'prelim'} -= $diff;

  # Update extreme node and its sum of modifiers.
  $tree->{'children'}->[ 0 ]->{'el'}   = $tree->{'children'}->[ $i ]->{'el'};
  $tree->{'children'}->[ 0 ]->{'msel'} = $tree->{'children'}->[ $i ]->{'msel'};
}

# Symmetrical to setLeftThread.
sub setRightThread
{
  my $tree     = shift;
  my $i        = shift;
  my $sr       = shift;
  my $modsumsr = shift;

  my $ri = $tree->{'children'}->[ $i ]->{'er'};
  $ri->{'tr'} = $sr;

  # Change mod so that the sum of modifier after following thread is correct.
  my $diff =
      ( $modsumsr - $sr->{'mod'} ) - $tree->{'children'}->[ $i ]->{'mser'};
  $ri->{'mod'} += $diff;

  # Change preliminary x coordinate so that the node does not move.
  $ri->{'prelim'} -= $diff;

  # Update extreme node and its sum of modifiers.
  $tree->{'children'}->[ $i ]->{'er'} = $tree->{'children'}->[ $i - 1 ]->{'er'};
  $tree->{'children'}->[ $i ]->{'mser'} =
      $tree->{'children'}->[ $i - 1 ]->{'mser'};
}

sub positionRoot
{
  my $tree = shift;

  # Position root between children, taking into account their mod.
  $tree->{'prelim'} =
      ( $tree->{'children'}->[ 0 ]->{'prelim'} +
        $tree->{'children'}->[ 0 ]->{'mod'} +
        $tree->{'children'}->[ $#{ $tree->{'children'} } ]->{'mod'} +
        $tree->{'children'}->[ $#{ $tree->{'children'} } ]->{'prelim'} +
        $tree->{'children'}->[ $#{ $tree->{'children'} } ]->{'w'} ) / 2 -
      ( $tree->{'w'} / 2 );
}

sub secondWalk
{
  my $tree   = shift;
  my $modsum = shift;
  my $minX   = shift;

  $modsum += $tree->{'mod'};

  # Set absolute (non-relative) horizontal coordinate.
  $tree->{'x'} = $tree->{'prelim'} + $modsum;
  $minX = $tree->{'x'} if ( !defined $minX || $minX > $tree->{'x'} );
  &addChildSpacing( $tree );
  for ( my $i = 0 ; $i <= $#{ $tree->{'children'} } ; $i++ )
  {
    $minX = &secondWalk( $tree->{'children'}->[ $i ], $modsum, $minX );
  }
  $minX = $tree->{'x'} if ( !defined $minX || $minX > $tree->{'x'} );
  return ( $minX );
}

sub distributeExtra
{
  my $tree = shift;
  my $i    = shift;
  my $si   = shift;
  my $dist = shift;

  # Are there intermediate children?
  if ( $si != $i - 1 )
  {
    my $nr = $i - $si;
    $tree->{'children'}->[ $si + 1 ]->{'shift'} += $dist / $nr;
    $tree->{'children'}->[ $i ]->{'shift'}  -= $dist / $nr;
    $tree->{'children'}->[ $i ]->{'change'} -= $dist - $dist / $nr;
  }
}

# Process change and shift to add intermediate spacing to mod.
sub addChildSpacing
{
  my $tree = shift;

  my $d           = 0;
  my $modsumdelta = 0;
  for ( my $i = 0 ; $i <= $#{ $tree->{'children'} } ; $i++ )
  {
    $d           += $tree->{'children'}->[ $i ]->{'shift'};
    $modsumdelta += $d + $tree->{'children'}->[ $i ]->{'change'};
    $tree->{'children'}->[ $i ]->{'mod'} += $modsumdelta;
  }
}

# A linked list of the indexes of left siblings and their lowest vertical coordinate.
# IYL = { 'lowY' : lowest vertical coordinate
#         'index': Index for lowest vertical coordinate
#         'next' : pointer to next IYL record
#       }
sub updateIYL
{
  my $minY = shift;
  my $i    = shift;
  my $ih   = shift;

  # Remove siblings that are hidden by the new subtree.
  while ( $ih && $minY >= $ih->{'lowY'} )
  {
    $ih = $ih->{'next'};
  }

  # Prepend the new subtree.
  return ( { 'lowY' => $minY, 'index' => $i, 'next' => $ih } );
}

1;
