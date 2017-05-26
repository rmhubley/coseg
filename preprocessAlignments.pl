#!/usr/local/bin/perl -w
##---------------------------------------------------------------------------##
##  File:
##      @(#) preprocessAlignments.pl
##  Author:
##      Robert M. Hubley   rhubley@systemsbiology.org
##  Description:
##     Read alignment output and convert to coseg aligned 
##     seq/ins/cons format.
##
#******************************************************************************
#*  This software is provided ``AS IS'' and any express or implied            *
#*  warranties, including, but not limited to, the implied warranties of      *
#*  merchantability and fitness for a particular purpose, are disclaimed.     *
#*  In no event shall the authors or the Institute for Systems Biology        *
#*  liable for any direct, indirect, incidental, special, exemplary, or       *
#*  consequential damages (including, but not limited to, procurement of      *
#*  substitute goods or services; loss of use, data, or profits; or           *
#*  business interruption) however caused and on any theory of liability,     *
#*  whether in contract, strict liability, or tort (including negligence      *
#*  or otherwise) arising in any way out of the use of this software, even    *
#*  if advised of the possibility of such damage.                             *
#*                                                                            *
#******************************************************************************
#
## ChangeLog
#
##     $Log: preprocessAlignments.pl,v $
##     Revision 1.6  2012/11/02 23:02:19  rhubley
##       - Cleaned up and consolidated EM routines.
##       - Moved final minCount check to prune routine.
##
##     Revision 1.5  2008/09/25 17:10:10  rhubley
##     - Improved code documentation
##     - Single mutation significance cutoff ( SIGMATHRESH ) was
##       pre-calculated for Alkes Alu analysis and hardcoded.  This
##       version calculates the correct sigma cutoff using the length
##       of the input sequence.
##     - Switched default pvalue method to Andy Siegel's method and
##       provided a new "-k" switch to use Alkes Price's method.
##     - Fixed bug where the program was exiting when calculations
##       fell below the precision of the machine ( epsilon ). Message
##       given was "Below epsilon..." and the runcoseg.pl script
##       moved on even though coseg failed.
##     - Begun to code CpG adjusted consensus routine
##
##     Revision 1.4  2008/08/07 18:35:07  rhubley
##      - More doc updates
##
##     Revision 1.3  2008/08/06 20:47:55  rhubley
##       - More docs
##
##     Revision 1.2  2008/07/24 22:17:58  rhubley
##       - Minor changes
##
##     Revision 1.1  2008/07/24 21:33:36  rhubley
##       - Improvements for preprocess script
## 
#
################################################################################
#
# To Do:
#
=head1 NAME

preprocessAlignments.pl - Convert alignments to coseg's aligned seq/ins format.

=head1 SYNOPSIS

  preprocessAlignments.pl [-version] 
                         [-w]
                         [-maxEdgeGap #]
                         [-minConsRange #]
                         [-maxConsRange #]
                         -consensus <consensusFile>
                         -alignments <alignmentFile>

=head1 DESCRIPTION

The options are:

=over 4

=item -w

Specify that the alignment data is in WUBlast format.  The
default is to assume the data is in cross_match format.

=item -maxEdgeGap #

The tolerance for including alignments that do not reach the end 
( as defined by minConsRange and maxConsRange ). Default is 5bp.

=item -minConsRange #

The start of the core analysis range of the alignment. Default is
the start of the consensus.

=item -maxConsRange #

The end of the core analysis range of the alignment. The default is
the end of the consensus.

=item -alignments <alignmentFile>

The alignment file in crossmatch format.

=item -consensus <consensusFile>

The consensus sequence used to create the alignments in FASTA format.

=item -version

Displays the version of the program

=back

=head1 SEE ALSO

=head1 COPYRIGHT

Copyright 2008 Robert Hubley, Institute for Systems Biology

=head1 AUTHOR

Robert Hubley <rhubley@systemsbiology.org>

=cut

#
# Module Dependence
#
use strict;
use FindBin;
use lib $FindBin::RealBin;
## NOTE: This must be set to point to your local RepeatMasker directory.
use lib "/usr/local/RepeatMasker";
use Getopt::Long;
use CrossmatchSearchEngine;
use WUBlastSearchEngine;
use Data::Dumper;

my $Version = 0.1;

#
# Option processing
#  e.g.
#   -t: Single letter binary option
#   -t=s: String parameters
#   -t=i: Number paramters
#
my @getopt_args = (
    '-version', # print out the version and exit
    '-maxEdgeGap=i',
    '-minConsRange=i',
    '-maxConsRange=i',
    '-consensus=s',
    '-alignments=s',
    '-w'
);

my %options = ();
Getopt::Long::config("noignorecase", "bundling_override");
unless (GetOptions(\%options, @getopt_args)) {
    usage();
}

sub usage {
  print "$0 - $Version\n";
  exec "pod2text $0";
  exit;
}

if ($options{'version'}) {
  print "$Version\n";
  exit;
}

if ( ! defined $options{'alignments'} )
{
  print "Missing -alignments parameter!\n";
  usage();
}

my $maxEdgeGap = 5;
if ( defined $options{'maxEdgeGap'} )
{
  $maxEdgeGap = $options{'maxEdgeGap'};
}
   
# The range is unknown at this time.  Need to open the results file  
# to determine the size.
my $maxConsRange = 1;

##
##  Read in the result set
##
my $resultsFile = $options{'alignments'};

my $RMResults;

if ( $options{'w'} )
{
  $RMResults = WUBlastSearchEngine::parseOutput( 
                                       searchOutput => $resultsFile );
}else {
  $RMResults = CrossmatchSearchEngine::parseOutput( 
                                       searchOutput => $resultsFile );
}

##
## Post process the results
##
if ( $RMResults->size() ) 
{
  my $conSize = $RMResults->get(0)->getSubjEnd() + 
                $RMResults->get(0)->getSubjRemaining();
  my $conName = $RMResults->get(0)->getSubjName();

  # Determine the data set's actual min/max cons positions
  my $actualMinStart = $conSize;
  my $actualMaxEnd = 0;
  for ( my $k = 0 ; $k < $RMResults->size() ; $k++ ) 
  {
    my $start = $RMResults->get( $k )->getSubjStart();
    my $end = $RMResults->get( $k )->getSubjEnd();
    $actualMinStart = $start if ( $start < $actualMinStart );
    $actualMaxEnd = $end if ( $end > $actualMaxEnd );
  }
  print "Consensus range: 1 - $conSize\n";
  print "Aligned data consensus range: $actualMinStart - $actualMaxEnd\n";

  # min/max cons range
  my $maxConsRange = $actualMaxEnd;
  $maxConsRange = $options{'maxConsRange'} if ( defined $options{'maxConsRange'} );
  if ( $actualMaxEnd < $maxConsRange )
  { 
    warn "WARNING: maxConsRange specified is $maxConsRange and the largest cons position\n"
       . "in the data set is $actualMaxEnd.  Setting maxConsRange to $actualMaxEnd\n";
    $maxConsRange = $actualMaxEnd;
  }

  my $minConsRange = $actualMinStart;
  $minConsRange = $options{'minConsRange'} if ( defined $options{'minConsRange'} );
  if ( $actualMinStart > $minConsRange )
  { 
    warn "WARNING: minConsRange specified is $minConsRange and the smallest cons position\n"
       . "in the data set is $actualMinStart.  Setting minConsRange to $actualMinStart\n";
    $minConsRange = $actualMinStart;
  }
  print "Using consensus range: $minConsRange - $maxConsRange\n";


  #
  # Open output files
  #
  open ALIGN, ">$resultsFile.seqs";
  open INSERT, ">$resultsFile.ins";
  open FASTA, ">$resultsFile.fasta";
  open OUTL, ">$resultsFile.outliers";
  my $index = 0;


  #my $overallConsensus = "."x(($maxConsRange - $minConsRange)+1);
  my $overallConsensus = "."x($conSize);

  my @coverage = ();
  my $numFiltered = 0;

  for ( my $k = 0 ; $k < $RMResults->size() ; $k++ ) 
  {
    #print "Considering sequence $k\n";

    #
    # Query = sequence
    # Subject = consensus
    #
    my $result = $RMResults->get( $k );

    #
    # Double check that subject is the consensus
    #
    if ( $result->getSubjName() ne $conName ) 
    {
      die "ERROR: The subject may not be the consensus sequence!  Two\n"
         ."different id values identified ( $result-getSubjName() and\n"
         ."$conName ).\n";
    }

    my $conStart = $result->getSubjStart();
    my $conEnd = $result->getSubjEnd();

    for( my $l = $conStart; $l <= $conEnd; $l++ )
    {
      $coverage[$l]++;
    }
  
    #
    # First check that we cover the range
    #
    if ( $conStart <= $minConsRange + $maxEdgeGap &&
         $conEnd >= ( $maxConsRange - $maxEdgeGap ) )
    {
      $index++;

      my $querySeq = $result->getQueryString();
      my $subjSeq = $result->getSubjString();

      #
      # Align with subject in forward direction
      #
      if ( $result->getOrientation() eq "C" ) 
      {
        $querySeq = reverse( $querySeq );
        $querySeq =~ tr/ACGTRYWSKMNXBDHV/TGCAYRSWMKNXVHDB/;
        $subjSeq = reverse( $subjSeq );
        $subjSeq =~ tr/ACGTRYWSKMNXBDHV/TGCAYRSWMKNXVHDB/;
      }  

      my $outQuery = $querySeq;
      $outQuery =~ s/-//g;
      print FASTA ">Seq$index " . $result->getQueryName() . ":" 
                  . $result->getQueryStart() . "-" 
                  . $result->getQueryEnd() . " alignment=$k\n$outQuery\n";
      $outQuery = "";
      
      # 
      # Crudely build a consensus sequence
      #
      my $tmpConsSeq = $subjSeq;
      $tmpConsSeq =~ s/-//g;
      substr( $overallConsensus, $conStart-1, 
              length( $tmpConsSeq ) ) = $tmpConsSeq;

      #
      # Fill in missing edges 
      #
      if ( $conStart > $minConsRange )
      {
        # Add gaps to begining
        $querySeq = "-"x($conStart-$minConsRange) . $querySeq;
        $subjSeq = "N"x($conStart-$minConsRange) . $subjSeq;
        $conStart = $minConsRange;
      }
      if ( $conEnd < $maxConsRange )
      {
        # Add gaps to end 
        $querySeq .= "-"x($maxConsRange-$conEnd);
        $subjSeq .= "N"x($maxConsRange-$conEnd);
        $conEnd = $maxConsRange;
      }

      #
      # Determine the seq position of the min/max consensus
      # base.
      #
      my $minConsPos = 0;
      my $tmpCons = $subjSeq;
      my $baseCounter = $conStart - 1;
      while ( $baseCounter < $minConsRange &&
              $minConsPos <= length( $tmpCons ) )
      {
        if ( substr( $tmpCons, $minConsPos, 1 ) ne "-" )
        {
          $baseCounter++;
        }
        $minConsPos++;
      }
      $minConsPos--;
  

      my $maxConsPos = length( $tmpCons );
      $baseCounter = $conEnd;
      if ( $baseCounter > $maxConsRange ) 
      { 
        do
        {
          $maxConsPos--;
          if ( substr( $tmpCons, $maxConsPos, 1 ) ne "-" )
          {
            $baseCounter--;
          }
        }while ( $baseCounter >= $maxConsRange &&
                 $maxConsPos >= 0 )
      }
   
      #print "Pos = $minConsPos, $maxConsPos\n";
      
      $querySeq = substr( $querySeq, $minConsPos, 
                          ( $maxConsPos - $minConsPos + 1 ) );
      $subjSeq = substr( $subjSeq, $minConsPos, 
                         ( $maxConsPos - $minConsPos + 1 ) );
      $conStart = $minConsRange;
      $conEnd = $maxConsRange;

      # Locate all inserts and replace with single "+" in query seq
      $querySeq =~ s/X/A/ig;
      $querySeq =~ s/N/A/ig;
      my $subjEnd = $conEnd;
     
      my $outSeq = $querySeq;
      my $refSeq = reverse( $subjSeq );
      my @inserts = ();
      my $totalInsertLength = 0;
      #print "RefSeq=$refSeq\n";
      while ( $refSeq =~ /[^-](-+)(?=[^-])/ig )
      {
        #print "Found insert: $1\n";
        $totalInsertLength += length( $1 );
        unshift @inserts,
           ($subjEnd - pos( $refSeq ) + $totalInsertLength ) - $minConsRange .
           ":" .
           substr( $querySeq, length( $querySeq ) - pos( $refSeq ),
                   length( $1 ) );
        substr( $outSeq, length( $querySeq ) - pos( $refSeq ),
                length( $1 ) ) = "+";
      }

      ## Sanity check
      my $numInserts = ($outSeq =~ tr/+/+/ );
      if ( ( length( $outSeq ) - $numInserts ) != 
           ( $maxConsRange - $minConsRange + 1 ) )
      {
        die "Hmm..something went wrong with the sequence trimming " . 
            "final length ( excl + ): " . (length( $outSeq ) - $numInserts) . 
            "\nExpected: " . ( $maxConsRange - $minConsRange + 1 ) . "\n" .
            "query = " . $result->getQueryString() . "\n" .
            "subj  = " . $result->getSubjString() . "\n" .
            "final = $outSeq\n";
      }

      print ALIGN "$outSeq\n";
      #print join( " ", @inserts ) . "\n";
      print INSERT join( " ", @inserts ) . "\n";
    }else {
      $numFiltered++;
      print OUTL "" . $result->toStringFormatted( SearchResult::AlignWithQuerySeq ) . "\n";
    }
  }
  if ( $options{'consensus'} )
  {
    open IN, "<$options{'consensus'}" || 
           die "Error: Could not open $options{'consensus'} file!\n";
    my $seq = "";
    while (<IN>)
    {
      next if ( /^>/ );
      $seq .= $_; 
    }
    close IN;
    $seq =~ s/[\s\r\n]//g;
    if ( length( $seq ) != $conSize )
    {
      die "Error: provided consensus sequence size ( " . length( $seq ) 
           . ")\ndoes not agree with the alignment file ( $conSize )\n";
    }
    $seq = substr( $seq, $minConsRange-1, ( $maxConsRange-$minConsRange+1 ));
    open CONS, ">$resultsFile.cons";
    print CONS ">AlignmentConsensus - from file $options{'consensus'}:"
              ."$minConsRange-$maxConsRange\n" . uc($seq) ."\n";
    close CONS;
  }else {
    my $seq = substr( $overallConsensus, $minConsRange-1, ( $maxConsRange-$minConsRange+1 ));
    if ( $seq =~ /[BDHVRYKMSW]/i )
    {
      print "  Warning: The consensus file being generated ($resultsFile.cons) \n"
           ."           contains IUB codes.  These sites will be treated as\n"
           ."           \"n\" by coseg.  You may want to edit the file and\n"
           ."           change these to a known base before runing coseg.\n";
    }
    open CONS, ">$resultsFile.cons";
    print CONS ">AlignmentConsensus - from file $resultsFile:"
             . "$minConsRange-$maxConsRange\n" . uc($seq) . "\n";
    close CONS;
  }
  print "Total alignments = " . $RMResults->size() . "\n";
  print "Alignments filtered out = $numFiltered\n";
  print "Remaining = " . ( $RMResults->size() - $numFiltered ) . "\n";

  # Testing: Histogram of coverage...
  #open HIST, ">$resultsFile-hist";
  #print HIST join(",", @coverage);
  #close HIST;

  close ALIGN;
  close INSERT;
  close OUTL;
}

1;
