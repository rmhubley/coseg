#!/usr/bin/perl -w
##---------------------------------------------------------------------------##
##  File:
##      @(#) extractSubSeqs.pl
##  Author:
##      Robert M. Hubley   rhubley@systemsbiology.org
##  Description:
##      Given a set of coseg files extract the source sequences
##      from the fasta file for a given subfamily.
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
# ChangeLog
#
#     $Log: extractSubSeqs.pl,v $
#     Revision 1.7  2016/11/11 23:46:04  rhubley
#       - Cleaning up perl scripts
#
#     Revision 1.6  2012/11/02 23:02:19  rhubley
#       - Cleaned up and consolidated EM routines.
#       - Moved final minCount check to prune routine.
#
#     Revision 1.5  2011/06/16 15:59:28  rhubley
#       - Fixed extractSubSeqs.pl
#
#     Revision 1.4  2008/09/25 17:10:10  rhubley
#     - Improved code documentation
#     - Single mutation significance cutoff ( SIGMATHRESH ) was
#       pre-calculated for Alkes Alu analysis and hardcoded.  This
#       version calculates the correct sigma cutoff using the length
#       of the input sequence.
#     - Switched default pvalue method to Andy Siegel's method and
#       provided a new "-k" switch to use Alkes Price's method.
#     - Fixed bug where the program was exiting when calculations
#       fell below the precision of the machine ( epsilon ). Message
#       given was "Below epsilon..." and the runcoseg.pl script
#       moved on even though coseg failed.
#     - Begun to code CpG adjusted consensus routine
#
#     Revision 1.2  2008/07/24 22:17:58  rhubley
#       - Minor changes
# 
#
###############################################################################
#
# To Do:
#
=head1 NAME

extractSubSeqs.pl - Extract fasta elements for a given subfamily.

=head1 SYNOPSIS

  extractSubSeq.pl [-version] -seqFile <coseg seq file> -subFam # 
                              [-fastaDB <fasta file>] 
                              [-flankBases #]

=head1 DESCRIPTION

The options are:

=over 4

=item -seqFile <coseg seq file>

The name of the input sequence file to the coseg program.  NOTE: This is
not the output file created by coseg.  

=item -subFam #

The subfamily number for which you would like sequences extracted.

=item -fastaDB <fasta file>

The fasta file used to create the alignment file fed to coseg

=item -flankBases #

The number of flanking bases to capture from the fastaDB when
extracting the sequences. Default 100 -- only when fastaDB option
is given.

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
use Getopt::Long;
use Data::Dumper;
use FindBin;
use lib $FindBin::RealBin;
use CosegConfig;
use lib $CosegConfig::REPEATMASKER_DIR;
if ( ! -s "$CosegConfig::REPEATMASKER_DIR/FastaDB.pm" )
{
  print "\nError: Could not locate the required module \"FastaDB.pm\" from\n";
  print   "       the RepeatMasker directory ( $CosegConfig::REPEATMASKER_DIR ).\n";
  print   "       Perhaps this path is old or incorrect. Please check the\n";
  print   "       configuration file CosegConfig.pm before continuing.\n\n";
  exit;
}
use FastaDB;


#
# Magic numbers/constants here
#  ie. my $PI = 3.14159;
#   
my $DEBUG = 0;
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
    '-subFam=n',
    '-flankBases=n',
    '-seqFile=s',
    '-fastaDB=s'
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

#
# ARGV Processing
#
if ( ! defined $options{'subFam'} || 
     ! defined $options{'seqFile'}  ) {
  usage();
}


my $prefix = $options{'seqFile'};
my $subnum = $options{'subFam'};

my $seqDB = undef;
if ( defined $options{'fastaDB'} && -f $options{'fastaDB'} )
{
  $seqDB = FastaDB->new( fileName => $options{'fastaDB'},
                         openMode => SeqDBI::ReadOnly );
}

my $numFlankingBases = 100;
if ( $seqDB && $options{'flankBases'} )
{
  $numFlankingBases = $options{'flankBases'};
}

open ASSIGN, "<$prefix.assign" or
   die "Error: Could not open assignment file $options{'seqFile'}.assign\n";

# Create hash of sequences assigned to the subfamily
my %seqs = ();
while ( <ASSIGN> )
{
  if ( /^(\d+)\s+(\d+)/ ) 
  {
    my $seqNum = $1;
    my $subFamNum = $2;
    if ( $subFamNum == $subnum )
    {
      $seqs{ $seqNum } = 1;
    }
  }
}
close ASSIGN;


$prefix =~ s/\.seqs/.fasta/g;
open IN,"<$prefix" or 
  die "Error: Could not open fasta file $prefix\n";
my $seqIndex = -1;
my $header = "";
my $seq = "";
my $cosegID = "";
my $fastaDBID="";
my $fastaDBStart=0;
my $fastaDBEnd=0;
while (<IN>)
{
  #>Seq4 chr1_428953_431249:913-1804 alignment=3
  if ( /^>Seq(\d+)\s+(\S+):(\d+)-(\d+)\s+alignment=\d+/ )
  {
    $cosegID = $1;
    $fastaDBID = $2;
    $fastaDBStart = $3;
    $fastaDBEnd = $4;
    # Is this sequence in the requested subfamily?
    if ( $seqs{$seqIndex} )
    { 
      # Should we give the sequence from the .fasta file
      # ( only as long as the coseg sequence ) or should
      # we give sequence from the fastaDB file ( with 
      # flanking sequences attached )?
      if ( $seqDB )
      {
        # sequence from fastaDB
        my $seqLen = $seqDB->getSeqLength( $fastaDBID );
        # $fastaDBStart is 1-based
        my $start = $fastaDBStart - $numFlankingBases - 1;
        $start = 0 if ( $start < 0 );
        my $end = $fastaDBEnd + $numFlankingBases;
        $end = $seqLen if ( $end > $seqLen );
        $seq = $seqDB->getSubstr( $fastaDBID, $start, $end-$start+1);
        print ">$fastaDBID:$start-$end Seq$cosegID\n$seq\n";
      }else 
      {
        # .fasta sequences only
        $seq =~ s/[\s\n\r]//g;
        print "$header\n$seq\n";
      }
    }
    $seqIndex++;
    $header = $_;
    $header =~ s/[\n\r]//g;
    $seq = "";
  }
  else {
    $seq .= $_;
  }
}
close IN;
# Last entry
if ( $seqs{$seqIndex} )
{ 
  if ( $seqDB )
  {
    # sequence from fastaDB
  }else 
  {
    # .fasta sequences only
    $seq =~ s/[\s\n\r]//g;
    print "$header\n$seq\n";
  }
}


1;
