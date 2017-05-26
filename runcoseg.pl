#!/usr/bin/perl -w
##---------------------------------------------------------------------------##
##  File:
##      @(#) coseg.pl
##  Author:
##      Robert M. Hubley   rhubley@systemsbiology.org
##  Description:
##      This is a driver script for coseg subfamily generator.
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
#     $Log: runcoseg.pl,v $
#     Revision 1.7  2015/11/04 22:24:35  rhubley
#       Getting ready for a new release
#
#     Revision 1.6  2012/11/06 22:19:01  rhubley
#       - Cleanup
#
#     Revision 1.5  2012/11/02 23:02:19  rhubley
#       - Cleaned up and consolidated EM routines.
#       - Moved final minCount check to prune routine.
#
#     Revision 1.4  2011/06/16 15:59:28  rhubley
#       - Fixed extractSubSeqs.pl
#
#     Revision 1.3  2008/09/25 17:10:10  rhubley
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
#     Revision 1.2  2008/08/06 20:00:40  rhubley
#       - renamed some files
#       - Updated -d parameter in runcoseg ... doc changes etc.
#
#     Revision 1.1  2008/07/23 22:00:24  rhubley
#       - Renamed alkes_scaffold to coseg
#
#     Revision 1.5  2008/07/23 21:57:05  rhubley
#       - Cleanup code
#
#     Revision 1.4  2008/07/22 23:04:44  rhubley
#       - Started implementing disabledSites feature
#       - Cleaned up more code / printing statements
#
#     Revision 1.3  2008/07/22 19:08:49  rhubley
#       - Major changes.  Coded andy's tripple p-value calculation.  Merged
#         coseg_singlemut.c into coseg_scaffold.c so now there is only one C
#         program to maintain.
#
#     Revision 1.2  2008/07/01 23:42:18  rhubley
#       - Added Andy's pValue method
#
#     Revision 1.1  2008/06/11 23:08:10  rhubley
#       - ISB Changes:
#
#           - Added script coseg.pl to assist in running the suite of tools
#             on a dataset.
#
#           - Attempted to disable ALU specific functions.  Much to do here.
#             I want to generalize this "ignore msa region" function so that
#             we can use this on many different repeat classes.
#
#           - Coded a 3-mutation extension to coseg statistical analysis.
#
#           - Created command line processing flags
# 
#
###############################################################################
#
# To Do:
#
=head1 NAME

coseg - Run the coseg subfamily generator

=head1 SYNOPSIS

  coseg [-version] [-m #] [-t] [-k] [-d] [-u #]
        -c <file> -s <file> -i <file> 
        
or

  coseg [-version] [-m #] [-t] [-k] [-d] [-u #]
        -filePrefix <prefix>

=head1 DESCRIPTION

The options are:

=over 4

=item -version

Displays the version of the program

=item -m #

Set the minimum number of elements for the creation of a new subfamily.
The default is 50.

=item -u #

Set the minimum distance between diagnostic sites.  The
default is 10.  Minimum is 1, indicating adjacent sites
will be considered.

=item -filePrefix <prefix>

This is a shortcut to specifying the sequences, consensus and
insertion files separately using the -s, -c, and -i flags.  If
the -filePrefix flag is uses then the value supplied will be
used to automatically generate these parameters as <prefix>.seqs,
<prefix>.cons, and <prefix>.ins.

=item -s <file>

A file containing the aligned sequences ( one per line ).

=item -c <file>

A fasta file containing the consensus sequence used in the alignments.

=item -i <file>

A file containing the removed insertions from the aligned sequences.

=item -k

Use Alkes Price et al. pvalue calculation method.  
[ default is to use Andy Siegel's new pvalue method ]

=item -d

Treat lowercase bases in the consensus sequence as disabled sites
for the entire analysis.  This is useful if your sequences contain
a highly variable internal region which can confound identification
of true subfamilies.

=item -t

Use tri-segregating mutations in calculation. [ default doubles only ]

=back

=head1 SEE ALSO

=head1 COPYRIGHT

Copyright 2007-2012 Robert Hubley, Institute for Systems Biology

=head1 AUTHOR

Robert Hubley <rhubley@systemsbiology.org>

=cut

#
# Module Dependence
#
use strict;
use FindBin;
use lib $FindBin::RealBin;
use lib "$FindBin::RealBin/../";
use Getopt::Long;
use Data::Dumper;

#
# Version
#
#  This is a neat trick.  CVS allows you to tag
#  files in a repository ( i.e. cvs tag "2003/12/03" ).
#  If you check out that release into a new
#  directory with "cvs co -r "2003/12/03" it will
#  place this string into the $Name:  $ space below
#  automatically.  This will help us keep track
#  of which release we are using.  If we simply
#  check out the code as "cvs co Program" the
#  $Name:  $ macro will be blank so we should default
#  to what the ID tag for this file contains.
#
my $CVSNameTag = '$Name:  $';
my $CVSIdTag = '$Id: runcoseg.pl,v 1.7 2015/11/04 22:24:35 rhubley Exp $';
my $Version = $CVSNameTag;
$Version = $CVSIdTag if ( $Version eq "" );

#
# Option processing
#  e.g.
#   -t: Single letter binary option
#   -t=s: String parameters
#   -t=i: Number paramters
#
my @getopt_args = (
    '-version', # print out the version and exit
    '-filePrefix=s',
    '-m=i',
    '-s=s',
    '-c=s',
    '-i=s',
    '-u=s',
    '-d',
    '-t',
    '-k',
    '-v=s',
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

if ( ! $options{'filePrefix'} && 
     (  ! $options{'s'} || ! $options{'c'} || ! $options{'i'} ) )
{
  usage();
}

my $cmd = "$FindBin::RealBin/coseg ";
$cmd .= "-m " . $options{'m'} if ( $options{'m'} );
$cmd .= " -u " . $options{'u'} if ( $options{'u'} );
if ( $options{'filePrefix'} )
{
  my $pre = $options{'filePrefix'};
  $cmd .= " -c $pre.cons -s $pre.seqs";
}else
{
  $cmd .= " -c " . $options{'c'} .  " -s " . $options{'s'}; 
}
$cmd .= " -t " if ( $options{'t'} );
$cmd .= " -k " if ( $options{'k'} );
$cmd .= " -d " if ( $options{'d'} );
$cmd .= " -v $options{'v'} " if ( $options{'v'} );
print "Running: $cmd\n";
system( $cmd );

# Evaluate the run
my $log;
if ( $options{'filePrefix'} )
{
  $log = $options{'filePrefix'} . ".seqs.log";
}else
{
  $log = $options{'s'} . ".log";
}
open IN,"<$log" || 
   die   "ERROR: $log file is missing!  Coseg probably faulted " 
       . "in the middle of the run and didn't create an output file.";
my $goodRun = 0;
while (<IN>)
{
  if ( /^(\d+) subfamilies overall \(\d+ leaves in tree\)/ )
  {
    $goodRun++;
  }
  if ( /^Program duration is.*/ )
  {
    $goodRun++;
  }
}
close IN;
if ( $goodRun < 2 )
{
  die "ERROR: Coseg probably faulted in the middle of the run.  The log file " 
     ."is missing data.  Please report this to the developers.\n";
}
$cmd = "$FindBin::RealBin/postprocess.pl -l sub ";
if ( $options{'filePrefix'} )
{
  my $pre = $options{'filePrefix'};
  $cmd .= " -c $pre.cons -s $pre.seqs";
  $cmd .= " -i $pre.ins";
}else 
{
  $cmd .= " -c " . $options{'c'} . " -s " . $options{'s'};
  $cmd .= " -i " . $options{'i'};
}
print "Running: $cmd\n";
system( $cmd );


1;
