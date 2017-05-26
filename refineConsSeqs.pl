#!/usr/local/bin/perl -w
##---------------------------------------------------------------------------##
##  File:
##      @(#) refineConsSeqs
##  Author:
##      Robert M. Hubley   rhubley@systemsbiology.org
##  Description:
##      A script to refine subfamily consensi derived by coseg
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
#     $Log: refineConsSeqs.pl,v $
#     Revision 1.1  2015/11/04 22:24:35  rhubley
#       Getting ready for a new release
# 
#
###############################################################################
#
# To Do:
#
=head1 NAME

refineConsSeqs - Refine subfamily consensi derived by coseg

=head1 SYNOPSIS

  refineConsSeqs [-version] -subConsFile <*.subfamilies.seq>

=head1 DESCRIPTION

  Create rootCons.fa : The subfamily with the highest average Kimura divergence
         refined_subs.fa : The refined consensus for each of the coseg subfamilies

=over 4

=item -version

Displays the version of the program

=back

=head1 SEE ALSO

=head1 COPYRIGHT

Copyright 2012 Robert Hubley, Institute for Systems Biology

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
use Cwd;

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
my $CVSIdTag = '$Id: refineConsSeqs.pl,v 1.1 2015/11/04 22:24:35 rhubley Exp $';
my $Version = $CVSNameTag;
$Version = $CVSIdTag if ( $Version eq "" );

##----------------------------------------------------------------------##
##      S I T E   S P E C I F I C   C O N F I G U R A T I O N
##
##  If you must include site specific variables in the program
##  itself put them here.
##
##  ie. my $blastPrgrmDir = "/user/local/blast/bin";
##      my $indelPenalty = 30;
##
##  END OF SITE SPECIFIC CONFIGURATION -- DO NOT EDIT BELOW THIS LINE
##----------------------------------------------------------------------##

#
# Magic numbers/constants here
#  ie. my $PI = 3.14159;
#   
my $DEBUG = 0;
my $RepeatModelerDir = "/usr/local/RepeatModeler";

#
# Option processing
#  e.g.
#   -t: Single letter binary option
#   -t=s: String parameters
#   -t=i: Number paramters
#
my @getopt_args = (
    '-version', # print out the version and exit
    '-subConsFile=s',
);

my %options = ();
Getopt::Long::config("noignorecase", "bundling_override");
unless (GetOptions(\%options, @getopt_args)) {
    usage();
}

sub usage {
  print "$0 - $Version\n\n";
  exec "pod2text $0";
  exit;
}

if ($options{'version'}) {
  print "$Version\n";
  exit;
}

if ( ! $options{'subConsFile'} )
{
  print "\n\nMissing option subConsFile!\n\n";
  usage();
}

my $date = localtime( time() );
$date =~ s/[ ,\t,\n:]//g;
my $scratchDir =  cwd() . "/REF_$$.$date";
mkdir( $scratchDir );

##
## Gather subfamilies defined by Bi/Tri-segregating mutations
##
my $conFile = $options{'subConsFile'};
my $seqFile;
if ( $conFile =~ /(.*)\.subfamilies.seq/ )
{
   $seqFile = $1;
}else
{
  die "Could not locate *.seqs file corresponding to the *.subfamilies.seq file: $seqFile!\n"; 
}

open IN,"<$conFile" or 
 die "Could not open $conFile: $!\n";
my @sigSubs = ( 0 );
while ( <IN> )
{
  # >subfamily1 count=100 pvalue=3e-133 ( parent_pvalue=3e-40 ) parent=0 muts= 152:t/- 142:c/t 80:a/g 
  if ( /^>subfamily(\d+).*muts=\s*(\S.*)/ )
  {
    my $subNum = $1;
    my $muts = $2;
    my $numMuts = ($muts =~ tr/:/:/);
    push @sigSubs, $subNum if ( $numMuts > 1 );
  }
}
close IN;

##
##  Create Fasta Files For Each Subfamily
##
my %consSeqs = ();
my %consDiv = ();
foreach my $sn ( @sigSubs )
{
  system("$FindBin::Bin/extractSubSeqs.pl -seqFile $seqFile -subFam $sn > $scratchDir/sub$sn.fa");
  system("cd $scratchDir; $RepeatModelerDir/Refiner sub$sn.fa");
  open CONSIN,"<$scratchDir/sub$sn.fa.refiner_cons";
  open CONSOUT,">>$scratchDir/refined_subs.fa";
  
  my $seq = "";
  while ( <CONSIN> )
  {
    #>family ( Final Multiple Alignment Size = 153 , Avg Kimura = 15.8 )
    s/family/sub$sn/g;
    print CONSOUT $_;
    if ( /Avg Kimura = ([\d\.]+)/ )
    {
      $consDiv{ "sub$sn" } = $1;
      next;
    }
    s/[\s\n\r]//g;
    $seq .= $_;
  }
  $consSeqs{ "sub$sn" } = $seq;
  close CONSIN;
  close CONSOUT;
}
system("cp $scratchDir/RM_*/*.html $scratchDir");

my $highestDivSub = (sort { $consDiv{ $b } <=> $consDiv{ $a } } keys( %consDiv ))[0];

open OUT,">$scratchDir/rootCons.fa" or die "Could not open rootCons.fa for writing!\n";
print OUT ">rootSub $highestDivSub\n";
print OUT "$consSeqs{ $highestDivSub }\n";
close OUT;

system("cp $scratchDir/rootCons.fa ./");
system("cp $scratchDir/refined_subs.fa ./");

unlink("$scratchDir/*.nhr");
unlink("$scratchDir/*.nin");
unlink("$scratchDir/*.nsq");

print "Done!\n\n";
print "A scratch directory was created ( $scratchDir ) that would\n";
print "typically be removed, however it remains in this release for\n";
print "debugging purposes.\n\n";
print "Output Files: ./rootCons.fa\n";
print "              ./refined_subs.fa\n";
print "\n\n";

