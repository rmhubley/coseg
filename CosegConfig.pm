#!/usr/bin/perl
##---------------------------------------------------------------------------##
##  File:
##      @(#) CosegConfig.pm
##  Author:
##      Robert Hubley <rhubley@systemsbiology.org>
##  Description:
##      This is the main configuration file for the Coseg
##      perl programs.  Before you can run the programs included
##      in this package you will need to edit this file and
##      configure for your site. 
##
#******************************************************************************
#* Copyright (C) Institute for Systems Biology 2016 Developed by
#* Arian Smit and Robert Hubley.
#*
#* This work is licensed under the Open Source License v2.1.  To view a copy
#* of this license, visit http://www.opensource.org/licenses/osl-2.1.php or
#* see the license.txt file contained in this distribution.
#*
###############################################################################
package CosegConfig;
use FindBin;
require Exporter;
@EXPORT_OK = qw( $REPEATMASKER_DIR );

%EXPORT_TAGS = ( all => [ @EXPORT_OK ] );
@ISA         = qw(Exporter);

BEGIN {
##----------------------------------------------------------------------##
##     CONFIGURE THE FOLLOWING PARAMETERS FOR YOUR INSTALLATION         ##
##                                                                      ##
##
## RepeatMasker Location
## ======================
## The path to the RepeatMasker programs and support files
## This is the directory with this file as well as
## the ProcessRepeats and Library/ and Matrices/ subdirectories
## reside.
##
##    i.e. Typical UNIX installation
##     $REPEATMASKER_DIR = "/usr/local/RepeatMasker";
##
  $REPEATMASKER_DIR          = "/usr/local/RepeatMasker";

##                                                                      ##
##                      END CONFIGURATION AREA                          ##
##----------------------------------------------------------------------##
}

sub standalone_entry_point
{
  print "Enter location of the RepeatMasker program: ";
  my $answer = <STDIN>;
  $answer =~ s/[\n\r]+//g;
  # TODO Validate
  
  open IN,"<CosegConfig.pm" or die;
  open OUT,">CosegConfig.new" or die;
  while ( <IN> )
  {
    if ( /^\s*\$REPEATMASKER_DIR\s*\=/ )
    {
      print OUT "  \$REPEATMASKER_DIR          = \"$answer\";\n";
    }else
    {
      print OUT;
    }
  }
  close IN;
  close OUT;

  system("mv CosegConfig.new CosegConfig.pm");
  exit;
}

## Allow this module to be called as a standalone script
__PACKAGE__->standalone_entry_point() unless caller;

1;
