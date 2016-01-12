#!/usr/bin/perl -w

use strict;
use lib ("$ENV{GEMC}/api/perl");
use utils;
use materials;

# Help Message
sub help()
{
	print "\n Usage: \n";
	print "   materials.pl <configuration filename>\n";
 	print "   Will create the electronics readout board PCB materials\n";
 	print "   Note: The passport and .visa files must be present to connect to MYSQL. \n\n";
	exit;
}

# Make sure the argument list is correct
# If not pring the help
if( scalar @ARGV != 1)
{
	help();
	exit;
}

# Loading configuration file and paramters
our %configuration = load_configuration($ARGV[0]);

# One can change the "variation" here if one is desired different from the config.dat
# $configuration{"variation"} = "myvar";

sub print_materials
{
	# uploading the mat definition
	
	# BONUS Sensitive region He (user define the density)
	my %mat = init_mat();
	$mat{"name"}          = "GasHe";
	$mat{"description"}   = "bonus Sensitive He gas region";
	$mat{"density"}       = "0.0001786";
	$mat{"ncomponents"}   = "1";
	$mat{"components"}    = "He4_1atm 1";
	print_mat(\%configuration, \%mat);

	
	
}

print_materials();

