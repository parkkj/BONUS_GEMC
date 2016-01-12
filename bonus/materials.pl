#!/usr/bin/perl -w

use strict;
use lib ("$ENV{GEMC}/io");
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
# From Jixie
#    ultem = C37H24O6N2
#    density = 1.27*g/cm3;
#    ultem = new G4Material(name="Ultem", density, nElem=4);
#    ultem->AddElement(elH, nAtoms=37);
#    ultem->AddElement(elC, nAtoms=24);
#    ultem->AddElement(elO, nAtoms=6);
#    ultem->AddElement(elN, nAtoms=2);

sub print_materials
{
	# uploading the mat definition
	
	# BONUS electronics readout board
	my %mat = init_mat();
	$mat{"name"}          = "pcb";
	$mat{"description"}   = "bonus electronics readout board material";
        $mat{"density"}       = "1.27";     
	$mat{"ncomponents"}   = "4";
	$mat{"components"}    = "H 37 C 24 O 6 N 2";
	print_mat(\%configuration, \%mat);

	# PC Board
#	%mat = init_mat();
#	$mat{"name"}          = "pcb";
#	$mat{"description"}   = "BONUS (same as BST) pc board material";
#	$mat{"density"}       = "1.860";
#	$mat{"ncomponents"}   = "3";
#	$mat{"components"}    = "G4_Fe 0.3 G4_C 0.4 G4_Si 0.3";
#	print_mat(\%configuration, \%mat);


	    
	# Rohacell
	%mat = init_mat();
	$mat{"name"}          = "rohacell";
	$mat{"description"}   = "bst rohacell material";
	$mat{"density"}       = "0.11";
	$mat{"ncomponents"}   = "4";
	$mat{"components"}    = "G4_C 0.6465 G4_H 0.0784 G4_N 0.0839 G4_O 0.1912";
	print_mat(\%configuration, \%mat);
	



	# BONUS Sensitive region He (user define the density)
        %mat = init_mat();
	$mat{"name"}          = "GasHe";
	$mat{"description"}   = "bonus buffer region He4 gas 1atm";
	$mat{"density"}       = "0.0001786";
	$mat{"ncomponents"}   = "1";
	$mat{"components"}    = "He4_1atm 1";
	print_mat(\%configuration, \%mat);

	
}

print_materials();

