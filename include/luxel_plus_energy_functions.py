# Written by: John A. Boyd 3
# email: john.boydiii@landauer-fr.com
# Original Publish Date: 24 February, 2025

# Luxel+ Algorithm for Whole Body Dosimeters
# Revision 04601

# NOTE: This script is built upon the ROOT analysis framework.
# Certain variables and libraries may be ROOT-specific.
# ROOT Version: 6.24/08, built for win32

import numpy as np
import math
from .luxel_plus_radiation_quality_functions import *

def calc_pure_photon_energy(OW, PL, Al, Cu, standardErrorConditions = True):
	# Function to approximate the energy of a Pure Photon Source
	# Based on a fit determined from a plot of Energy vs (OW/Cu)(Al/PL)
	# -----------------------

	if not standardErrorConditions:
		check_for_negative_elements(OW, PL, Al, Cu)

    # energy = c0*math.log((c1*Al_Cu - c2)) + c3
	a = 4.0
	b = -0.7
	c = 3.0
	d = 12.0
	f = -1.3

	elem_factor = (OW/Cu)*(Al/PL)

	energy = 0.0

	energy = math.exp( a + b*(math.log(elem_factor - f)) + c/(elem_factor - f)) + d

	if elem_factor <= 1.055:
		energy = 662
	else:
		if energy > 662:
			energy = 662
		elif energy < 16:
			energy = 16
		else:
			energy = energy 

	return energy 

def calc_pure_beta_energy(OW, PL, Al, Cu, beta_rq_input = None, standardErrorConditions = True):
	# Function to evaluate the energy of a given Pure Beta Source.
	# If no beta_rq_input is provided, this function will evaluate a Beta RQ and proceed.
	# -----------------------

	if not standardErrorConditions:
		check_for_negative_elements(OW, PL, Al, Cu)

	energy = 0.0 
	beta_rq = "" 

	if beta_rq_input is None:
		beta_rq = beta_rad_quality(OW, PL, Al, Cu)
	else:
		beta_rq = beta_rq_input

	if beta_rq == "BH_DU":
    	# Uranium Slab
		energy = 620
	elif beta_rq == "BL":
    	# BL: Kr85
		energy = 687
	else:
		# BH: Sr90
		energy = 2284

	return energy 

def calcGaussianInInterval(x, xmin, xmax, Amplitude):
	# Creates a gaussian that spans the xmin and xmax interval.
	# Given an x-value, the y-value of this gaussian is returned. 
	# The amplitude of the gaussian is specified by the user as "Amplitude"
    # Compute mean (mu) and standard deviation (sigma)
	# -----------------------
	gaus_y_min = 16
	Amplitude = Amplitude - gaus_y_min

	mu = (xmin + xmax) / 2.0
	sigma = (xmax - xmin) / (2.0 * math.sqrt(-2.0 * math.log(0.01))) # Derived from boundary conditions

	# Compute Gaussian value
	gaussian = Amplitude * math.exp(-math.pow(x - mu, 2) / (2 * math.pow(sigma, 2))) + gaus_y_min

	return gaussian

def calcHalfGaussianInInterval(x, xmin, xmax, Amplitude):
	# Creates the Left Hand Side of a Gaussian that would typically span twice the 
	# interval of xmin to xmax.
	# Given an x-value, the y-value of this half-gaussian is then returned.
	# The amplitude of the Gaussian is specified by the user as "Amplitude".
	# The mean (mu) and width (sigma or StdDev) are determined herewithin. 
	# -----------------------
	gaus_y_min = 16
	Amplitude = Amplitude - gaus_y_min

    # Compute mean (mu) at xmax (the peak of the left-half Gaussian)
	mu = xmax

    # Compute standard deviation (sigma) such that the Gaussian is 1% of Amplitude at xmin
	sigma = (xmax - xmin) / math.math.sqrt(-2.0 * math.log(0.01))

    # Compute and return Gaussian value
	gaussian = Amplitude * exp(-math.pow(x - mu, 2) / (2 * math.pow(sigma, 2))) + gaus_y_min

	return gaussian

def calc_continuous_energy(OW, PL, Al, Cu, rq_input = None, standardErrorConditions = True):
	# Function for evaluating the energy of a source as a "continuous" energy function.
	# We use a mix of various functions and combine them together along various ranges. 
	# We will stitch together calc_pure_photon_energy, calcGaussianInInterval,
	# calcHalfGaussianInInterval, and the individual Pure Beta Source energies
	# at given intervals of x {(OW/Cu)(Al/PL)}.
	# -----------------------

	if not standardErrorConditions:
		check_for_negative_elements(OW, PL, Al, Cu)

	energy = 0.0 
	elem_factor = (OW/Cu)*(Al/PL)

	source_type = 0

	if rq_input is not None:
		if rq_input.startswith("BL") or rq_input.startswith("BH"):
			source_type = 1
	else:
		source_type = source_type_indicator(OW, PL, Al, Cu)

	# We have 3 specific ranges where we can try to pick out high energy betas

	# Ranges: Uranium Slab: (OW/Cu)*(Al/PL) = 12.89 +/- 0.157099 -> Energy = 620
	#        Sr90: (OW/Cu)*(Al/PL) = 17.018 +/- 0.5 -> Energy = 2284
	#        Kr85: (OW/Cu)*(Al/PL) > 45 -> Energy = 687

	if source_type == 1:
		beta_rq = beta_rad_quality(OW, PL, Al, Cu)

		if beta_rq == "BH_DU":
			# Uranium Slab
			energy = 620
		elif beta_rq == "BL":
			# BL - Kr85
			energy = 687
		else:
			# BH - Sr90
			energy = 2284
	else:
		# DU energy is centered about 12.89.
		# We define a range about this mean, empirically determined as: +/- 0.157099
		DU_xmin = 12.89 - 0.157099
		DU_xmax = 12.89 + 0.157099

		# Sr90 energy is centered about 17.018
		# The range about this was empirically determined as: +/- 0.5
		Sr90_xmin = 17.018 - 0.5
		Sr90_xmax = 17.018 + 0.5

		if (elem_factor > DU_xmin) and (elem_factor < DU_xmax):
			# Uranium Slab
			DU_energy = 620.0
			energy = calcGaussianInInterval(elem_factor, DU_xmin, DU_xmax, DU_energy)
		elif (elem_factor > Sr90_xmin) and (elem_factor < Sr90_xmax):
			# Sr90
			Sr90_energy = 2284.0
			energy = calcGaussianInInterval(elem_factor, Sr90_xmin, Sr90_xmax, Sr90_energy)
		elif (elem_factor > 45) and (elem_factor < 50):
			# Transition to Kr85
			# This transition is defined as the Left Hand Gaussian "leading up" to the Kr85 Energy
			energy = calcHalfGaussianInInterval(elem_factor, 45, 50, 687)
		elif elem_factor >= 50:
			# Kr85
			Kr85_energy = 687.0
			energy = Kr85_energy
		else:
			# Continuous energy function
			energy = calc_pure_photon_energy(OW, PL, Al, Cu)

		photon_rq = photon_rad_quality(OW, PL, Al, Cu)

		if (photon_rq == "PH") and (energy < 205):
			energy = 662
	
	return energy