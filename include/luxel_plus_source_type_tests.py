# Written by: John A. Boyd 3
# email: john.boydiii@landauer-fr.com
# Original Publish Date: 24 February, 2025
# Last modification: 10 October 2025

# Luxel+ Algorithm for Whole Body Dosimeters
# Revision 04601

# NOTE: This script is built upon the ROOT analysis framework.
# Certain variables and libraries may be ROOT-specific.
# ROOT Version: 6.32/02, built for win32

from .luxel_plus_error_handling import check_for_negative_elements
from .luxel_plus_boundary_and_fit_params import *
import math

def _testIfWithinEllipsoid(x, y, z, x0, y0, z0, a, b, c_ellipsoid):
	withinEllipsoid = False

	ellipsoid = (pow( (x - x0), 2))/pow(a, 2) + (pow( (y - y0), 2))/pow(b, 2) + (pow( (z - z0), 2))/pow(c_ellipsoid, 2)

	if ellipsoid <= 1:
		withinEllipsoid = True

	return withinEllipsoid

def pure_BL_test(OW, PL, Al, Cu, standardErrorConditions = True):
	# Function to determine if a source is a Pure BL Beta (Kr85).
	# These fall within a unique region far beyond all other sources and so they are relatively "easy" to pick off
	# The ratio of OW to any other element (OW/PL, OW/Al, OW/Cu) is huge and so we can simply call anything with 
	# a large ratio of, say, OW/PL, a Pure BL.
	# To prevent any "potentially unforeseen" contaminatrion, we also apply a bound on Cu/OW to just provide an upper limit

	pure_BL_beta = False

	if not standardErrorConditions:
		OW, PL, Al, Cu = check_for_negative_elements(OW, PL, Al, Cu)

	OWPL = OW/PL 
	CuOW = Cu/OW

	BL_Beta_OWPL_Limit = 4.0 
	CuOW_Beta_Upper_Limit = 0.125 

	if( (OWPL > BL_Beta_OWPL_Limit) and (CuOW < CuOW_Beta_Upper_Limit) ):
		# BL Beta
		pure_BL_beta = True 

	return pure_BL_beta

def pure_BH_test(OW, PL, Al, Cu, standardErrorConditions = True):
	# Function to check if a given source is a Pure High-Energy Beta Source (Sr90).  

	pure_BH_beta = False 

	if not standardErrorConditions:
		OW, PL, Al, Cu = check_for_negative_elements(OW, PL, Al, Cu)

	# /////
	# From a 3D plot of (OW/Cu) [x-axis], (OW/PL) [y-axis], and (OW/Al) [z-axis] we can define an elliptical
	# boundary around the Pure BH region.
	# This region is defined by 3 semi-major axes and 3 center locations.
	# We can then check if a given combination of (OW, PL, Al, Cu) lies within the elliptical boundary. 

	# x: OW/Cu
	# y: OW/PL
	# z: OW/Al

	x0 = 19.829787
	a = 9.5083431

	y0 = 1.7589584
	b = 0.26849493

	z0 = 2.6139127
	c = 0.49276144

	x = OW/Cu
	y = OW/PL
	z = OW/Al

	owpl = OW/PL
	owal = OW/Al
	owcu = OW/Cu

	# We can do a simple upper bound cut to grab any of the PureBH that are not near the boundary of the MixedBH sources
	# The max (OW/PL) value for the MixedBH sources is 1.8177637 (conversative cut above this: 1.85).
	# the max (OW/PL) value for the PureBH sources is 2.1437745 (conservative upper cut limit for this: 2.5);
	# Just to be sure about not going too far "above" the PureBH in "z" (OW/AL), the max (OW/AL) is: 3.3635650;
	# We know that PureBL are way above this but we can go ahead and cap at 4.0

	# We want to make sure that we can still allow some identifcation of PureBH beyond our ellipsoid
	# Particularly in the +y and +/-z axex directions beyond the regions of the MixedBH
	# So, letting the (OW/PL) boundary for MixedBH end at 1.85 and arbitrarily capping it at 2.5, along with an (OW/Al) cut at 4.0:
	# We also need to make sure that we don't push this "volume" into the Mixed_BL volume----->
	# ----> Mixed_BL: We apply a Cu/OW (1.0/x) cut above 0.15 to prevent picking out Mixed_BL

	# We have that:
	if ((1.0/x) < 0.15) and (1.85 < y) and (y < 2.5) and (z < 4.0):
		pure_BH_beta = True
	else:
		isWithinTestVolume = _testIfWithinEllipsoid(x, y, z, x0, y0, z0, a, b, c)

		if isWithinTestVolume:
			pure_BH_beta = True

	return pure_BH_beta

def pure_NS20_test(OW, PL, Al, Cu, standardErrorConditions = True):
	# Function to check if a given source is a Pure NS20 Beta Source. 

	pure_NS20 = False 

	if not standardErrorConditions:
		OW, PL, Al, Cu = check_for_negative_elements(OW, PL, Al, Cu)

	owcu = OW/Cu 
	owpl = OW/PL 
	owal = OW/Al 

	# /////
	# From a 3D plot of (OW/Cu) [x-axis], (OW/PL) [y-axis], and (OW/Al) [z-axis] we can define an elliptical
	# boundary around the Pure BH region.
	# This region is defined by 3 semi-major axes and 3 center locations.
	# We can then check if a given combination of (OW, PL, Al, Cu) lies within the elliptical boundary. 

	# x: OW/Cu
	# y: OW/PL
	# z: OW/Al	

	## ------------

	# Updated 9 September 2025
	# AFter further analysis, we can perform a simple cut on OW/Cu to determine Pure NS20 sources.
	# We only lose (mis-identify) a very small portion of Pure NS20 and MixedBH sources

	# Number of MixedBH entries with (OW/Cu) > 59.5 ----> 27 (out of 477909 total entries);
	# Number of PureNS20 entries with (OW/Cu) < 59.5 ----> 4 (out of 10000 total entries);

	# We also do a cut on OW/AL as well to make sure we don't accidentally select any Pure M30
	# (OW/Al) limits on PureNS20: (2.7707636 <= OW/Al <= 3.9755065)
	# (OW/Al) limits on PureM30: (1.8497425 <= OW/Al <=  2.4073312)

	# So, we can set a boundary at the mid-point between M30 and NS20: (2.4073312 + 2.7707636)/2.0 = 2.5890474
	# We go ahead and pick some upper bound at 5.0 to be safe

	NS20_owcu_limit = 59.5

	NS20_owal_lower_limit = 2.5890474 
	NS20_owal_upper_limit = 5.0 

	if( owcu >= NS20_owcu_limit ):
		if (NS20_owal_lower_limit <= owal) and (owal <= NS20_owal_upper_limit):
			pure_NS20 = True 

	return pure_NS20

def pure_M30_test(OW, PL, Al, Cu, standardErrorConditions = True):
	# Function to check if a given source is a Pure M30 Beta Source. 

	pure_M30 = False

	if not standardErrorConditions:
		OW, PL, Al, Cu = check_for_negative_elements(OW, PL, Al, Cu)

	owcu = OW/Cu
	owpl = OW/PL 
	owal = OW/Al 

# 	/////
# 	From a 3D plot of (OW/Cu) [x-axis], (OW/PL) [y-axis], and (OW/Al) [z-axis] we can define an elliptical
# 	boundary around the Pure BH region.
# 	This region is defined by 3 semi-major axes and 3 center locations.
# 	We can then check if a given combination of (OW, PL, Al, Cu) lies within the elliptical boundary. 

# 	x: OW/Cu
# 	y: OW/PL
# 	z: OW/Al

# 	Updated 10 September 2025
# 	We can do a combination of tests which involve an ellipsoid to allow to get close to the MixedBH boundary
# 	But also, we can allow a wider region beyond that (in +x (OW/Cu) and +/- z (OW/Al) )

# 	M30 wide, flat, general regions "far" beyond MixedBH but also considering other Pure Photon sources. 

	if (owcu > 30.5) and (owpl < 1.05) and ( (1.83 < owal) and (owal < 2.41) ):
		pure_M30 = True 
	else:
		# M30 volume
		x0 = 48.408094
		a = 25.591335

		y0 = 1.0582187
		b = 0.11519215

		z0 = 2.1337529
		c = 0.28750167

		x = OW/Cu
		y = OW/PL
		z = OW/Al

		isWithinTestVolume = _testIfWithinEllipsoid(x, y, z, x0, y0, z0, a, b, c)

		if isWithinTestVolume:
			pure_M30 = True 

	return pure_M30

def MixedBH_M30_test(OW, PL, Al, Cu, standardErrorConditions = True):
	# Function to check if a given source is a MixedBH M30:Sr90 Source. 

	is_MixedBH_M30 = False 

	# lower_line_slope = -1.75
	# lower_line_intercept = 5.1 
	# lower_line = lower_line_slope*(PL/Al) + lower_line_intercept 
	lower_line = get_OWAl_yVal_for_boundary_line_between_MixedBHM30_and_MixedBHHigherEnergies_from_PLAl(PL, Al)

	# upper_line_slope = -0.75
	# upper_line_intercept = 3.975 
	# upper_line = upper_line_slope*(PL/Al) + upper_line_intercept
	upper_line = get_OWAl_yVal_for_boundary_line_between_MixedBHNS20_and_MixedBHM30_from_PLAl(PL, Al)

	# Middle region is where we are greater than the lower line AND less than the upper line.
	# Therefore, we are in the lower-right overlapping region of the two lines. 
	# Here we only consider values in the middle region

	if ((OW/Al) >= lower_line) and ((OW/Al) <= upper_line):
		is_pure_BH = pure_BH_test(OW, PL, Al, Cu, standardErrorConditions)
		is_pure_M30 = pure_M30_test(OW, PL, Al, Cu, standardErrorConditions)

		if is_pure_BH or is_pure_M30:
			is_MixedBH_M30 = False
		else:
			is_MixedBH_M30 = True 
	else:
		is_MixedBH_M30 = False

	return is_MixedBH_M30

def MixedBH_NS20_test(OW, PL, Al, Cu, standardErrorConditions = True):
	# Function to check if a given source is a MixedBH NS20:Sr90 Source.

	is_MixedBH_NS20 = False 

	# lower_line_slope = -1.75
	# lower_line_intercept = 5.1
	lower_line = get_OWAl_yVal_for_boundary_line_between_MixedBHM30_and_MixedBHHigherEnergies_from_PLAl(PL, Al)

	# upper_line_slope = -0.75
	# upper_line_intercept = 3.975
	upper_line = get_OWAl_yVal_for_boundary_line_between_MixedBHNS20_and_MixedBHM30_from_PLAl(PL, Al)

	# middle region is where we are greater than the lower line AND less than the upper line.
	# Therefore, we are in the lower-right overlapping region of the two lines. 
	# Here we only consider values in the middle region

	if ((OW/Al) >= lower_line) and ((OW/Al) <= upper_line):
		is_pure_BH = pure_BH_test(OW, PL, Al, Cu, standardErrorConditions)
		is_pure_NS20 = pure_NS20_test(OW, PL, Al, Cu, standardErrorConditions)

		if is_pure_BH or is_pure_NS20:
			is_MixedBH_NS20 = False 
		else:
			is_MixedBH_NS20 = True 
	else:
		is_MixedBH_NS20 = False 

	return is_MixedBH_NS20

def pure_beta_test(OW, PL, Al, Cu, return_DU_test = False, b_rad_env_class_DU = False, standardErrorConditions = True):
	# check if a given source is a Pure Beta Source. 
	# Also contains a special sub-test for determining if it is a Pure Depleted Uranium source
	# Returns the result of the Pure DU test if return_DU_test = true. 	

	pure_beta_test_output = False 

	#### ----------
	if not standardErrorConditions:
		OW, PL, Al, Cu = check_for_negative_elements(OW, PL, Al, Cu)

	#### ----------

	pure_beta = False 
	pure_DU = False 

	OWPL = OW/PL 
	CuOW = Cu/OW 

	DU_OWPL_Lower_Limit = None

	BL_Beta_OWPL_Limit = 4.0
	BH_Beta_Upper_OWPL_Limit = 4.0
	BH_Beta_CuOW_Lower_Limit = 0.001
	BH_Beta_CuOW_Upper_Limit = 0.1
	BH_OWPL_Lower_Limit = 1.44

	CuOW_Beta_Upper_Limit = 0.125

    # We can adjust the limits of the DU limit if we are dealing with a UBeta Rad ENV Class.
    # We know that we have to have a DU to pass as a pure beta, so, we can "open" up or "ease" the limit a bit
    # From our test dataset, I averaged the values of OW/PL for UBeta for OW/PL < 2.0 and the average came to 1.94

    # Therefore, we can set the value to that. 

	if b_rad_env_class_DU:
		# Special DU_OWPL limit for UBeta Rad ENV Class
		DU_OWPL_Lower_Limit = 1.94
	else:
		# Standard Condition
		DU_OWPL_Lower_Limit = 2.15

	# We need a special consideration for the Uranium Beta Environment Class
	# For this case, we only report that we have a beta IFF that beta is a Pure DU.
	# If it is not a pure DU then we move on.
	# We don't allow it be to a BL or BH. If NOT DU, then not Pure Beta. 

	if b_rad_env_class_DU:
		if (DU_OWPL_Lower_Limit < OWPL) and (OWPL < BL_Beta_OWPL_Limit) and (CuOW < CuOW_Beta_Upper_Limit):
			pure_DU = True

		if pure_DU is False:
			pure_beta_test_output = False
		else:
			pure_beta_test_output = True 
	else:
		if pure_BL_test(OW, PL, Al, Cu, standardErrorConditions):
			# BL Beta
			pure_beta = True 
		else:
			# BH Beta

			pure_beta = pure_BH_test(OW, PL, Al, Cu, standardErrorConditions)

		if return_DU_test:
			if (DU_OWPL_Lower_Limit < OWPL) and (OWPL < BL_Beta_OWPL_Limit) and (CuOW < CuOW_Beta_Upper_Limit):
				pure_DU = True

	if return_DU_test:
		pure_beta_test_output = pure_DU 
	else:
		pure_beta_test_output = pure_beta 

	return pure_beta_test_output

def pure_photon_test(OW, PL, Al, Cu, standardErrorConditions = True):
	# Function to check if a given source is a Pure Photon Source
	# -------

	if not standardErrorConditions:
		OW, PL, Al, Cu = check_for_negative_elements(OW, PL, Al, Cu)

	pure_photon_test_output = False 
	pure_photon_test = False 

	# Here we draw a boundary to separate Pure Photons from other sources
	# In order to be as selective as possible without erroneously grabbing too many MixedBH we can use a composite line/boundary. 

	# Considering a plot of (OW/PL) vs (OW/Cu) we can draw two parabolas (opening upwards), a negatively-sloped line and a positively-sloped line

	# Define region separation/boundary lines

	# Boundary between the two parabolas:

	parabola_owcu_intersection = 8.50

	# Right boundary of the right parabola (to the right of this boundary is the straight line).
	right_parab_owcu_boundary = 21.8

	# Parameters for the left-most parabola: 
	left_parab_c0 = 0.0094191
	left_parab_c1 = -0.098501
	left_parab_c2 = 1.33

	# Parameters for the right parabola: 
	right_parab_c0 = 0.0013201
	right_parab_c1 = -0.041246
	right_parab_c2 = 1.4282

	# Negatively-sloped line parameter
	# starts at x = 21.8 (y = 1.1564)
	# ends at x = 30.7476 (y = 1.091618)
	# ----------------
	negative_linear_boundary_right_endpoint = 30.7476
	negative_linear_boundary_slope = -0.0072
	negative_linear_boundary_intercept = 1.3142

	# Positively-sloped line parameters: 
	# Previous: 
	linear_boundary_slope = 0.0016
	linear_boundary_intercept = 1.1207

	# starts at x = 30.7476 (y = 1.091618)
	positive_linear_boundary_slope = 0.0033
	positive_linear_boundary_intercept = 0.9907


	OWPL = OW/PL 
	OWCu = OW/Cu 

	owpl_boundary_value = 0.0 

	# Test (boolean) values for Pure NS20 and Pure M30:
	is_pure_NS20 = False 
	is_pure_M30 = False 

	# Special tests for the two lowest energy Pure Photon Sources: NS20 and M30
	is_pure_NS20 = pure_NS20_test(OW, PL, Al, Cu, standardErrorConditions)
	is_pure_M30 = pure_M30_test(OW, PL, Al, Cu, standardErrorConditions)

    # If we look at a plot of (OW/Al) vs (PL/Al) we can draw boundary lines to identify NS20 and M30

    # Boundary line between separating (NS20 & M30) from all other (higher-energy) pure photons:
	# M30_NS20_boundary_line_slope = -1.75
	# M30_NS20_boundary_line_intercept = 5.1
	# M30_NS20_boundary_line = M30_NS20_boundary_line_slope*(PL/Al) + M30_NS20_boundary_line_intercept
	M30_NS20_boundary_line = get_OWAl_yVal_for_boundary_line_between_MixedBHM30_and_MixedBHHigherEnergies_from_PLAl(PL, Al)

	# We can go ahead and say that if we pass the Pure NS20 and Pure M30 tests we have a Pure Photon:
	if is_pure_NS20 or is_pure_M30:
		pure_photon_test = True 
	else:
		if( OW/Al >= M30_NS20_boundary_line):
			# We've already separated Pure NS20 and Pure M30 but we can check for other Photon sources
			# Here we need to still allow for other Pure Photons (or Photon:Photon Mixtures) in this region to be detected
			# y = OW/Al, x = PL/Al

			pure_photon_plal_owal_line_slope = 1.235
			pure_photon_plal_owal_line_intercept = -0.20
			pure_photon_plal_owal_line = pure_photon_plal_owal_line_slope*(PL/Al) + pure_photon_plal_owal_line_intercept

			if( OW/Al < pure_photon_plal_owal_line ):
				pure_photon_test = True 
			else:
				pure_photon_test = False
		else:
			# All other (non NS20 & M30) pure photons:

			if( OWCu < parabola_owcu_intersection ):
				# Left parabola (highest energy photons)
				owpl_boundary_value = left_parab_c0*math.pow(OWCu, 2) + left_parab_c1*OWCu + left_parab_c2 
			elif (parabola_owcu_intersection <= OWCu) and (OWCu < right_parab_owcu_boundary):
				# Right Parabola
				owpl_boundary_value = right_parab_c0*math.pow(OWCu, 2) + right_parab_c1*OWCu + right_parab_c2 
			# elif (right_parab_owcu_boundary <= OWCu) and (OWCu < negative_linear_boundary_right_endpoint):
			# 	# Negatively-sloped line to the right of the right parabola:
			# 	owpl_boundary_value = negative_linear_boundary_slope*OWCu + negative_linear_boundary_intercept
			# else:
			# 	# Postively-sloped line to the right of the negatively sloped line
			# 	owpl_boundary_value = positive_linear_boundary_slope*OWCu + positive_linear_boundary_intercept
			else:
			# Previous
				# Single positively-sloped linear boundary:
				owpl_boundary_value = linear_boundary_slope*OWCu + linear_boundary_intercept

			if( OWPL <= owpl_boundary_value ):
				pure_photon_test = True 
			else:
				pure_photon_test = False 

	return pure_photon_test

def mixed_source_test( beta_test, photon_test ):
	# Function to check if a given source is a Mixed Source.
	# ---------------

	mixed_source_output = False 

	# If it isn't a Pure Beta and not a Pure Photon then it is a Mixed BetaPhoton
	# If it is both a Pure Beta AND a Pure Photon, this is an error, and it is a Mixed BetaPhoton
	if (beta_test is False and photon_test is False) or (beta_test and photon_test):
		mixed_source_output = True 

	return mixed_source_output	