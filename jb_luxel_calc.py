# Written by: John A. Boyd 3
# email: john.boydiii@landauer-fr.com
# Original Publish Date: 24 February, 2025
# Last modification: 12 August 2025

# Luxel+ Algorithm for Whole Body Dosimeters

# NOTE: This script is built upon the ROOT analysis framework.
# Certain variables and libraries may be ROOT-specific.
# ROOT Version: 6.24/08, built for win32

# This script is intended to be used with as wrapper script.
# Said wrapper script should access the following functions, feed them element values,
# and receive, as returns, the desired variables. 
# At its simplest, one can simply provide 4 element values for OW, PL, Al, and CU
# to any given function and receive the calculated result. 
# --------------
# The two primary calculations are calc_total_SDE and calc_total_DDE.
# These two functions will calculate the SDE and DDE, respectively, of any given 4 inputs.
# They will access all other functions, as necessary, in determining the proper downstream
# handling and routing of the element values for processing and final calculation. 
# --------------
# Other commonly reported variables are: radiation_quality, report_radiation_quality,
# beta_indicator, and calc_continuous_energy. 

import numpy as np
import math
from include.luxel_plus_boundary_and_fit_params import *
from include.luxel_plus_dose_calc_coefficients import *
from include.luxel_plus_energy_functions import *
from include.luxel_plus_error_handling import *
from include.luxel_plus_radiation_quality_functions import *
from include.luxel_plus_reporting_functions import *
from include.luxel_plus_source_type_indicators import *
from include.luxel_plus_source_type_tests import *

def calc_photon_total_SDE(OW, PL, Al, Cu, known_photon_rq_input = None, standardErrorConditions = True):
	# Function to calculate the SDE for Pure Photons
	# -----------------------------

	if not standardErrorConditions:
		check_for_negative_elements(OW, PL, Al, Cu)

	photon_total_SDE = 0.0

	if known_photon_rq_input is None:
		photon_rq = photon_rad_quality(OW, PL, Al, Cu)
	else:
		photon_rq = known_photon_rq_input

	if photon_rq == "PH":

		# 10 Sep 2025
		# ***Use these coefficients***
		c1 = 4.23273e-12
		c2 = 4.8378e-13
		c3 = 7.05436e-12
		c4 = 0.982287
		#Refined mean squared error: 8.38826e+07 (9158.75)

	elif photon_rq == "PL" or photon_rq == "PM":

		# 10 Sep 2025
		# ***Use these coefficients***
		c1 = 0.983565
		c2 = -0.677245
		c3 = -0.131354
		c4 = 0.76856
		#Refined mean squared error: 1.70627e+08 (13062.4)

		# # 17 Oct 2025
		# c1 = 1.38627
		# c2 = -1.08848
		# c3 = -0.101479
		# c4 = 0.750755
		# #Refined mean squared error: 3.29266e+08 (18145.7)
	else:
		# Mixed BH
		v_coeffs = get_Mixed_BH_SDE_coefficients(OW, PL, Al, Cu)
		c1 = v_coeffs[0]
		c2 = v_coeffs[1]
		c3 = v_coeffs[2]
		c4 = v_coeffs[3]

	photon_total_SDE = (c1*OW) + (c2*PL) + (c3*Al) + (c4*Cu)

	return photon_total_SDE

def calc_beta_total_SDE(OW, PL, Al, Cu, known_beta_rq_input = None, v_b_rad_ENV_classes = [False, False, False], standardErrorConditions = True):
	# Function to calculate the SDE for Pure, NON-DU, Beta Sources
	# DU has its own SDE function, however, this function will get used in case of mis-identification. 

	# If beta_rq_input is empty, then the Beta RQ will first be evaluated.
	# Otherwise, the value of Beta RQ is taken from beta_rq_input.
	# -----------------------

	if not standardErrorConditions:
		check_for_negative_elements(OW, PL, Al, Cu)

	beta_total_SDE = 0.0	

	b_rad_env_class_DU = v_b_rad_ENV_classes[1]

	if known_beta_rq_input is None:
		# beta_rq_input possibilities: BL, BH
		beta_rq = beta_rad_quality(OW, PL, Al, Cu, b_rad_env_class_DU)
	else:
		beta_rq = known_beta_rq_input

	if beta_rq == "BH":
		# 14 Oct 2025
		# ***Use these coefficients***
		c1 = 0.277212
		c2 = 0.736014
		c3 = 0.485456
		c4 = 0.411977
		#/Refined mean squared error: 66896.9 (258.644)

	elif beta_rq == "BL":
		# 10 Sep 2025
		# ***Use these coefficients***
		c1 = 2.20836
		c2 = 0.000999875
		c3 = 0.000999875
		c4 = 0.000999953
		#Refined mean squared error: 91085.7 (301.804)
	else:
		# Mixed BH
		v_coeffs = get_Mixed_BH_SDE_coefficients(OW, PL, Al, Cu)
		c1 = v_coeffs[0]
		c2 = v_coeffs[1]
		c3 = v_coeffs[2]
		c4 = v_coeffs[3]

	beta_total_SDE = (c1*OW) + (c2*PL) + (c3*Al) + (c4*Cu)

	return beta_total_SDE

def calc_DU_total_SDE(OW, PL, Al, Cu, standardErrorConditions = True):
	# Function to calculate the SDE for Pure DU sources
	# It is possible this function may get used for Pure BH that are mis-identified as DU
	# -----------------------

	if not standardErrorConditions:
		check_for_negative_elements(OW, PL, Al, Cu)

	DU_SDE = 0.0

	# 10 Sep 2025
	# ***Use these coefficients***
	c1 = 1.43564
	c2 = -0.83034
	c3 = -3.3665
	c4 = 27.0087
	#Refined mean squared error: 1.0906e+06 (1044.32)

	DU_SDE = (c1*OW) + (c2*PL) + (c3*Al) + (c4*Cu)

	return DU_SDE	

def calc_mixed_total_SDE(OW, PL, Al, Cu, known_mixed_rq_input = None, known_source = None, standardErrorConditions = True):
	# Function to calculate the SDE of Mixed Sources
	# -----------------------

	if not standardErrorConditions:
		check_for_negative_elements(OW, PL, Al, Cu)

	mixed_SDE = 0.0

	if known_mixed_rq_input is None:
		# beta_rq_input_possibilities: BL, BH
		mixed_RQ = mixed_BetaPhotons_radiation_quality(OW, PL, Al, Cu)
	else:
		mixed_RQ = known_mixed_rq_input

	if mixed_RQ == "Mixed_BL":
		v_coeffs = get_Mixed_BL_SDE_coefficients()
		c1 = v_coeffs[0]
		c2 = v_coeffs[1]
		c3 = v_coeffs[2]
		c4 = v_coeffs[3]	

	else:
		v_coeffs = get_Mixed_BH_SDE_coefficients(OW, PL, Al, Cu, known_source)
		c1 = v_coeffs[0]
		c2 = v_coeffs[1]
		c3 = v_coeffs[2]
		c4 = v_coeffs[3]
	

	mixed_SDE = (c1*OW) + (c2*PL) + (c3*Al) + (c4*Cu)

	return mixed_SDE

def calc_total_SDE(OW, PL, Al, Cu, branch_rq, standardErrorConditions = True):
	# Function to calculate the SDE for a given set of Element Values
	# The function will determine/evaluate source type and RQ and send the element values to the appropriate SDE Fn.
	# -----------------------

	# Vectors to hold ENV Class and PGC Zone booleans
	# 	ENV Class: {Photon only, UBeta, Xray Only}
	# 	PGC Zones: {Photon Only, Col & Waist, Extremity}

	if not standardErrorConditions:
		check_for_negative_elements(OW, PL, Al, Cu)

	total_SDE = 0.0

	if branch_rq == "BL":
		total_SDE = calc_beta_total_SDE(OW, PL, Al, Cu, "BL")
	elif branch_rq == "BH":
		total_SDE = calc_beta_total_SDE(OW, PL, Al, Cu, "BH")
	elif branch_rq == "BH_DU":
		total_SDE = calc_DU_total_SDE(OW, PL, Al, Cu)
	elif branch_rq == "PL" or branch_rq == "PM" or branch_rq == "PH":
		total_SDE = calc_photon_total_SDE(OW, PL, Al, Cu)
	else:
		if branch_rq == "Mixed_BL":
			total_SDE = calc_mixed_total_SDE(OW, PL, Al, Cu, "Mixed_BL")
		else:
			is_MixedBH_M30 = MixedBH_M30_test(OW, PL, Al, Cu, standardErrorConditions)
			is_MixedBH_NS20 = MixedBH_NS20_test(OW, PL, Al, Cu, standardErrorConditions)

			if is_MixedBH_M30:
				total_SDE = calc_mixed_total_SDE(OW, PL, Al, Cu, "Mixed_BH", "M30")
			elif is_MixedBH_NS20:
				total_SDE = calc_mixed_total_SDE(OW, PL, Al, Cu, "Mixed_BH", "NS20")
			else:
				total_SDE = calc_mixed_total_SDE(OW, PL, Al, Cu, "Mixed_BH")

	return total_SDE

def calc_total_DDE(OW, PL, Al, Cu, branch_rq, standardErrorConditions = True):
	# Function for calculating the DDE for a given set of element values
	# The function will determine/evaluate source type and proceed as necessary
	# -----------------------

	# Vectors to hold ENV Class and PGC Zone booleans
	# 	ENV Class: {Photon only, UBeta, Xray Only}
	# 	PGC Zones: {Photon Only, Col & Waist, Extremity}

	if not standardErrorConditions:
		check_for_negative_elements(OW, PL, Al, Cu)

	total_DDE = 0.0
	source_type = 0

	CuOW_Limit = 0.75

	if branch_rq == "BL" or branch_rq == "BH" or branch_rq == "BH_DU":
		source_type = 1

	# Source type = 1 --> Pure Beta
	if source_type == 1:
		total_DDE = Cu 
	else:
		# NOT Pure Beta --> Pure Photon or Mixed

		# Mixed
		# Updated Post the May 2025 failure of DOELAP & NVLAP tests.
		# We now add a branch to separate the Mixed BL:Photons and Mixed BH:Photons	

		if branch_rq.startswith("Mixed_B"):
			if branch_rq == "Mixed_BL":
				v_coeffs = get_Mixed_BL_DDE_coefficients()
				c1 = v_coeffs[0]
				c2 = v_coeffs[1]
				c3 = v_coeffs[2]
				c4 = v_coeffs[3]
			else:
				# Mixed_BH
				# vector<double> v_coeffs = get_MixedBH_UpperEnergies_DDE_coefficients(OW, PL, Al, Cu);
				# 14 Oct:				
				v_coeffs = get_Mixed_BH_DDE_coefficients(OW, PL, Al, Cu)
				c1 = v_coeffs[0]
				c2 = v_coeffs[1]
				c3 = v_coeffs[2]
				c4 = v_coeffs[3]
		else:
			# Pure Photons
			photon_rq = photon_rad_quality(OW, PL, Al, Cu)

			# Higher Energy
			if photon_rq == "PH":
				# 10 Sep 2025
				# ***Use these coefficients***
				c1 = 0.00213403
				c2 = 0.0
				c3 = 0.0
				c4 = 1.02008
				#Refined mean squared error: 4.97376e+07 (7052.49)
			elif photon_rq == "PL" or photon_rq == "PM":
				# 10 Sep 2025
				# ***Use these coefficients***
				c1 = 1.34875e-11
				c2 = 0.0487997
				c3 = 0.150351
				c4 = 0.600787
				#Refined mean squared error: 3.23585e+07 (5688.45)

				# # 17 Oct 2025
				# c1 = 1.20551e-11
				# c2 = 0.0631591
				# c3 = 0.0906002
				# c4 = 0.724383
				# #Refined mean squared error: 6.35652e+07 (7972.78)
			else:
				v_coeffs = get_Mixed_BH_DDE_coefficients(OW, PL, Al, Cu)
				c1 = v_coeffs[0]
				c2 = v_coeffs[1]
				c3 = v_coeffs[2]
				c4 = v_coeffs[3]

		total_DDE = c1*OW + c2*PL + c3*Al + c4*Cu

	return total_DDE

def calc_LDE(DDE, SDE):
	LDE_calc = 0.0

	if SDE < 1.0:
		LDE_calc = SDE 
	else: 
		LDE_SdeDde_limit = SDE*(1.4 - (1.04*(np.exp(-DDE/SDE))))

		if LDE_SdeDde_limit > SDE:
			LDE_calc = SDE 
		else: 
			LDE_calc = SDE*(1.4 - (1.04*(np.exp(-DDE/SDE))))

	return LDE_calc

def calc_beta_dose(beta_indicator, Cu, DDE, SDE, Rad_Qual):

	# This function will determine a "Beta DDE" for a given source which contains a Known Beta
	# We should be able to determine a Beta DDE for all Pure Beta Sources and Mixed Beta:Photons
	# If we know the source type (Rad_Qual) already, that is great and will save calculation time
	# But we can determine that in any case.

	# This function will return a pair with the value of beta_dose and a note/remark.

	beta_dose = 0.0
	remark = ""

	# If we have a beta in the source, i.e. beta_indicator = 1, then we should return a beta_dose
	# Otherwise, we will check the type of pure source that it is (either pure photon or pure beta)

	if beta_indicator:
		# Check for SDE = 0:
		# If we have a beta_indicator, then we can't have a zero-value in either SDE or beta_dose
		# If DDE = 0 then we use the value of Cu as DDE

		if DDE == 0:
			DDE = Cu 
		
		if Rad_Qual == "BL" or Rad_Qual == "BH" or Rad_Qual == "BH_DU":
			# For these we can use the Cu element value with some considerations.
			# Some Pure Betas can reach the 10mm depth (DDE) and should give a dose in the range of 10% of the SDE
			# We can use this as an upper limit

			Max_Pure_beta_dose = 0.1*SDE 

			if Cu >= Max_Pure_beta_dose:
				beta_dose = Max_Pure_beta_dose
			else:
				beta_dose = SDE - DDE 
		else:
			# This should only be reached by:
			# Mixed_PL, Mixed_PM, Mixed_PL, or Mixed_BL_PH only

			beta_dose = SDE - DDE 

		if beta_dose <= 1.0:
			beta_dose = 1.0 

			# If we have a value of beta_dose, or Cu, less than 1 we report this as "Minimal" or "M"
			remark = "M"
	else:
		# We don't really need a check here for the Pure Photon Radiation Qualities but.....
		# It is a good safe guard in case we missed some configuration/combination of values and/or inputs		
		if Rad_Qual == "PL" or Rad_Qual == "PM" or Rad_Qual == "PH":
			# In case we are trying to check this for a Pure Photon source.... we rpoert nothing and remark on it.
			beta_dose = 0.0
			remark = Rad_qual + "_Photons_Only"
		else:
			# Not sure we should ever reach this branch as this indicates:
			# No beta_indicator, and not a pure photon.
			# Therefore, in case we missed some combination of RQs and beta_ind, this can account for it. 
			# Set beta_dose to the error value: -1
			beta_dose = -1
			remark = "Beta_Dose_Calc_ERR_" + Rad_Qual

	beta_dose_with_remark = [beta_dose, remark]

	return beta_dose_with_remark

def calc_beta_dose_and_adjustedSDE(beta_indicator, Cu, DDE, SDE, Rad_Qual, b_extremity_zone = False, b_PO_ENV_Class_and_UBeta_Zone = False, v_b_rad_ENV_classes = [False, False, False]):
	
	# The return will be of the form: { {beta_dose, adjustedSDE}, "Remark"}

	# This function will determine a "Beta DDE" for a given source which contains a Known Beta
	# We should be able to determine a Beta DDE for all Pure Beta Sources and Mixed Beta:Photons
	# If we know the source type (Rad_Qual) already, that is great and will save calculation time
	# But we can determine that in any case.

	# This function will return a pair with the value of beta_dose and a note/remark.

	beta_dose = 1.0
	remark = ""
	SDE_beta_dose_adjusted = -1.0

	# If we have an "extremity", i.e. b_extremity_zone = true, then DDE (and LDE) are zero.
	# Therefore, we need to separate these accordingly.
	# If we have a beta_source then we set SDE to Cu, and beta_dose to the calculated value (with SDE = Cu and DDE = 0).
	# Otherwise, if there is no beta present then beta_dose is zero, and SDE is the calculated value.
	
	# If we have a beta in the source, i.e. beta_indicator = 1, then we should return a beta_dose
	# Otherwise, we will check the type of pure source that it is (either pure photon or pure beta)	

	if beta_indicator:
		# Check for SDE = 0:
		# If we have a beta_indicator, then we can't have a zero-value in either SDE or beta_dose
		# If DDE = 0 then we use the value of Cu as DDE

		if DDE == 0:
			DDE = Cu

		if Rad_Qual == "BL" or Rad_Qual == "BH" or Rad_Qual == "BH_DU":
			# For these we can use the Cu element value with some considerations.
			# Some Pure Betas can reach the 10mm depth (DDE) and should give a dose in the range of 10% of the SDE
			# We can use this as an upper limit.

			Max_Pure_beta_dose = 0.1*SDE 

			if Cu >= Max_Pure_beta_dose:
				betadose = Max_Pure_beta_dose
				SDE_beta_dose_adjusted = SDE - beta_dose
			else:
				beta_dose = SDE - DDE 
				SDE_beta_dose_adjusted = DDE 
		else:
			# This should only be reached by:
			# MixedPL, Mixed_PM, Mixed_PL, or Mixed_BL_PH only
			beta_dose = SDE - DDE
			SDE_beta_dose_adjusted = DDE 

		# We set the minimum reported beta_dose here. Anything less than 1 is set to 1
		if beta_dose <= 1.0:
			beta_dose = 1.0 
			SDE_beta_dose_adjusted = SDE - 1.0 
			# If we have a value for beta_dose, or Cu, less than 1, we report this as "Minimal" or "M".
			remark = "M" 
	else:
		# We don't really need a check here for the Pure Photon Radiation Qualities but.....
		# It is a good safe guard in case we missed some configuration/combination of values and/or inputs

		if Rad_Qual == "PL" or Rad_Qual == "PM" or Rad_Qual == "PH":
			# In case we are trying to check this for a Pure Photon Source.... we report nothing and remark on it.
			beta_dose = 0.0 
			SDE_beta_dose_adjusted = SDE 
			remark = Rad_Qual + "_Photons_Only" 
		else: 
			# Here we are left with the "default" option of: Mixed_XX
			# No beta_indicator, and not a pure photon.

			beta_dose = 0 
			SDE_beta_dose_adjusted = SDE 
			remark = Rad_Qual 

	b_rad_env_class_DU = v_b_rad_ENV_classes[1] 

	if (b_PO_ENV_Class_and_UBeta_Zone or b_rad_env_class_DU) and (Rad_Qual == "P"):
		beta_dose = 0.0
		SDE_beta_dose_adjusted = SDE 

	beta_dose_with_remark = [[beta_dose, SDE_beta_dose_adjusted], remark]					

	return beta_dose_with_remark

# def main():
#     # your program logic goes here
#     print("Running jb_luxel_calc...")
#     OW = 155
#     PL = 1
#     Al = 1.5
#     Cu = 2
#     branch_RQ = "BL"
#     SDE = calc_total_SDE(OW, PL, Al, Cu, branch_RQ)

#     print(f"SDE: {SDE}")


# if __name__ == "__main__":
#     main()