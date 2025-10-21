# Written by: John A. Boyd 3
# email: john.boydiii@landauer-fr.com
# Original Publish Date: 24 February, 2025

# Luxel+ Algorithm for Whole Body Dosimeters
# Revision 04601

# NOTE: This script is built upon the ROOT analysis framework.
# Certain variables and libraries may be ROOT-specific.
# ROOT Version: 6.24/08, built for win32

# This script is intended to be used as wrapper script.
# This will access functions in jb_luxel_calc.C, feed them the element values,
# and receive, as returns, the desired variables. 
# At its simplest, one can simply provide 4 element values for OW, PL, Al, and CU
# to any given function and receive the calculated result. 
# --------------

from dataclasses import dataclass

import jb_luxel_calc
import math

from jb_luxel_calc import check_error_conditions
from jb_luxel_calc import test_for_negative_elements
from jb_luxel_calc import radiation_quality
from jb_luxel_calc import convert_RQ_string_to_int
from jb_luxel_calc import calc_total_DDE
from jb_luxel_calc import calc_total_SDE
from jb_luxel_calc import calc_LDE
from jb_luxel_calc import assign_beta_indicator
from jb_luxel_calc import report_radiation_quality
from jb_luxel_calc import calc_beta_dose_and_adjustedSDE
from jb_luxel_calc import calc_continuous_energy

@dataclass
class return_list_value_and_descr:
	# Inner pair of the return list
	value: float
	descr: str

@dataclass
class return_list_name_and_data:
	# Outer pair which contains the name of the value, and the value & desc
	value_name: str
	data: return_list_value_and_descr

@dataclass
class luxel_plus_calculations:
	# data class which holds all of the calculations and return values for the algorithm
	Error: return_list_value_and_descr
	DDE: return_list_value_and_descr
	SDE: return_list_value_and_descr
	LDE: return_list_value_and_descr
	Branch_RQ: return_list_value_and_descr
	Reported_RQ: return_list_value_and_descr
	Beta_Indicator: return_list_value_and_descr
	Beta_Dose: return_list_value_and_descr
	Energy: return_list_value_and_descr
	SDE_preadjusted: return_list_value_and_descr
	Known_OW: float
	Known_PL: float
	Known_Al: float
	Known_Cu: float
	Known_DDE: float
	Known_SDE: float
	Sample_ID: int
	Source: str
	Source_Type: str
	Ratio: None

def luxel_plus_algorithm(OW, PL, Al, Cu, PGC_RAD_ENV_Class = -99, PGC_ZONE = None, controlValues = [-1, -1, -1, -1], standardErrorConditions = True, print_detailed_err = False):
	# First component of pair: Descriptor string of the variable, i.e., "SDE", "DDE", "RQ", etc.
	# Second component: 
	# 		Nested Pair containing a: Numerical value and TString as necessary by variable type. 
	# 		In the case of variables like RQ which return a string, the numerical value will indicate a specified string value to be described as necessary.

	# Description of elements in "return_vector:
	# [0] Main error condition. 
	# 		0 NO ERRORS
	# 		1 Errors exist BUT we didn't have to skip the data due to the error (err2 and err4).
	# 		2 Errors exist AND we SHOULD SKIP the data due to the error (err1 and err3). Populates all other array elements to err condition. 
	# [1] DDE
	# [2] SDE
	# [3] LDE
	# [4] BranchRQ
	# 		0 --> BL
	# 		1 --> BH_DU
	# 		2 --> BH
	# 		3 --> PL
	# 		4 --> PM
	# 		5 --> PH
	# 		6 --> Mixed_PL
	# 		7 --> Mixed_PM
	# 		8 --> Mixed_PH
	# [5] ReportedRQ	
	# 		Same as above with an additional possibility:
	# 		9 --> Mixed_BL_PH
	# [6] Beta Indicator
	# [7] Beta Dose
	# [8] Energy

	print_individual_calculations = False 

	# Returns a list of lists
	num_return_values = 10
	def make_default_descr():
		return return_list_value_and_descr(value=None, descr=None)
		
	calculation_results = luxel_plus_calculations(
		Error=make_default_descr(),
		DDE=make_default_descr(),
		SDE=make_default_descr(),
		LDE=make_default_descr(),
		Branch_RQ=make_default_descr(),
		Reported_RQ=make_default_descr(),
		Beta_Indicator=make_default_descr(),
		Beta_Dose=make_default_descr(),
		Energy=make_default_descr(),
		SDE_preadjusted=make_default_descr(),
		Known_OW=OW,
		Known_PL=PL,
		Known_Al=Al,
		Known_Cu=Cu,
		Known_DDE=-1,
		Known_SDE=-1,
		Sample_ID = None,
		Source = None,
		Source_Type = None,
		Ratio = None,
	)

	# 	We will process all possible variables/calculations and store them in the single vector.

	# This portion will assess the "Radiation Environment Class"
	# Rad Env Classes:
	# 	Photon Only: 	3
	# 	Uranium Beta: 	6
	# 	X-Ray Only:		7

	# 	int known_source_type_input = 0;
	# 	Known Source Type Inputs:
	# 	Pure Beta: 			1
	# 	Pure Photon: 		2
	# 	Mixed Beta:Photon:	-1

	# Define variables and set defaults:
	DDE = 0.0
	SDE = 0.0
	SDE_beta_dose_adjusted = 0.0
	LDE = 0.0 
	Beta_Dose = [[0.0, 0.0], ""]
	Beta_Indicator = None 
	Energy = 0.0

	Branch_RQ = None
	Reported_RQ = None 

	b_rad_env_class_PO = False
	b_rad_env_class_DU = False
	b_rad_env_class_XRO = False 
	b_UBeta_ENV_Class_and_PO_Zone = False
	b_PO_ENV_Class = False
	b_UBeta_ENV_Class = False 
	b_PO_Zone = False 

	# List to hold ENV Class and PGC Zone booleans:

	# ENV Class: [Photon only, UBeta, Xray Only]

	v_b_rad_ENV_classes = [False, False, False]

	if PGC_RAD_ENV_Class == 3:
		# Pure Photon Only
		# known_source_type_input = 2
		b_rad_env_class_PO = True 

		b_PO_ENV_Class = True 
	elif PGC_RAD_ENV_Class == 6:
		# Pure Uranium Beta
		# known_source_type_input = 1
		# known_rq_input = "BH_DU"
		b_rad_env_class_DU = True 

		# Beta_Indicator = 1
		b_UBeta_ENV_Class = True 
	elif PGC_RAD_ENV_Class == 7:
		# X-Ray Only
		# Taken from DoseRulesLuxel:
			# Handle x-ray only environment class - assign PL as radiation quality code
		# known_source_type_input = 2
		# known_rq_input = "PL"
		b_rad_env_class_XRO = True

		b_PO_ENV_Class = True 
	else:
		# Really we do nothing here....

		# known_source_type_input = known_source_type_input;
		# known_rq_input = known_rq_input;
		# b_rad_env_class_XRO = false;
		place_holder = False 

	v_b_rad_ENV_classes = [b_rad_env_class_PO, b_rad_env_class_DU, b_rad_env_class_XRO]

	# This portion will assess the "Zone"
	# Zones: 
	# 	Whole Body:			"01"
	# 	Right Wrist: 		"05"
	# 	Left Wrist: 			"06"
	# 	Right Ankle:			"11"
	# 	Left Ankle: 			"12"
	# 	Collar:				"13"
	# 	Waist:				"14"
	# 	Outer Extremity:		"20"	

	b_collar_and_waist_zone = False
	b_extremity_zone = False
	b_collar_waist_extremity_with_UBeta = False

	# Collar and Waist Zones
	if PGC_ZONE == "13" or PGC_ZONE == "14":
		# Collar and Waist are considered as "Photon Only" and have Zero Beta Component.
		# Taken from: enmCalcLuxelAlumVer0313Nvlap
		# 	'Photon Only Data Check
	    #         If vintRadEnvClassID = PGC_RAD_ENV_CLASS_PO Or _
	    #            vintRadEnvClassID = PGC_RAD_ENV_CLASS_XRO Or _
	    #            vintRadEnvClassID = PGC_RAD_ENV_CLASS_H150 Or _
	    #            vintRadEnvClassID = PGC_RAD_ENV_CLASS_AM241 Or _
	    #            vstrZoneCode = PGC_ZONE_CLR Or _
	    #            vstrZoneCode = PGC_ZONE_WST Then		
		
		# Pure Photon Only
		# known_source_type_input = 2;

		b_collar_and_waist_zone = True 

		if PGC_RAD_ENV_Class == 6:
			b_collar_waist_extremity_with_UBeta = True 

		b_PO_Zone = True 
	# Extremities Ankle, wrist, outer:
	elif PGC_ZONE == "05" or PGC_ZONE == "06" or PGC_ZONE == "11" or PGC_ZONE == "12" or PGC_ZONE == "20":
		# We calculate the extremities in the same way, EXCEPT
		# There are no DDE or LDE calculations provided for these zones.

		b_extremity_zone = True 

		if PGC_RAD_ENV_Class == 6:
			b_collar_waist_extremity_with_UBeta = True 
	else:
		# The only thing left is Whole Body: "01"
		b_extremity_zone = False 
		b_collar_and_waist_zone = False 

	# We need to account for the special case that we have a Pure Photon Env (or extremity) AND a UBeta Zone
	if (b_collar_and_waist_zone or b_extremity_zone or b_PO_Zone) and b_UBeta_ENV_Class:
		b_UBeta_ENV_Class_and_PO_Zone = True 

# ***************************
# 	I have provided an exception to the "Standard Error Conditions"
# 	This exception will automatically adjust all 4 control values if any of them are negative or zero-valued.
# 	The "Standard Error Conditions" will typically skip calculations which don't meet the required conditions.
# 	The exception allows for de-bugging and developmental work.
# 	We check for whether or not the user specified to use "Standard Error Conditions" or not. 

	if standardErrorConditions:
		errEffectNumber = -1
		errFlagString = "" 
		RQ_place_holder = "" 
		errConditions = check_error_conditions(OW, PL, Al, Cu, RQ_place_holder, controlValues, print_detailed_err)

		isErr = errConditions.isErr
		err1 = errConditions.errors_list[0]
		err2 = errConditions.errors_list[1]
		err3 = errConditions.errors_list[2]
		err4 = errConditions.errors_list[3]

		if isErr:
			# This means an error flag exists and so we should report their effect.
			# Due to indexing, errConditions.second[0] is actually Error Type 1, and
			# errConditions.second[2] is Error Type 3.
			# Error 1 & 3 require us to ignore the data. 

			if err1 or err3:
				errEffectNumber = 2
			else:
				# If Error Types 1 or 3 weren't thrown, we don't need to skip the data.
				# This conditional should reflect that we have either Error Type 2 or 4.
				errEffectNumber = 1

			# We can build the error flag string
			if err1: errFlagString + "1"
			if err2: errFlagString + "2"
			if err3: errFlagString + "3"
			if err4: errFlagString + "4"
		else:
			errEffectNumber = 0 

		calculation_results.Error.value = errEffectNumber
		calculation_results.Error.descr = errFlagString 
	else:
		# we do a test to see if any values need to be shifted
		calculation_results.Error = test_for_negative_elements(OW, PL, Al, Cu)

		# Could also use this method:
		# calculation_results.Error = reassign_negative_and_zerovalued_element_values(OW, PL, Al, Cu).isErr

		# This Standard Error Condition Exception doesn't skip any data, so we store 1.
		calculation_results.Error.value = 1
		calculation_results.Error.descr = "NonStdErr"


# ***************************
# Branch Radiation Quality		
	Branch_RQ = radiation_quality(OW, PL, Al, Cu, b_UBeta_ENV_Class_and_PO_Zone, v_b_rad_ENV_classes, b_PO_Zone, standardErrorConditions)

	if print_individual_calculations:
		print(f"Branch_RQ= {Branch_RQ}")

	# Here we need to get an integer for the Radiation Quality as defined above.
	Branch_RQ_int = -1
	Branch_RQ_int = convert_RQ_string_to_int(Branch_RQ)

	calculation_results.Branch_RQ.value = Branch_RQ_int
	calculation_results.Branch_RQ.descr = Branch_RQ 

# ***************************
# DDE

# 	We need to consider if this is for an extremity zone
	if b_extremity_zone:
		DDE = 0.0 
	else:
		DDE = calc_total_DDE(OW, PL, Al, Cu, Branch_RQ, standardErrorConditions)

	if print_individual_calculations:
		print(f"DDE = {DDE}")

	calculation_results.DDE.value = DDE 

# ***************************
# SDE

	SDE = calc_total_SDE(OW, PL, Al, Cu, Branch_RQ, standardErrorConditions)

	calculation_results.SDE.value = SDE
	
	if print_individual_calculations:
		print(f"SDE = {SDE}")

# ***************************
# LDE
# 	We need to consider if this is for an extremity zone

	if b_extremity_zone:
		LDE = 0.0
	else:
		LDE = calc_LDE(DDE, SDE)
	
	if print_individual_calculations:
		print(f"LDE = {LDE}")	

	calculation_results.LDE.value = LDE

# ***************************
# Beta Indicator
# 	We have some custom inputs which force us down certain branches.
# 	If we have a forced beta source input then we know that the default value for Beta_Indicator has changed
# 	If the value of Beta_Indicator IS the default value (-1), then we know we have not already accounted for a custom input
	
	if Beta_Indicator is None:
		Beta_Indicator = assign_beta_indicator(OW, PL, Al, Cu, Branch_RQ)
	
	if print_individual_calculations:
		print(f"Beta_Indicator = {Beta_Indicator}")

	calculation_results.Beta_Indicator.value = Beta_Indicator

# ***************************
# Reported Radiation Quality
	Reported_RQ = report_radiation_quality(OW, PL, Al, Cu, Beta_Indicator, b_UBeta_ENV_Class_and_PO_Zone, v_b_rad_ENV_classes, b_PO_Zone, standardErrorConditions)
	
	if print_individual_calculations:
		print(f"Reported_RQ = {Reported_RQ}")

	# Here we need to get an integer for the Radiation Quality as defined above.
	Reported_RQ_int = -1
	Reported_RQ_int = convert_RQ_string_to_int(Reported_RQ)

	calculation_results.Reported_RQ.value = Reported_RQ_int
	calculation_results.Reported_RQ.descr = Reported_RQ


# ***************************
# Beta Dose
# 	Standard section for non-extremities
	Beta_Dose = calc_beta_dose_and_adjustedSDE(Beta_Indicator, Cu, DDE, SDE, Reported_RQ, b_extremity_zone, b_UBeta_ENV_Class_and_PO_Zone, v_b_rad_ENV_classes)

# Here we need to make an exception and check on how we report the Beta Dose when we have a Mixed_* source
# The issue arises when we truncated "Mixed_PL" to simply "PL". We lose information about the fact that the
# source had any beta in it. So, when we simply see "PL" we assume there is no beta in it.
# This can be confusing, so we need to adjust how we represent it. 

# 	if( Branch_RQ(0, 5) == "Mixed" ){
# 		// We konw that any "Mixed" RQ will be Reported without the "Mixed" label. So, we lose some relevant info
# 		// and we need to account for this in the reporting aspect of the doses. 
# 		// We are reading the RQ as if it has no beta components so... beta_dose -> 0
# 		// With beta_dose = 0, the entire reported (adjusted) SDE is just the total (pre)adjusted SDE
# 		Beta_Dose.first = {0.0, SDE};
# 	}

	calculation_results.Beta_Dose.value = Beta_Dose[0][0]
	calculation_results.Beta_Dose.descr = Beta_Dose[1]

	  # -*-*-*-*-*-*-*-*
	# 	We calculate an adjusted SDE by considering the beta dose (if there is any)
	# 	We also need to consider if this is for an extremity zone

	SDE_beta_dose_adjusted = Beta_Dose[0][1]

	# calculation_results.SDE.value = SDE_beta_dose_adjusted
	#    -*-*-*-*-*-*-*-*

# ***************************
# Energy
	Energy = calc_continuous_energy(OW, PL, Al, Cu, Reported_RQ, standardErrorConditions)

	calculation_results.Energy.value = Energy 

# ***************************

# ***************************
# Pre-Adjusted SDE
	calculation_results.SDE_preadjusted.value = SDE 

# ***************************
# ***************************

	return calculation_results

def main():
	print("Running full luxel algorithm...")
	calculation_results = luxel_plus_algorithm(155, 148, 137, 122)

	print_individual_calculations = True

	if print_individual_calculations:
		print(f"Error: {calculation_results.Error.value}")
		print(f"DDE: {calculation_results.DDE.value}")
		print(f"SDE: {calculation_results.SDE.value}")
		print(f"LDE: {calculation_results.LDE.value}")
		print(f"Branch_RQ: {calculation_results.Branch_RQ.descr}")
		print(f"Reported_RQ: {calculation_results.Reported_RQ.descr}")
		print(f"Beta_Indicator: {calculation_results.Beta_Indicator.value}")
		print(f"Beta_Dose: {calculation_results.Beta_Dose.value}")
		print(f"Energy: {calculation_results.Energy.value}")
		print(f"SDE_preadjusted: {calculation_results.SDE_preadjusted.value}")

if __name__ == "__main__":
    main()