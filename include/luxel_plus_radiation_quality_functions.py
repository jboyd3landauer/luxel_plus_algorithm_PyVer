# Written by: John A. Boyd 3
# email: john.boydiii@landauer-fr.com
# Original Publish Date: 24 February, 2025

# Luxel+ Algorithm for Whole Body Dosimeters
# Revision 04601

# NOTE: This script is built upon the ROOT analysis framework.
# Certain variables and libraries may be ROOT-specific.
# ROOT Version: 6.24/08, built for win32

from .luxel_plus_source_type_tests import *
from .luxel_plus_source_type_indicators import *
from .luxel_plus_boundary_and_fit_params import *
from .luxel_plus_error_handling import check_for_negative_elements

def beta_rad_quality(OW, PL, Al, Cu, b_rad_env_class_DU = False, standardErrorConditions = True):
	# Function to determine the Beta Radiation Quality of a source.

	# Beta Radiation Quality:
	# BL: Low-energy Beta, Kr85
	# BH: High-energy Beta, Sr90
	# BH_DU: Depleted Uranium (or Uranium Slab), High-energy Beta
	# -----------------------
	if not standardErrorConditions:
		OW, PL, Al, Cu = check_for_negative_elements(OW, PL, Al, Cu)

	OWPL = None
	OWPL_limit = None

	beta_rq = ""

	# We need a special consideration for the Uranium Beta Environment Class
	# For this case, we only report that we have a beta IFF that beta is a Pure DU.
	# If it is not a pure DU then we move on.
	# We don't allow it be to a BL or BH. If NOT DU, then not Pure Beta. 

	if b_rad_env_class_DU:
		pure_DU_test = pure_beta_test(OW, PL, Al, Cu, True, b_rad_env_class_DU)

		if pure_DU_test is False:
			beta_rq = "-1"
		else:
			beta_rq = "BH_DU"
	else:
		pure_DU_test = pure_beta_test(OW, PL, Al, Cu, True, b_rad_env_class_DU)

		OWPL = OW/PL 
		OWPL_limit = 4

		BH_OWPL_Lower_Limit = 1.44

		CuOW = Cu/OW 
		CuOW_Beta_Upper_Limit = 0.125

		is_Pure_Beta = pure_beta_test(OW, PL, Al, Cu, False, b_rad_env_class_DU)

		if is_Pure_Beta:
			if pure_BH_test(OW, PL, Al, Cu):
				beta_rq = "BH"
			elif pure_DU_test:
				beta_rq = "BH_DU"
			else:
				beta_rq = "BL"
		else:
			beta_rq = "Mixed_BH"

	return beta_rq

def photon_rad_quality(OW, PL, Al, Cu, standardErrorConditions = True):
	# Function to determine the Photon Radiation Quality of a source

	# Photon Radiation Quality:
	# PL: Low-energy Pure Photons
	# PM: Medium-energy Pure Photons
	# PH: High-energy Pure Photons
	# -----------------------	

	if not standardErrorConditions:
		OW, PL, Al, Cu = check_for_negative_elements(OW, PL, Al, Cu)

	photon_rq = ""

	CuOW = Cu/OW 
	CuOW_Lower_limit = 0.155
	CuOW_Upper_Limit = 0.75

	OWPL = OW/PL 
	OWPL_Lower_limit = 0.25 
	OWPL_Upper_limit = 1.1 

	# Pure NS20 and M30 tests:

	if pure_photon_test(OW, PL, Al, Cu):
		# Added 10 Sep 2025
		# We can do a special check for Pure NS20 and Pure M30 that let's us jump straight to M30

		is_Pure_NS20 = pure_NS20_test(OW, PL, Al, Cu)
		is_Pure_M30 = pure_M30_test(OW, PL, Al, Cu)

		if is_Pure_NS20 or is_Pure_M30:
			# We can immediately define it as "PL" if it passes either of these tests
			photon_rq = "PL"
		else:
			if( CuOW_Upper_Limit <= CuOW ):
				photon_rq = "PH"
			elif (CuOW_Lower_limit < CuOW) and (CuOW < CuOW_Upper_Limit):
				photon_rq = "PM"
			else:
				#CuOW < CuOW_Lower_limit
				photon_rq = "PL"	
	else:
		photon_rq = "Mixed_BH"

	return photon_rq

def photon_or_mixed_beta_photon_rad_quality(OW, PL, Al, Cu):
	# Here is basically our standard RQ function
	# Priority of RQ ordering:
	# 1) Pure Beta: BH_DU, BL, BH
	# 2) Pure Photons: PH, PM, PL
	# 3) Mixed BetaPhoton: Mixed_BL, Mixed_BH;

	# The "fall-back"/default RQ is Mixed_BH;

	pure_beta_source = pure_beta_test(OW, PL, Al, Cu)

	if pure_beta_source:
		return beta_rad_quality(OW, PL, Al, Cu)
	else:
		pure_photon_source = pure_photon_test(OW, PL, Al, Cu)

		if pure_photon_source:
			return photon_rad_quality(OW, PL, Al, Cu)
		else:
			# Mixed Source
			OWPL = OW/PL
			CuOW = Cu/OW 
			boundary_line = -2.6786*CuOW + 2.8304 

			return "Mixed_BL" if OWPL > boundary_line else "Mixed_BH"


def radiation_quality(
    OW,
    PL,
    Al,
    Cu,
    b_PO_ENV_Class_and_UBeta_Zone=False,
    v_b_rad_ENV_classes=[False, False, False],
    b_PO_Zone=False,
    standardErrorConditions=True
):
	# Updated 26 June 2025

	# General function for determining the Radiation Quality of a source.
	# If source_type_input = 0, then this function will first classify the source as:
	# Pure Beta, Pure Photon, or Mixed.

	# Radiation Quality:
	# 	Pure Betas
	# 		BL: Low-energy Beta, Kr85
	# 		BH: High-energy Beta, Sr90
	# 		BH_DU: Depleted Uranium (or Uranium Slab), High-energy Beta	
	# 	Pure Photons
	# 		PL: Low-energy Pure Photons
	# 		PM: Medium-energy Pure Photons
	# 		PH: High-energy Pure Photons
	# 	Mixed Sources
	# 		Mixed_BL: Mixed Kr85 Betas with Pure Photons (typically Cs137)
	# 		Mixed_BH: Mixed Sr90 Betas with Pure Photons

	# -----------------------

	# Vectors to hold ENV Class and PGC Zone booleans
	# 	ENV Class: {Photon only, UBeta, Xray Only}
	# 	PGC Zones: {Photon Only, Col & Waist, Extremity}

    # Default parameter (cannot use mutable list as default directly)
    # if v_b_rad_ENV_classes is None:
    #     v_b_rad_ENV_classes = [False, False, False]

    # Optional error checks
    if not standardErrorConditions:
        OW, PL, Al, Cu = check_for_negative_elements(OW, PL, Al, Cu)

    rq_output = "Mixed_BH"
    source_type = 0
    pure_beta_source = False
    pure_photon_source = False

    b_rad_env_class_DU = v_b_rad_ENV_classes[1]
    b_force_photon_only = False

	# We first check the special case that we have a Pure Photon ENV Class AND a UBeta Zone:
    # --- Case 1: Pure Photon ENV Class and UBeta Zone
    if b_PO_ENV_Class_and_UBeta_Zone:
        beta_rq = beta_rad_quality(OW, PL, Al, Cu, b_rad_env_class_DU)
        if beta_rq == "-1":
            # Equivalent to 'goto photon_or_mixed'
            return determine_radiation_quality(OW, PL, Al, Cu)

	# If source_type_input != 0 then we take the value of source_type_input as the "source_type_indicator" value.
	# 	Source type as indicator:
	# 	-1: mixed
	# 	1: beta
	# 	2: pure photon


	# Here we check for any Photon Only ENV Classes and/or Zones.
	# We only consider these as decision points IFF we don't also have a UBeta ENV Class
	# We must exclusively have a PO ENV Class or PG Zone in order to force this condition

    # --- Case 2: Force photon only if conditions met
    if (v_b_rad_ENV_classes[0] or v_b_rad_ENV_classes[2] or b_PO_Zone) and (not v_b_rad_ENV_classes[1]):
        b_force_photon_only = True
        source_type = 2
        pure_beta_source = False
        pure_photon_source = True
    else:
        pure_beta_source = pure_beta_test(OW, PL, Al, Cu)

    # From here we can pass our CVs and booleans to our standard radiation determination function:
    rq_output = determine_radiation_quality(OW, PL, Al, Cu, pure_beta_source, pure_photon_source)

    # --- photon_or_mixed block (standard RQ logic)
    # if pure_beta_source:
    #     rq_output = beta_rad_quality(OW, PL, Al, Cu)
    # else:
    #     pure_photon_source = pure_photon_test(OW, PL, Al, Cu)
    #     if pure_photon_source:
    #         rq_output = photon_rad_quality(OW, PL, Al, Cu)
    #     else:
    #         # Mixed source region determination
    #         OWPL = OW / PL
    #         boundary_line = get_OWPL_yVal_for_boundary_line_between_MixedBH_and_MixedBL_from_CuOW(Cu, OW)
    #         if OWPL > boundary_line:
    #             rq_output = "Mixed_BL"
    #         else:
    #             rq_output = "Mixed_BH"

    return rq_output


def determine_radiation_quality(OW, PL, Al, Cu, pure_beta_source = False, pure_photon_source = False):
    """ Basic function for determining Radiation Quality
		This basically feeds the input CVs into the 3 possible main source types and returns the RQ from those
		We can feed some pre-determined booleans as well to help nudge the logic along....
    """
    pure_beta_source = pure_beta_test(OW, PL, Al, Cu)

    if pure_beta_source:
        return beta_rad_quality(OW, PL, Al, Cu)
    else:
        pure_photon_source = pure_photon_test(OW, PL, Al, Cu)
        if pure_photon_source == 1:
            return photon_rad_quality(OW, PL, Al, Cu)
        else:
            OWPL = OW / PL
            boundary_line = get_OWPL_yVal_for_boundary_line_between_MixedBH_and_MixedBL_from_CuOW(Cu, OW)
            return "Mixed_BL" if OWPL > boundary_line else "Mixed_BH"

def mixed_BetasPhotons_radiation_quality(OW, PL, Al, Cu, standardErrorConditions):
	# Function to determine the Radiation Quality for Mixed Betas and Photons
	
	# Mixed BetasPhotons RQ:
	# Mixed_BL: Mixed BL (Kr85) with Photons (typically Cesium)
	# Mixed_BH: Mixed BH (Sr90) with Photons (All)
	# -----------------------

    if not standardErrorConditions:
        OW, PL, Al, Cu = check_for_negative_elements(OW, PL, Al, Cu)

    mixed_rad_qual = ""

    MixedBH_source = False
    MixedBL_source = False 

    owpl = OW/PL 
    cuow = Cu/OW

    Mixed_BH_BL_boundary_line = get_OWPL_yVal_for_boundary_line_between_MixedBH_and_MixedBL_from_CuOW(Cu, OW)

    if( owpl > Mixed_BH_BL_boundary_line ):
    	MixedBL_source = True 
    	mixed_rad_qual = "Mixed_BL"
    else:
    	MixedBH_source = True 
    	mixed_rad_qual = "Mixed_BH"

    return mixed_rad_qual

def convert_RQ_string_to_int(RQ_string):
	# Convert RQ Strings to "Double Integers"

	# 0 <--> BL
	# 1 <--> BH_DU
	# 2 <--> BH
	# 3 <--> PL
	# 4 <--> PM
	# 5 <--> PH
	# 6 <--> Mixed_PL
	# 7 <--> Mixed_PM	

	RQ_num = -1 

	if RQ_string == "BL": RQ_num = 0
	if RQ_string == "BH_DU": RQ_num = 1
	if RQ_string == "BH": RQ_num = 2
	if RQ_string == "PL": RQ_num = 3
	if RQ_string == "PM": RQ_num = 4
	if RQ_string == "PH": RQ_num = 5
	if RQ_string == "Mixed_BH": RQ_num = 6
	if RQ_string == "Mixed_BL": RQ_num = 7
	if RQ_string == "P": RQ_num = 8

	return RQ_num 

def convert_RQ_int_to_string(RQ_int):
	# Convert RQ Ints to TStrings

	# 0 <--> BL
	# 1 <--> BH_DU
	# 2 <--> BH
	# 3 <--> PL
	# 4 <--> PM
	# 5 <--> PH
	# 6 <--> Mixed_BH
	# 7 <--> Mixed_BL

	RQ_string = ""

	if RQ_int == 0: RQ_string = "BL"
	if RQ_int == 1: RQ_string = "BH_DU"
	if RQ_int == 2: RQ_string = "BH"
	if RQ_int == 3: RQ_string = "PL"
	if RQ_int == 4: RQ_string = "PM"
	if RQ_int == 5: RQ_string = "PH"
	if RQ_int == 6: RQ_string = "Mixed_BH"
	if RQ_int == 7: RQ_string = "Mixed_BL"	
	if RQ_int == 8: RQ_string = "P"

	return RQ_string