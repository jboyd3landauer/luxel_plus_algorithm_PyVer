# Written by: John A. Boyd 3
# email: john.boydiii@landauer-fr.com
# Original Publish Date: 24 February, 2025

# Luxel+ Algorithm for Whole Body Dosimeters
# Revision 04601

# NOTE: This script is built upon the ROOT analysis framework.
# Certain variables and libraries may be ROOT-specific.
# ROOT Version: 6.24/08, built for win32

from .luxel_plus_source_type_indicators import *
from .luxel_plus_energy_functions import *
from .luxel_plus_radiation_quality_functions import *

def report_radiation_quality(OW, PL, Al, Cu, beta_indicator, b_PO_ENV_Class_and_UBeta_Zone = False, v_b_rad_ENV_classes = [False, False, False], b_PO_Zone = False, standardErrorConditions = True, return_special_rq = False):
	# Function to determine the Reporting Radiation Quality
	# This Reporting Radiation Quality is not used for any branching or decision purposes.
	# This is strictly used for reporting purposes which may provide valuable information.
	# -----------------------	

	# Vectors to hold ENV Class and PGC Zone booleans
	# 	ENV Class: {Photon only, UBeta, Xray Only}
	# 	PGC Zones: {Photon Only, Col & Waist, Extremity}

	if not standardErrorConditions:
		OW, PL, Al, Cu = check_for_negative_elements(OW, PL, Al, Cu)

	report_rq_output = ""
	beta_rq = ""
	photon_rq = ""
	mixed_rq = ""

	b_rad_env_class_DU = v_b_rad_ENV_classes[1]
	b_force_photon_only = False

	source_type = source_type_indicator(OW, PL, Al, Cu, b_rad_env_class_DU)

	# We first check the special case that we have a Pure Photon ENV Class AND UBeta Zone:
	if b_PO_ENV_Class_and_UBeta_Zone:
		beta_rq = beta_rad_quality(OW, PL, Al, CU, b_rad_env_class_DU)

		if( beta_rq == "-1" ):
			report_rq_values = determine_report_radiation_quality(OW, PL, Al, Cu, source_type, b_PO_ENV_Class_and_UBeta_Zone, b_rad_env_class_DU, beta_indicator)

	# Standard algorithm operation when an rq_input is not given.

	# Here we check for any Photon Only ENV Classes and/or Zones.
	# We only consider these as decision points IFF we don't also have a UBeta ENV Class
	# We must exclusively have a PO ENV Class or PG Zone in order to force this condition

	if (v_b_rad_ENV_classes[0] or v_b_rad_ENV_classes[2] or b_PO_Zone) and not v_b_rad_ENV_classes[1]:
		b_force_photon_only = False 
		source_type = 2
	else: 
		source_type = source_type_indicator(OW, PL, Al, Cu, b_rad_env_class_DU)

	# Lets put a check in for the case that the source_type is still 0 (no value) for some reason? 
	# If that is the case, we set the default value as a Mixed Source.	
	if source_type == 0:
		source_type = -1  # default to mixed
	
	if source_type == 1:
		beta_rq = beta_rad_quality(OW, PL, Al, Cu, b_rad_env_class_DU)

		if beta_rq == "BH_DU":
			beta_rq = "BU"
		if beta_rq == "Mixed_BH":
			beta_rq = "BH"
		report_rq_output = beta_rq
	else:
		report_rq_values = determine_report_radiation_quality(OW, PL, Al, Cu, source_type, b_PO_ENV_Class_and_UBeta_Zone, b_rad_env_class_DU, beta_indicator)
		report_rq_output = report_rq_values[0]
		beta_rq = report_rq_values[1]
		photon_rq = report_rq_values[2]

	if return_special_rq:
		report_rq_output = beta_rq

	return report_rq_output

def determine_report_radiation_quality(OW, PL, Al, Cu, source_type, b_PO_ENV_Class_and_UBeta_Zone, b_rad_env_class_DU, beta_indicator):

	photon_energy = calc_pure_photon_energy(OW, PL, Al, Cu)

	beta_rq = ""

	if source_type == 2:
		if photon_energy >= 205:
			photon_rq = "PH"
			report_rq_output = "PH"
		elif (photon_energy > 43) and (photon_energy < 204):
			photon_rq = "PM"
			report_rq_output = "PM"
		else:
			# PL: energy <= 40
			photon_rq = "PL"
			report_rq_output = "PL"
	else:
		# Here we have our Mixed Sources
		# To avoid any downstream and identification issues, we will rename "Mixed_XX" to simply "BH"

		# The idea is that when this hits Data Review, the presence of a DDE will kick in further logic
		# that passes the output along where it will be re-directed as a "PB"
		# "PB" stands for photon, beta mixture and so this is what we want.
		# *** We lose a bit of extra information about the Photon Only Rad ENV Class with a UBeta zone BUT...
		# that is fine because we can only access a limited set of RQ's into the Data Review Process: BL, BH, PL, PM, and PH	
	
		report_rq_output = "BH"

		# ***** For Dose Calc algorithm internal-use, only:
		if b_PO_ENV_Class_and_UBeta_Zone or b_rad_env_class_DU or beta_indicator:
			# Here we have the special case where we have a Photon Only Env Class but a UBeta Zone
			# Here we have detected a beta (beta_indicator = 1) BUT, it is not a Pure_DU
			# Therefore, we don't consider it to be a beta.
			# But we still need to consider this for the dose calc, so we calculate it as a Mixed Source
			# but we give it a new RQ of "P" for photon 
			beta_rq = "P"
			photon_rq = "P"
		else:
			beta_rq = "PB"
			photon_rq = "PB"

	rq_list = [report_rq_output, beta_rq, photon_rq]

	return rq_list