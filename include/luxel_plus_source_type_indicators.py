# Written by: John A. Boyd 3
# email: john.boydiii@landauer-fr.com
# Original Publish Date: 24 February, 2025

# Luxel+ Algorithm for Whole Body Dosimeters
# Revision 04601

# NOTE: This script is built upon the ROOT analysis framework.
# Certain variables and libraries may be ROOT-specific.
# ROOT Version: 6.24/08, built for win32

from .luxel_plus_error_handling import check_for_negative_elements
from .luxel_plus_source_type_tests import pure_BL_test, pure_BH_test, pure_NS20_test, pure_M30_test, MixedBH_M30_test, MixedBH_NS20_test, pure_beta_test, pure_photon_test, mixed_source_test

def assign_beta_indicator(OW, PL, Al, Cu, known_rq_input = "", standardErrorConditions = True):
	# Assigns a Beta Indicator based on a known RQ
	# Will be true (1) for Pure Beta and Mixed Sources
	# -----------------------

	if not standardErrorConditions:
		OW, PL, Al, Cu = check_for_negative_elements(OW, PL, Al, Cu)

	betaIndicator = False 

	if known_rq_input.startswith("PL") or known_rq_input.startswith("PM") or known_rq_input.startswith("PH"):
		betaIndicator = False 
	if known_rq_input.startswith("Mixed"):
		betaIndicator = True 
	if known_rq_input.startswith("BL") or known_rq_input.startswith("BH"):
		betaIndicator = True

	return betaIndicator

def beta_indicator(OW, PL, Al, Cu, known_rq_input = "", b_rad_env_class_DU = False, standardErrorConditions = True):
	# General Beta Indicator
	# Will be true (1) for Pure Beta and Mixed Sources

	if not standardErrorConditions:
		OW, PL, Al, Cu = check_for_negative_elements(OW, PL, Al, Cu)

	betaIndicator = False

	if known_rq_input.startswith("PL") or known_rq_input.startswith("PM") or known_rq_input.startswith("PH"):
		betaIndicator = False
	if known_rq_input.startswith("Mixed"):
		betaIndicator = True

	# We need a special consideration for the Uranium Beta Environment Class
	# For this case, we only report that we have a beta IFF that beta is a Pure DU.
	# If it is not a pure DU then we move on.
	# We don't allow it be to a BL or BH. If NOT DU, then not Pure Beta. 

	beta_test = False 
	if b_rad_env_class_DU:
		pure_DU_test = pure_beta_test(OW, PL, Al, Cu, True, b_rad_env_class_DU)

		if pure_DU_test is False:
			beta_test = False 
		else:
			beta_test = pure_DU_test
	else:
		beta_test = pure_beta_test(OW, PL, Al, Cu, False, b_rad_env_class_DU)

	photon_test = pure_photon_test(OW, PL, Al, Cu)
	mixed_test = mixed_source_test(OW, PL, Al, Cu)

	if beta_test or mixed_test:
		betaIndicator = True

	return betaIndicator

def source_type_indicator(OW, PL, Al, Cu, b_rad_env_class_DU = False, standardErrorConditions = True):
	# Function to classify a given source as Pure Beta, Pure Photon, or Mixed Source.
	# -----------------------	
	if not standardErrorConditions:
		OW, PL, Al, Cu = check_for_negative_elements(OW, PL, Al, Cu)

	source_type_indicator_output = 0

	beta_test = False
	photon_test = False
	mixed_test = False 
	pure_DU_test = False

	# We need a special consideration for the Uranium Beta Environment Class
	# For this case, we only report that we have a beta IFF that beta is a Pure DU.
	# If it is not a pure DU then we move on.
	# We don't allow it be to a BL or BH. If NOT DU, then not Pure Beta. 

	if b_rad_env_class_DU:
		pure_DU_test = pure_beta_test(OW, PL, Al, Cu, True, b_rad_env_class_DU)

		if pure_DU_test is False:
			beta_test = False
		else:
			beta_test = True
	else:
		beta_test = pure_beta_test(OW, PL, Al, Cu, False, b_rad_env_class_DU)

	photon_test = pure_photon_test(OW, PL, Al, Cu)
	mixed_test = mixed_source_test(beta_test, photon_test)

	if beta_test:
		source_type_indicator_output = 1
	else:
		if photon_test:
			source_type_indicator_output = 2
		else:
			source_type_indicator_output = -1

	return source_type_indicator_output