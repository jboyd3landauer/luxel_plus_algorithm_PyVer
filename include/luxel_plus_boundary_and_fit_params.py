def get_OWAl_yVal_for_boundary_line_between_MixedBHM30_and_MixedBHHigherEnergies_from_PLAl(pl, al):
	PLAl = pl/al

	# Line between Mixed_M30 and higher energy Mixed_BH
	M30_and_HigherEnergies_boundary_line_slope = -1.75
	M30_and_HigherEnergies_boundary_line_intercept = 5.05
	M30_and_HigherEnergies_boundary_line = M30_and_HigherEnergies_boundary_line_slope*(PLAl) + M30_and_HigherEnergies_boundary_line_intercept

	return M30_and_HigherEnergies_boundary_line

def get_OWAl_yVal_for_boundary_line_between_MixedBHNS20_and_MixedBHM30_from_PLAl(pl, al):
	PLAl = pl/al 

	# Updated 13 October 2025 after updating mixture calculations/methods
	MixedBHNS20_and_MixedBHM30_boundary_line_slope = -0.3405
	MixedBHNS20_and_MixedBHM30_boundary_line_intercept = 3.2191
	MixedBHNS20_and_MixedBHM30_boundary_line = MixedBHNS20_and_MixedBHM30_boundary_line_slope*(PLAl) + MixedBHNS20_and_MixedBHM30_boundary_line_intercept

	return MixedBHNS20_and_MixedBHM30_boundary_line

def get_OWPL_yVal_for_boundary_line_between_MixedBH_and_MixedBL_from_CuOW(cu, ow):
	CuOW = cu/ow 

	MixedBH_and_MixedBL_boundary_line_slope = -2.6786
	MixedBH_and_MixedBL_boundary_line_intercept = 2.8304
	MixedBH_and_MixedBL_boundary_line = MixedBH_and_MixedBL_boundary_line_slope*(CuOW) + MixedBH_and_MixedBL_boundary_line_intercept

	return MixedBH_and_MixedBL_boundary_line