global fit_MixedBH_NS20, fit_MixedBH_M30, fit_MixedBH_UpperEnergies
fit_MixedBH_NS20 = False
fit_MixedBH_M30 = False
fit_MixedBH_UpperEnergies = False

global separate_BH_and_PureDU
separate_BH_and_PureDU = True

global use_v2, pure_beta_cut, mixed_beta_cut, mixed_BL_BH_Cut, cut_MixedBL, mixed_beta_low_cut, mixed_beta_high_cut, mixed_BetaPhotons_super_cuts
use_v2, pure_beta_cut, mixed_beta_cut, mixed_BL_BH_Cut, cut_MixedBL, mixed_beta_low_cut, mixed_beta_high_cut, mixed_BetaPhotons_super_cuts = [False]*8

global n_all_sources, n_pure_photons, n_pure_Kr85, n_pure_Sr90, n_MixedBL_SDE, n_MixedBL_DDE, n_MixedBH_SDE, n_MixedBH_DDE, n_pure_DU, n_MixedBL, n_MixedBH
n_all_sources, n_pure_photons, n_pure_Kr85, n_pure_Sr90, n_MixedBL_SDE, n_MixedBL_DDE, n_MixedBH_SDE, n_MixedBH_DDE, n_pure_DU, n_MixedBL, n_MixedBH = [False]*11

global n_MixedPhotons_Gen, n_MixedPhotons_HighE, n_MixedPhotons_MedE
n_MixedPhotons_Gen, n_MixedPhotons_HighE, n_MixedPhotons_MedE = [False]*3

global n_all_sources, n_pure_photons, n_pure_Kr85, n_pure_Sr90, n_MixedBL_SDE, n_MixedBL_DDE, n_MixedBH_SDE, n_MixedBH_DDE, n_pure_DU
n_all_sources, n_pure_photons, n_pure_Kr85, n_pure_Sr90, n_MixedBL_SDE, n_MixedBL_DDE, n_MixedBH_SDE, n_MixedBH_DDE, n_pure_DU = [False]*9

global n_MixedBL, n_MixedBH, n_MixedPhotons_Gen, n_MixedPhotons_HighE, n_MixedPhotons_MedE
n_MixedBL, n_MixedBH, n_MixedPhotons_Gen, n_MixedPhotons_HighE, n_MixedPhotons_MedE = [False]*5


def within_valid_dose_range(source_type, dose_value, dose_type):
	# Dose values in mrem

	within_range = False 

	if source_type == "Pure_Photons":
		# 5 to 500 rad
		dose_min = 5000
		dose_max = 500000

	if source_type == "Photon_Photon_Mixtures":
		# 0.05 to 5 rem
		dose_min = 50
		dose_max = 5000

	if source_type == "Photon_Beta_Mixtures":
		# SDE: 0.3 to 30 rem
		if dose_type == "SDE":
			dose_min = 300
			dose_max = 30000

		# DDE: 0.05 to 5 rem
		if dose_type == "DDE":
			dose_min = 50
			dose_max = 5000
	if (dose_value >= dose_min) and (dose_value <= dose_max):
		within_range = True 

	return within_range
