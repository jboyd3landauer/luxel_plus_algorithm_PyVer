import jb_luxel_plus_algorithm
import pandas as pd
import numpy as np
import csv
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
from jb_luxel_calc import calc_continuous_energy
from jb_luxel_calc import reassign_negative_and_zerovalued_element_values
from jb_luxel_plus_algorithm import luxel_plus_algorithm

input_filename = ""
output_cpp_filename = ""
output_excel_filename = ""

def fr_number(x):
    """Format float with comma as decimal separator."""
    if isinstance(x, float):
        return f"{x:.6f}".replace('.', ',')  # adjust precision as needed
    return str(x)  # leave strings unchanged

def load_data( filename, outfile, excel_outfile):

	# read the specified file:
	# data_file = pd.read_csv(filename, sep=',')
	data_file = pd.read_csv(filename, sep='\s+')
	data = data_file.to_numpy()

	# Convert all columns to NumPy arrays at once:

	# sample_id  = data[:, 0]   
	# source_type = data[:, 1]
	# source = data[:, 2]
	# ratio = data[:, 3]
	# ow = data[:, 4]
	# pl = data[:, 5]
	# al = data[:, 6]
	# cu = data[:, 7]
	# dde = data[:, 8]
	# sde = data[:, 9]

	# Test1_SimTestConVertVal_WB_SR90.csv
	sample_id  = data[:, 0]   
	source = data[:, 1]
	energy = data[:, 2]
	source_type = data[:, 3]
	sde = data[:, 4]
	dde = data[:, 5]
	ow = data[:, 6]
	pl = data[:, 7]
	al = data[:, 8]
	cu = data[:, 9]

	# CVs = np.array([reassign_negative_and_zerovalued_element_values(o, p, a, c) for o, p, a, c in zip(ow, pl, al, cu)])

	dose_calcs = np.array([luxel_plus_algorithm(o, p, a, c) for o, p, a, c in zip(ow, pl, al, cu)])

	for DoseCalc, samp_id, src_type, src, en, deep, shallow, in zip(dose_calcs, sample_id, source_type, source, energy, dde, sde):
		DoseCalc.Sample_ID = samp_id
		DoseCalc.Source_Type = src_type
		DoseCalc.Source = src 
		# DoseCalc.Ratio = rat
		DoseCalc.Energy = [en,""]
		DoseCalc.Known_DDE = deep
		DoseCalc.Known_SDE = shallow

	return dose_calcs

def main():
	# input_filename = ".\\test_input_data\\Test1_SimTestConVertVal_WB_SR90.csv"
	# input_filename = ".\\test_input_data\\Test2_SimTestConVertVal_WB_UBETA.csv"
	# input_filename = ".\\test_input_data\\Test3_PhotonOnlyConVertVal_WB_Photon_Only.csv"
	input_filename = ".\\test_input_data\\MSalasky_low_dose_test_data.csv"

	print(f"Input filename: {input_filename}")

	# output_cpp_filename = "./cpp_NS20_M30_luxel_plus_data_analysis_14Oct20255_pyVer.csv"
	# output_cpp_filename = "./cpp_Test1_SimTestConVertVal_WB_SR90_pyVer.csv"
	output_cpp_filename = "./cpp_MSalasky_low_dose_test_data_pyVer.csv"

	# output_excel_filename = ".\\output_data\\excel_Test1_SimTestConVertVal_WB_SR90_pyVer_PyVer.csv"
	# output_excel_filename = ".\\output_data\\excel_Test2_SimTestConVertVal_WB_UBETA_PyVer.csv"
	output_excel_filename = ".\\output_data\\excel_MSalasky_low_dose_test_data_PyVer.csv"

	dose_calcs = load_data(input_filename, output_cpp_filename, output_excel_filename)

	headers = ["Sample_ID", "Source_Type", "Source", "OW","PL","Al","Cu","Known_DDE","DDE_calc","Reported_DDE","Known_SDE","SDE_calc","Reported_SDE","LDE_calc","Reported_LDE","Branch_RQ","Reported_RQ"]

	with open(output_excel_filename, "w", newline='') as f:
		writer = csv.writer(f, delimiter=',')
		writer.writerow(headers)

		for dc in dose_calcs:
			# writer.writerow([dc.Sample_ID, dc.Source_Type, dc.Source, dc.Ratio, fr_number(dc.Known_OW), fr_number(dc.Known_PL), fr_number(dc.Known_Al), fr_number(dc.Known_Cu), fr_number(dc.Known_DDE), fr_number(dc.DDE.value), fr_number(dc.Known_SDE), fr_number(dc.SDE_preadjusted.value), dc.Branch_RQ.descr, dc.Reported_RQ.descr])
			writer.writerow([dc.Sample_ID, dc.Source_Type, dc.Source, dc.Known_OW, dc.Known_PL, dc.Known_Al, dc.Known_Cu, dc.Known_DDE, dc.DDE.value, dc.Reported_DDE.value, dc.Known_SDE, dc.SDE_preadjusted.value, dc.Reported_SDE.value, dc.LDE.value, dc.Reported_LDE.value, dc.Branch_RQ.descr, dc.Reported_RQ.value])
	
	print(f"Known_DDE = {dose_calcs[0].Known_DDE}")
	print(f"DDE = {dose_calcs[0].DDE.value}")

	print(f"Known_SDE = {dose_calcs[0].Known_SDE}")
	print(f"SDE = {dose_calcs[0].SDE_preadjusted.value}")

	print(f"Branch_RQ = {dose_calcs[0].Branch_RQ.descr}")
	print(f"Reported_RQ = {dose_calcs[0].Reported_RQ.descr}")

	print(f"----------------")
	print(f"Excel output filename: {output_excel_filename}")
	print(f"----------------")

if __name__ == "__main__":
    main()
