This is the Python3 version of the Luxel+ Algorithm which was originally written in c++ (ROOT c++)

The main wrapper script is: jb_luxel_plus_algorithm.py

The primary calculation script is: jb_luxel_calc.py

These all utilized modules and functions in include/

A sample test script is provided: luxel_plus_data_analysis_script.py

This test script will read in data from a given file (input_filename), analyze it, and create an output file (output_excel_filename)
**NOTE that the format of the output file is French so, the delimiter is a semi-colon (;)

You will need to modify the load_data function such that the data file you read in has the same format.

For instance, in its current form, the load_data function reads in data with the following column structure:

sample_id source energy source_type sde dde ow pl al cu

i.e.,

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

where "data" is the values read in from the datafile. 

Depending on the format, and what needs to be changed in order to make this match your input data format, you may need to modify the following lines as well:

	for DoseCalc, samp_id, src_type, src, en, deep, shallow, in zip(dose_calcs, sample_id, source_type, source, energy, dde, sde):
		DoseCalc.Sample_ID = samp_id
		DoseCalc.Source_Type = src_type
		DoseCalc.Source = src 
		# DoseCalc.Ratio = rat
		DoseCalc.Energy = [en,""]
		DoseCalc.Known_DDE = deep
		DoseCalc.Known_SDE = shallow


In its current state, you can run: python luxel_plus_data_analysis_script.py

This will analyze test_input_data/Test3_PhotonOnlyConVertVal_WB_Photon_Only.csv

and create the following output: output_data/excel_Test3_PhotonOnlyConVertVal_WB_Photon_Only_PyVer.csv