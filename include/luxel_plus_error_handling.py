from dataclasses import dataclass

@dataclass
class return_CVs_with_errFlag:
	CVs: float
	errFlag: bool

@dataclass
class error_check_values:
	isErr: bool
	errors_list: bool

def check_for_negative_elements(v1, v2, v3, v4):
	
	values = [v1, v2, v3, v4] 
	default_value = 1e-8 

	min_value = min(values)
	shift_value = 0.0

	if( min_value <= 0 ):
		shift_value = abs(min_value) + default_value
	elif( min_value == 0 ):
		shift_value = default_value 
	
	if( shift_value > 0 ):
		values = [v + shift_value for v in values]

	return tuple(values)

def test_for_negative_elements(v1, v2, v3, v4):
	# Function similar to "check_for_negative_elements" except that it won't modify the values.
	# This is just the test function of it only.

	negValsExist = False 
	values = [v1, v2, v3, v4] # Store values in this list
	minValue = min(values) # Find the minimum value

	if minValue <= 0:
		negValsExist = True 

	return negValsExist

def reassign_negative_and_zerovalued_element_values(v1, v2, v3, v4):
	# This is a test/check/adjustment of any element values which are negative or zero-valued
	# In this method we only shift the single element value which is negative or zero-valued
	# We shift the single value to some default value. All other values remain as they are. 

	hadErr = False 

	default_value = 1.0 

	# Only shift the individual values whcih need it to the default value
	if v1 <= 0:
		v1_return = default_value
		hadErr = True
	else:
		v1_return = v1 

	if v2 <= 0:
		v2_return = default_value
		hadErr = True
	else:
		v2_return = v2

	if v3 <= 0:
		v3_return = default_value
		hadErr = True
	else:
		v3_return = v3

	if v4 <= 0:
		v4_return = default_value
		hadErr = True 
	else: 
		v4_return = v4

	v_return = [v1_return, v2_return, v3_return, v4_return]

	# We are more or less superceding the previous/standard error_type1.
	# error_type1 would only return "true" if there was an error.
	# We do the same thing with "hadErr"
	return hadErr 

def check_error_type1(OW, PL, Al, Cu, print_detailed_err = False):
# *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
# 	1.	If a gross Converted value is <=0:
# 			set the gross Converted Value to 0 and treat as an error condition.  
# 			** NO DOSE CALCULATIONS ARE PERFORMED. **
	isErr = False 

	if (OW <= 0) or (PL <= 0) or (Al <= 0) or (Cu <= 0):
		isErr = True 

		if print_detailed_err:
			print("Error Condition 1.")

		if OW <= 0:
			OW = 0
			if print_detailed_err:
				print(f"OW = {OW}, ")
		if PL <= 0:
			PL = 0
			if print_detailed_err:
				print(f"PL = {PL}, ")
		if Al <= 0:
			Al = 0
			if print_detailed_err:
				print(f"Al = {Al}, ")
		if Cu <= 0:
			Cu = 0
			if print_detailed_err:
				print(f"Cu = {Cu}. ")

		if print_detailed_err:
			print("\n-----------------------------")
			print("Skipping calculation...")

		errFlag = True
	else:
		errFlag = False 

	CVs_with_errFlag = return_CVs_with_errFlag
	CVs_with_errFlag.CVs = [OW, PL, Al, Cu]
	CVs_with_errFlag.errFlag = errFlag 

	return CVs_with_errFlag

def check_error_type2(OW, PL, Al, Cu, controlValues, print_detailed_err = False):
# *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
# 2.	If a control value is subtracted from the gross Converted Value and the net Converted Value is <=1:
# 		set the net Converted Value to 1 and continue with the Luxel+ algorithm
# 		** Continue with the Luxel+ Algorithm. **

	isErr, errFlag = False, False
	errOW, errPL, errAl, errCu = False, False, False, False 

	default_value = 1.0

	# All of the Control Values should have already been checked to be non-zero/non-negative.
	# However, it can't hurt to perform an additional check and only proceed if all things look good.

	for CV in controlValues:
		if CV < 0: nonNegControlValues = False

	if nonNegControlValues:
		# We were given a non-negative vector of Control Values to check against the 4 given Converted Values
		if (OW - controlValues[0]) <= 1:
			isErr = True
			errOW = True
			OW = default_value
		if (PL - controlValues[1]) <= 1:
			isErr = True
			errPL = True
			PL = default_value
		if (Al - controlValues[2]) <= 1:
			isErr = True
			errPL = True
			Al = default_value
		if (Cu - controlValues[3]) <= 1:
			isErr = True
			errPL = True
			Cu = default_value

		if isErr and print_detailed_err:
			print("Error Condition 2.")
			if errOW:
				print(f"OW - CtrlVal err -- OW set to {default_value}, ")
			if errPL:
				print(f"PL - CtrlVal err -- PL set to {default_value}, ")
			if errAl:
				print(f"Al - CtrlVal err -- Al set to {default_value}, ")
			if errCu:
				print(f"Cu - CtrlVal err -- Cu set to {default_value}, ")

		errFlag = True


	CVs_with_errFlag = return_CVs_with_errFlag
	CVs_with_errFlag.CVs = [OW, PL, Al, Cu]
	CVs_with_errFlag.errFlag = errFlag 

	return CVs_with_errFlag

def check_error_type3(OW, PL, Al, Cu, print_detailed_err = False):
# *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
# 3.	We check various ratio limits and ranges and continue/kill the calculation accordingly:
# 		OW/Al<0.5 and abs(OW-Al) >10   or
# 		OW/Cu<0.5 and abs(OW-Cu) >10   or
# 		OW/PL<0.5 and abs(OW-PL) >10   or
# 		PL/Cu<0.5 and abs(PL-Cu) >10   or
# 		PL/Al<0.5 and abs(PL-Al) >10   or
# 		Cu/Al<0.5 and abs(Cu-Al) >10          
# 		------------------------------
# 		If any are true then: Attach a ratio error note and NO DOSE CALCULATIONS ARE PERFORMED.

	isErr = False
	errFlag = False

	ratio_limit = 0.5 
	sub_limit = 10.0 

	b_OW_Al, b_OWAl, b_OW_min_Al = False, False, False
	b_OW_Cu, b_OWCu, b_OW_min_Cu = False, False, False
	b_OW_PL, b_OWPL, b_OW_min_PL = False, False, False
	b_PL_Cu, b_PLCu, b_PL_min_Cu = False, False, False
	b_PL_Al, b_PLAl, b_PL_min_Al = False, False, False
	b_Cu_Al, b_CuAl, b_Cu_min_Al = False, False, False

	b_OWAl = ((OW/Al) < ratio_limit)
	b_OW_min_Al = (abs(OW-Al) > sub_limit)
	b_OW_Al = b_OWAl and b_OW_min_Al

	b_OWCu = ((OW/Cu) < ratio_limit);
	b_OW_min_Cu = (abs(OW-Cu) > sub_limit)
	b_OW_Cu = b_OWCu and b_OW_min_Cu;

	b_OWPL = ((OW/PL) < ratio_limit);
	b_OW_min_PL	= (abs(OW-PL) > sub_limit);
	b_OW_PL = b_OWPL and b_OW_min_PL;

	b_PLCu = ((PL/Cu) < ratio_limit);
	b_PL_min_Cu = (abs(PL-Cu) > sub_limit);
	b_PL_Cu = b_PLCu and b_PL_min_Cu;

	b_PLAl = ((PL/Al) < ratio_limit);
	b_PL_min_Al = (abs(PL-Al) > sub_limit);
	b_PL_Al = b_PLAl and b_PL_min_Al;

	b_CuAl = ((Cu/Al) < ratio_limit);
	b_Cu_min_Al = (abs(Cu-Al) > sub_limit);
	b_Cu_Al = b_CuAl and b_Cu_min_Al;

	if b_OW_Al or b_OW_Cu or b_OW_PL or b_PL_Cu or b_PL_Al or b_Cu_Al:
		isErr = True 
		if print_detailed_err:
			print("Error Condition 3.")

		if print_detailed_err:
			if b_OWAL: print("Ratio of OW/Al.")
			if b_OW_min_Al: print("Difference of OW - Al.")

			if b_OWCu: print("Ratio of OW/Cu.")
			if b_OW_min_Cu: print("Difference of OW - Cu.")

			if b_PLCu: print("Ratio of PL/Cu.")
			if b_PL_min_Cu: print("Difference of PL - Cu.")

			if b_PLAl: print("Ratio of PL/Al.")
			if b_PL_min_Al: print("Difference of PL - Al.")

			if b_CuAl: print("Ratio of Cu/Al.")
			if b_Cu_min_Al: print("Difference of Cu - Al.")

			print("\n-----------------------------")
			print("Skipping calculation...")

		errFlag = True 

	CVs_with_errFlag = return_CVs_with_errFlag
	CVs_with_errFlag.CVs = [OW, PL, Al, Cu]
	CVs_with_errFlag.errFlag = errFlag 

def check_error_type4(OW, PL, Al, Cu, RadQual, print_detailed_err = False):
# *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
# 4.	This is a check on low doses. 
# 		For doses between 1 and 20 we can define a general RQ of "P".
# 		The previous algorithm reported: SDE = Photon-Average-SDE and DDE = Photon-Average-DDE
# 		This algorithm doesn't calculate Photon-Average-SDe or Photon-Average-DDE
# 		Therefore, we will report the SDE and DDE per the usual algorithm determination.

	lower_limit = 10
	upper_limit = 20 

	errFlag = False
	inRangeOW = False
	inRangePL = False
	inRangeAl = False
	inRangeCu = False 

	inRangeOW = (lower_limit < OW) and (OW <= upper_limit)
	inRangePL = (lower_limit < PL) and (PL <= upper_limit)
	inRangeAl = (lower_limit < Al) and (Al <= upper_limit)
	inRangeCu = (lower_limit < Cu) and (Cu <= upper_limit)

	if inRangeOW and inRangePL and inRangeAl and inRangeCu:
		RadQual = "P"

		if print_detailed_err:
			print("Error Condition 4.")

		errFlag = True

	return errFlag 

def check_error_conditions(OW, PL, Al, Cu, RadQual, controlValues = [-1, -1, -1, -1], print_detailed_err = False):
	isErr = False
	err1 = False
	err2 = False
	err3 = False
	err4 = False 

	# Here we decide if we want to apply a check/test to adjust any negative and/or zero-valued element values to some positive default value

	adjust_negand_zeroVal_elements = True 

	if adjust_negand_zeroVal_elements:
		err1 = reassign_negative_and_zerovalued_element_values
	else:
		err1 = check_error_type1(OW, PL, Al, Cu, print_detailed_err).errFlag 

	# For Error 2, we only check this if a set of Control Values were given to test this condition. 
	# If Control Values were given then ALL values of the vector controlValues should be non-negative.
	# This will be our test condition for whether or not we should even consider this error condition. 

	for CV in controlValues:
		if CV < 0: nonNegControlValues = False

	if nonNegControlValues:
		err2 = check_error_type2(OW, PL, Al, Cu, print_detailed_err).errFlag 

	err3 = check_error_type3(OW, PL, Al, Cu, print_detailed_err)

	# For Error 4, we can only implement it if a proper Radiation Quality was already given
	# This error will overwrite the RQ and so we need a variable in memory to overwrite.

	if (RadQual is not None) or (RadQual != ""):
		err4 = check_error_type4(OW, PL, Al, Cu, RadQual, print_detailed_err)

	if err1 or err2 or err3 or err4:
		isErr = True 

	errors_list = [err1, err2, err3, err4]

	error_values = error_check_values
	error_values.isErr = isErr 
	error_values.errors_list = errors_list

	return error_values

def isValidRQinput(rq_input):
	rq_is_valid = False 

	validOptions = ["BL", "BH_DU", "BH", "PL", "PM", "PH", "Mixed_BL", "Mixed_BH"]

	if rq_input in validOptions: rq_is_valid = True 

	return rq_is_valid

