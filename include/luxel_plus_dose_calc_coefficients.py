from .luxel_plus_boundary_and_fit_params import *
import math

# /////////////////////////////
# ////		SDE
# /////////////////////////////

def get_Mixed_BL_SDE_coefficients():
	# 14 Oct 2025
	# ***Use these coefficients***
	c1 = 2.19902
	c2 = -0.461181
	c3 = -0.390345
	c4 = -0.312412
	# Refined mean squared error: 159595 (399.493)

	Mixed_BL_SDE_coefficients = [c1, c2, c3, c4]

	return Mixed_BL_SDE_coefficients

def get_MixedBH_M30_SDE_coefficients():

	# # 14 Oct 2025
	# # ***Use these coefficients***
	# c1 = 1.31328
	# c2 = -1.27104
	# c3 = 0.240113
	# c4 = 3.79024
	# #Refined mean squared error: 2.44759e+06 (1564.48)

	# 17 Oct 2025
	c1 = 1.31327
	c2 = -1.27113
	c3 = 0.240387
	c4 = 3.78949
	#Refined mean squared error: 2.44788e+06 (1564.57)

	MixedBH_M30_SDE_coefficients = [c1, c2, c3, c4]

	return MixedBH_M30_SDE_coefficients

def get_MixedBH_NS20_SDE_coefficients():

	# 14 Oct 2025
	# ***Use these coefficients***
	c1 = 0.550254
	c2 = -0.771924
	c3 = 1.50246
	c4 = 3.79982
	#Refined mean squared error: 260018 (509.92)

	MixedBH_NS20_SDE_coefficients = [c1, c2, c3, c4]

	return MixedBH_NS20_SDE_coefficients 

def get_MixedBH_UpperEnergies_SDE_coefficients():

	# # 14 Oct 2025
	# # ***Use these coefficients***
	# c1 = 1.55362
	# c2 = -0.931387
	# c3 = -0.530041
	# c4 = 0.985584
	# #Refined mean squared error: 343887 (586.419)	

	# 17 Oct 2025
	c1 = 1.55768
	c2 = -0.943272
	c3 = -0.521712
	c4 = 0.984623
	#Refined mean squared error: 452230 (672.48)

	MixedBH_NonNS20_NonM30_SDE_coefficients = [c1, c2, c3, c4]

	return MixedBH_NonNS20_NonM30_SDE_coefficients

def get_Mixed_BH_SDE_coefficients(ow, pl, al, cu, known_source = None):

    # Using two boundaries to create 3 regions
    # We have 2 negatively sloped lines. 
    # The region to the left is for higher energy sources
    # The region between the lines is the "middle"
    # The region to the right is for lower energy sources.

    # y = ow/al and x = pl/al

    # The two lines are defined as:
    # Left (lower) line:
    # y = -1.75x + 5.1

    # Right (upper) line: 
    # y = -0.75x + 3.975

    # 1)
    # Here we are fitting to the left (higher energy) region
    # So, we cut everything to the right of the lower line.
    # It is a negatively sloped line so everything greater than the line is to the right

	# Line between Mixed_M30 and higher energy Mixed_BH
    # double lower_line_slope = -1.75;
    # // double lower_line_intercept = 5.1;
    # double lower_line_intercept = 5.05;
    # double lower_line = lower_line_slope*(pl/al) + lower_line_intercept;

	lower_line = get_OWAl_yVal_for_boundary_line_between_MixedBHM30_and_MixedBHHigherEnergies_from_PLAl(pl, al)

	if known_source is not None:
		if known_source == "M30":
			v_MixedBH_M30_coeffs = get_MixedBH_M30_SDE_coefficients()

			c1 = v_MixedBH_M30_coeffs[0]
			c2 = v_MixedBH_M30_coeffs[1]
			c3 = v_MixedBH_M30_coeffs[2]
			c4 = v_MixedBH_M30_coeffs[3]

		elif known_source == "NS20":
			v_MixedBH_NS20_coeffs = get_MixedBH_NS20_SDE_coefficients()
			c1 = v_MixedBH_NS20_coeffs[0]
			c2 = v_MixedBH_NS20_coeffs[1]
			c3 = v_MixedBH_NS20_coeffs[2]
			c4 = v_MixedBH_NS20_coeffs[3]	

		else:
			# This is the case where an "invalid" source type was provided.
			# In this case, we just default to the "Upper Energies" MixedBH coefficients:	

			v_MixedBH_UpperEnergies_coeffs = get_MixedBH_UpperEnergies_SDE_coefficients()
			c1 = v_MixedBH_UpperEnergies_coeffs[0]
			c2 = v_MixedBH_UpperEnergies_coeffs[1]
			c3 = v_MixedBH_UpperEnergies_coeffs[2]
			c4 = v_MixedBH_UpperEnergies_coeffs[3]
	else:
		if (ow/al) < lower_line:
        	# Values in lower region: Higher energy sources (lower meaning, below the line)
            
            # These are values to the left of this lower line
            
        	# We can also try and make a cut on another range of sources.
        	# Here we cut on the (OW/PL) vs (OW/Cu) plane 

			v_MixedBH_UpperEnergies_coeffs = get_MixedBH_UpperEnergies_SDE_coefficients()
			c1 = v_MixedBH_UpperEnergies_coeffs[0]
			c2 = v_MixedBH_UpperEnergies_coeffs[1]
			c3 = v_MixedBH_UpperEnergies_coeffs[2]
			c4 = v_MixedBH_UpperEnergies_coeffs[3]
		else:
            # Line between Mixed_NS20 and Mixed_M30
            # double upper_line_slope = -0.75;
            # double upper_line_intercept = 3.975;
            
            # Updated 13 October 2025 after updating mixture calculations/methods:
            # double upper_line_slope = -0.3405;
            # double upper_line_intercept = 3.2191;

            # double upper_line = upper_line_slope*(pl/al) + upper_line_intercept;

			upper_line = get_OWAl_yVal_for_boundary_line_between_MixedBHNS20_and_MixedBHM30_from_PLAl(pl, al)

            # middle region is where we are greater than the lower line AND less than the upper line.
            # Therefore, we are in the lower-right overlapping region of the two lines. 
            # Here we only consider values in the middle region

			if (ow/al) >= lower_line and (ow/al) <= upper_line:
            	# M30
				v_MixedBH_M30_coeffs = get_MixedBH_M30_SDE_coefficients()
				c1 = v_MixedBH_M30_coeffs[0]
				c2 = v_MixedBH_M30_coeffs[1]
				c3 = v_MixedBH_M30_coeffs[2]
				c4 = v_MixedBH_M30_coeffs[3]
			else:
				# NS20
                # Here we have values in the right region
                # We take everything above the upper line
                # We can assume we are meeting this condition: (ow/al) > upper_line )

				v_MixedBH_NS20_coeffs = get_MixedBH_NS20_SDE_coefficients()
				c1 = v_MixedBH_NS20_coeffs[0]
				c2 = v_MixedBH_NS20_coeffs[1]
				c3 = v_MixedBH_NS20_coeffs[2]
				c4 = v_MixedBH_NS20_coeffs[3]               

	Mixed_BH_SDE_coefficients = [c1, c2, c3, c4]

	return Mixed_BH_SDE_coefficients

# /////////////////////////////
# ////		DDE
# /////////////////////////////

def get_Mixed_BL_DDE_coefficients():
	# 14 Oct 2025
	# ***Use these coefficients***
	c1 = -0.00184938
	c2 = 0.401851
	c3 = 0.33279
	c4 = 0.267751
	#Refined mean squared error: 3133.31 (55.976)

	Mixed_BL_DDE_coefficients = [c1, c2, c3, c4]

	return Mixed_BL_DDE_coefficients

def get_MixedBH_M30_DDE_coefficients():
	# # 14 Oct 2025
	# # ***Use these coefficients***
	# c1 = -0.120131
	# c2 = 0.244997
	# c3 = 0.025124
	# c4 = -0.466047
	# #Refined mean squared error: 84546.2 (290.768)

	# 17 Oct 2025
	c1 = -0.120155
	c2 = 0.24501
	c3 = 0.0251505
	c4 = -0.465947
	#Refined mean squared error: 84558.8 (290.79)

	MixedBH_M30_DDE_coefficients = [c1, c2, c3, c4]

	return MixedBH_M30_DDE_coefficients

def get_MixedBH_NS20_DDE_coefficients():
	# 14 Oct 2025
	# ***Use these coefficients***
	c1 = -0.0254205
	c2 = 0.171487
	c3 = -0.128046
	c4 = -0.392067
	#Refined mean squared error: 4692.08 (68.4987)

	MixedBH_NS20_DDE_coefficients = [c1, c2, c3, c4]

	return MixedBH_NS20_DDE_coefficients

def get_MixedBH_UpperEnergies_DDE_coefficients():
	# # 14 Oct 2025
	# # ***Use these coefficients***
	# c1 = -0.192212
	# c2 = 0.34852
	# c3 = -0.0522534
	# c4 = 0.742162
	# #Refined mean squared error: 35132.6 (187.437)

	# 17 Oct 2025
	c1 = -0.193218
	c2 = 0.350681
	c3 = -0.0533068
	c4 = 0.742411
	#Refined mean squared error: 42214 (205.46)

	Mixed_BH_UpperEnergies_DDE_coefficients = [c1, c2, c3, c4]

	return Mixed_BH_UpperEnergies_DDE_coefficients

def get_Mixed_BH_DDE_coefficients(ow, pl, al, cu, known_source = None):
    # Using two boundaries to create 3 regions
    # We have 2 negatively sloped lines. 
    # The region to the left is for higher energy sources
    # The region between the lines is the "middle"
    # The region to the right is for lower energy sources.

    # y = ow/al and x = pl/al

    # The two lines are defined as:
    # Left (lower) line:
    # y = -1.75x + 5.1

    # Right (upper) line: 
    # y = -0.75x + 3.975

    # 1)
    # Here we are fitting to the left (higher energy) region
    # So, we cut everything to the right of the lower line.
    # It is a negatively sloped line so everything greater than the line is to the right

	# Line between Mixed_M30 and higher energy Mixed_BH
    # double lower_line_slope = -1.75;
    # // double lower_line_intercept = 5.1;
    # double lower_line_intercept = 5.05;
    # double lower_line = lower_line_slope*(pl/al) + lower_line_intercept;

	lower_line = get_OWAl_yVal_for_boundary_line_between_MixedBHM30_and_MixedBHHigherEnergies_from_PLAl(pl, al)

	if known_source is not None:
		if known_source == "M30":
			v_MixedBH_M30_DDE_coeffs = get_MixedBH_M30_DDE_coefficients()
			c1 = v_MixedBH_M30_DDE_coeffs[0]
			c2 = v_MixedBH_M30_DDE_coeffs[1]
			c3 = v_MixedBH_M30_DDE_coeffs[2]
			c4 = v_MixedBH_M30_DDE_coeffs[3]
		elif known_source == "NS20":
			v_MixedBH_NS20_DDE_coeffs = get_MixedBH_NS20_DDE_coefficients()
			c1 = v_MixedBH_NS20_DDE_coeffs[0]
			c2 = v_MixedBH_NS20_DDE_coeffs[1]
			c3 = v_MixedBH_NS20_DDE_coeffs[2]
			c4 = v_MixedBH_NS20_DDE_coeffs[3]
		else:
    		# This is the case where an "invalid" source type was provided.
    		# In this case, we just default to the "Upper Energies" MixedBH coefficients:			
			v_MixedBH_UpperEnergies_DDE_coeffs = get_MixedBH_UpperEnergies_DDE_coefficients()
			c1 = v_MixedBH_UpperEnergies_DDE_coeffs[0]
			c2 = v_MixedBH_UpperEnergies_DDE_coeffs[1]
			c3 = v_MixedBH_UpperEnergies_DDE_coeffs[2]
			c4 = v_MixedBH_UpperEnergies_DDE_coeffs[3]
	else:
		if (ow/al) < lower_line:
        	# Values in lower region: Higher energy sources (lower meaning, below the line)
            
            # These are values to the left of this lower line
            
        	# We can also try and make a cut on another range of sources.
        	# Here we cut on the (OW/PL) vs (OW/Cu) plane 
			v_MixedBH_UpperEnergies_DDE_coeffs = get_MixedBH_UpperEnergies_DDE_coefficients()
			c1 = v_MixedBH_UpperEnergies_DDE_coeffs[0]
			c2 = v_MixedBH_UpperEnergies_DDE_coeffs[1]
			c3 = v_MixedBH_UpperEnergies_DDE_coeffs[2]
			c4 = v_MixedBH_UpperEnergies_DDE_coeffs[3]
		else:
			# Line between Mixed_NS20 and Mixed_M30
			upper_line = get_OWAl_yVal_for_boundary_line_between_MixedBHNS20_and_MixedBHM30_from_PLAl(pl, al)

            # middle region is where we are greater than the lower line AND less than the upper line.
            # Therefore, we are in the lower-right overlapping region of the two lines. 
            # Here we only consider values in the middle region
			if (ow/al) >= lower_line and (ow/al) <= upper_line:
            	# M30
				v_MixedBH_M30_DDE_coeffs = get_MixedBH_M30_DDE_coefficients()
				c1 = v_MixedBH_M30_DDE_coeffs[0]
				c2 = v_MixedBH_M30_DDE_coeffs[1]
				c3 = v_MixedBH_M30_DDE_coeffs[2]
				c4 = v_MixedBH_M30_DDE_coeffs[3]
			else:
            	# NS20
                # Here we have values in the right region
                # We take everything above the upper line
                # We can assume we are meeting this condition: (ow/al) > upper_line )
				v_MixedBH_NS20_DDE_coeffs = get_MixedBH_NS20_DDE_coefficients()
				c1 = v_MixedBH_NS20_DDE_coeffs[0]
				c2 = v_MixedBH_NS20_DDE_coeffs[1]
				c3 = v_MixedBH_NS20_DDE_coeffs[2]
				c4 = v_MixedBH_NS20_DDE_coeffs[3]

	Mixed_BH_DDE_coefficients = [c1, c2, c3, c4]

	return Mixed_BH_DDE_coefficients