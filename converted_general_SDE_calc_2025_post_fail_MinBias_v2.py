#!/usr/bin/env python3
"""
Pure Python replacement for C++ ROOT TMinuit fitting script.

- Loads multiple space-delimited CSV files.
- Ignores inconsistent headers (skips first row).
- Assigns fixed columns: source, ow, pl, al, cu, dde, sde.
- Fits modelFunction to SDE (or DDE if desired).
- Plots known vs calculated dose (Nominal & Refined fits).
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from pathlib import Path

import jb_luxel_calc
import math

from include.luxel_plus_radiation_quality_functions import radiation_quality
from include.luxel_plus_source_type_tests import (
    pure_photon_test,
    pure_beta_test,
    pure_NS20_test,
    pure_M30_test,
    pure_BL_test,
    pure_BH_test,
    MixedBH_NS20_test,
    MixedBH_M30_test,
    mixed_source_test,
)

global known_source_type, lookup_RQ

# Initialize empty global dataset (will be populated by load_data)
GLOBAL_DATA = pd.DataFrame(
    columns=["source", "ow", "pl", "al", "cu", "dde", "sde", "rad_quality"]
)

def add_data(new_df):
    """Append new data to the global dataset."""
    global GLOBAL_DATA
    if new_df is None or new_df.empty:
        return
    GLOBAL_DATA = pd.concat([GLOBAL_DATA, new_df], ignore_index=True)

def get_data():
    """Return a reference to the global dataset."""
    global GLOBAL_DATA
    return GLOBAL_DATA

def clear_data():
    """Clear the global dataset."""
    global GLOBAL_DATA
    GLOBAL_DATA = GLOBAL_DATA.iloc[0:0]  # reset but keep structure

# ============================================================
# Fit bounds per radiation type
# ============================================================

def get_fit_bounds(known_type):
    """
    Returns bounds for parameters (A, B, C, D)
    based on the radiation source type.
    """
    known_type = known_type.upper()

    if known_type == "BL":     # Kr85 Beta Low
        return [(0.5, 5.0), (0.0, 0.001), (0.0, 0.001), (0.0, 0.001)]
    elif known_type == "BH":   # Sr90 Beta High
        return [(-10.0, 10.0), (-10.0, 10.0), (-10.0, 10.0), (-10.0, 10.0)]
    elif known_type == "DU":   # Depleted Uranium
        return [(0.0, 2.0), (-1.0, 1.0), (-5.0, 1.0), (0.0, 35.0)]
    elif known_type == "PH_SDE":   # PH SDE
        return [(0.0, 5.0), (0.0, 5.0), (0.0, 5.0), (0.0, 5.0)]
    elif known_type == "PH_DDE":   # PH DDE
        return [(0.0, 5.0), (0.0, 5.0), (0.0, 5.0), (0.0, 5.0)]
    elif known_type == "PL_PM_SDE":   # PL or PM SDE
        return [(-5.0, 5.0), (-5.0, 5.0), (-5.0, 5.0), (-5.0, 5.0)]
    elif known_type == "PL_PM_DDE":   # PL or PM DDE
        return [(0.0, 5.0), (0.0, 5.0), (0.0, 5.0), (0.0, 5.0)]
    elif known_type == "Mixed_BL_SDE":   # Mixed_BL SDE
        return [(-5.0, 5.0), (-5.0, 5.0), (-5.0, 5.0), (-5.0, 5.0)]
    elif known_type == "Mixed_BL_DDE":   # Mixed_BL DDE
        return [(-5.0, 5.0), (-5.0, 5.0), (-5.0, 5.0), (-5.0, 5.0)]
    elif known_type == "Mixed_BH_SDE":   # Mixed_BH SDE
        return [(-5.0, 5.0), (-5.0, 5.0), (-5.0, 5.0), (-5.0, 15.0)]
    elif known_type == "Mixed_BH_DDE":   # Mixed_BH DDE
        return [(-5.0, 5.0), (-5.0, 5.0), (-5.0, 5.0), (-5.0, 5.0)]
    else:  # Fallback default
        return [(-5.0, 5.0), (-5.0, 5.0), (-5.0, 5.0), (-5.0, 5.0)]


# ============================================================
# 1. Utility Functions
# ============================================================

def sort_vectors_by_primary(primary, *vectors):
    """Sort multiple vectors by the values in the primary vector."""
    arr = np.array(list(zip(primary, *vectors)))
    arr = arr[arr[:, 0].argsort()]
    return [arr[:, i] for i in range(arr.shape[1])]

def sort_by_SDE(OW, PL, Al, Cu, SDE):
    """
    Sorts all input vectors by ascending SDE value.
    Returns: OW_sorted, PL_sorted, Al_sorted, Cu_sorted, SDE_sorted
    """
    sort_idx = np.argsort(SDE)
    return OW[sort_idx], PL[sort_idx], Al[sort_idx], Cu[sort_idx], SDE[sort_idx]

# ============================================================
# 2. Model Functions
# ============================================================

def calc_std_SDE(OW, PL, Al, Cu, A, B, C, D):
    """Polynomial dose-response model."""
    return A*OW + B*PL + C*Al + D*Cu

def model_function(OW, PL, Al, Cu, params):
    """Wrapper for model function selection."""
    A, B, C, D = params[:4]

    return calc_std_SDE(OW, PL, Al, Cu, A, B, C, D)


# ============================================================
# 3. Bias and Optimization Functions
# ============================================================

def calc_bias(known_dose, calc_dose):
    """Absolute fractional bias."""
    return np.abs((known_dose - calc_dose) / known_dose)

# def bias_chi_square(params, OW, PL, Al, Cu, dose_vals):
#     """Objective function: minimize total bias."""
#     calc_doses = model_function(OW, PL, Al, Cu, params)
#     biases = np.abs((dose_vals - calc_doses) / dose_vals)
#     return np.sum(biases)

def bias_chi_square(params, OW, PL, Al, Cu, dose):
    try:
        predicted = model_function(OW, PL, Al, Cu, params)
        residuals = predicted - dose
        chi2 = np.mean(residuals ** 2)
        if not np.isfinite(chi2):
            return 1e12  # Penalize non-finite result
        return chi2
    except Exception as e:
        print(f"⚠️ Error in bias_chi_square: {e}")
        return 1e12

def fit_parameters(OW, PL, Al, Cu, dose, known_type):
    """
    Fit parameters A–D (E–H fixed to 1) with type-specific bounds.
    """
    p0 = np.array([1, 1, 1, 1])             # initial guess

    print(f"Getting bounds for {known_type}")
    bounds = get_fit_bounds(known_type)     # type-specific bounds
    print("Bounds: ")
    print(bounds)

    result = minimize(
        bias_chi_square,
        p0,
        args=(OW, PL, Al, Cu, dose),
        method="L-BFGS-B",   # supports bounds
        bounds=bounds,
        options={"maxiter": 10000, "disp": False},
    )

    if not result.success:
        print(f"⚠️ Optimization warning: {result.message}")
    else:
        print(f"✅ Fit converged: {known_type} parameters within bounds")

    return result.x


# ============================================================
# 4. Data Loading (skip inconsistent headers)
# ============================================================

def load_data(file_path, known_type):
    """
    Loads a space-delimited data file (ignores header) and only keeps rows
    whose radiation_quality() matches known_type.
    Adds the result to GLOBAL_DATA.
    """
    if known_type.startswith("Mixed_"):
        col_names = ["source", "ow", "pl", "al", "cu", "dde", "sde", "ratio"]
    else:
        col_names = ["source", "ow", "pl", "al", "cu", "dde", "sde"]

    file_path = Path(file_path)
    if not file_path.is_file():
        print(f"⚠️ Skipping missing file: {file_path}")
        return pd.DataFrame(columns=["source", "ow", "pl", "al", "cu", "dde", "sde", "rad_quality"])

    try:
        df = pd.read_csv(
            file_path,
            sep='\s+',
            header=None,
            skiprows=1,
            names=col_names,
        )
    except Exception as e:
        print(f"❌ Error reading {file_path.name}: {e}")
        return pd.DataFrame(columns=["source", "ow", "pl", "al", "cu", "dde", "sde", "rad_quality"])

    # Convert numeric columns
    for col in ["ow", "pl", "al", "cu", "dde", "sde"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    # Compute radiation quality
    df["rad_quality"] = df.apply(
        lambda row: radiation_quality(row["ow"], row["pl"], row["al"], row["cu"]),
        axis=1,
    )

    # Filter to desired type
    before = len(df)
    df = df[df["rad_quality"] == ((known_type.replace("_SDE","")).replace("_DDE",""))]
    after = len(df)

    print("----------------------------------")
    print(f"📂 Loading from {file_path.name}")
    print(f"Type: {known_type}")
    print(f"N data Loaded {after}/{before}") 

    # Add to global dataset
    add_data(df)

    return df

def load_multiple(files, known_type):
    """Loads and concatenates multiple files."""
    dfs = [load_data(f, known_type) for f in files if Path(f).is_file()]
    if not dfs:
        raise RuntimeError("No valid input data found.")
    return pd.concat(dfs, ignore_index=True)


# ============================================================
# 5. Plotting
# ============================================================

def plot_fit_comparison(known_dose, nominal, calc_refined, title="Dose Fit Comparison"):
    """Plot known vs calculated dose for two fits."""
    plt.figure(figsize=(7, 7))
    plt.plot(known_dose, calc_refined, "x", label="Refined Fit")
    plt.plot(known_dose, nominal, "o", label="Nominal Fit")
    plt.xlabel("Known Dose")
    plt.ylabel("Calculated Dose")
    plt.title(title)
    plt.legend()
    plt.grid(True, ls="--", alpha=0.6)
    plt.tight_layout()
    plt.show()

def plot_fit_comparison_by_index(known_dose, nominal, calc_refined, title="Dose Fit Comparison"):
    """
    Plot calculated dose for each data point index.
    X-axis: index (i-th measurement)
    Y-axis: calculated and nominal doses.
    """
    # Ensure numpy arrays
    known_dose = np.array(known_dose)
    nominal = np.array(nominal)
    calc_refined = np.array(calc_refined)

    # X-axis: index of each data point (1..N)
    x_vals = np.arange(len(known_dose))

    plt.figure(figsize=(10, 6))
    plt.plot(x_vals, calc_refined, "x", label="Refined Fit")
    plt.plot(x_vals, nominal, "o", label="Nominal Fit")
    plt.xlabel("Data Point Index (i)")
    plt.ylabel("Calculated Dose")
    plt.title(title)
    plt.legend()
    plt.grid(True, ls="--", alpha=0.6)
    plt.tight_layout()
    plt.show()

def plot_fit_comparison_by_index_sorted(known_dose, nominal, calc_refined, title="Dose Fit Comparison (Sorted by Nominal Dose)"):
    """
    Plot calculated dose for each data point, sorted by nominal (or known) dose.
    X-axis: index after sorting (0 = smallest nominal dose)
    Y-axis: dose values.
    """
    # Convert to numpy arrays
    known_dose = np.array(known_dose)
    nominal = np.array(nominal)
    calc_refined = np.array(calc_refined)

    # Sort indices by nominal dose (or known dose)
    sort_idx = np.argsort(nominal)  # or np.argsort(known_dose)
    known_sorted = known_dose[sort_idx]
    nominal_sorted = nominal[sort_idx]
    refined_sorted = calc_refined[sort_idx]

    # X-axis: sorted index
    x_vals = np.arange(len(known_dose))

    plt.figure(figsize=(10, 6))
    plt.plot(x_vals, refined_sorted, "x", label="Refined Fit")
    plt.plot(x_vals, nominal_sorted, "o", label="Nominal (SDE)")
    plt.xlabel("Sorted Data Point Index (by Nominal Dose)")
    plt.ylabel("Dose")
    plt.title(title)
    plt.legend()
    plt.grid(True, ls="--", alpha=0.6)
    plt.tight_layout()
    plt.show()

# ============================================================
# 6. Main Workflow
# ============================================================

def main():
    # --------------------------------------------------------
    # Define your input files (space-delimited)
    # --------------------------------------------------------

    pure_photon_file = "C:\\Users\\jboyd\\OneDrive - Fortive\\jboyd_luxel+\\generated_response_matrices\\All_Pure_Photons\\cpp_All_Pure_Photons_genMeasFromNormRespMat_Nsamples10000_Ncopies1_Dose5000_to_500000_Step250_09m_11h_07_07_2025_v4.csv"
    pure_photon_file1 = "C:\\Users\\jboyd\\OneDrive - Fortive\\jboyd_luxel+\\generated_response_matrices\\All_Pure_Photons\\All_Pure_Photons_genMeasFromNormRespMat_Nsamples20000_Ncopies1_Dose5000_to_500000_Step250_41m_14h_30_07_2025_v4.csv"

    pure_Kr85_file = "C:\\Users\\jboyd\\OneDrive - Fortive\\jboyd_luxel+\\generated_response_matrices\\Kr85\\cpp_Kr85_genMeasFromNormRespMat_Nsamples25000_Ncopies1_Dose250_to_25000_Step25_36m_17h_02_07_2025_v4.csv"
    pure_Sr90_file = "C:\\Users\\jboyd\\OneDrive - Fortive\\jboyd_luxel+\\generated_response_matrices\\Sr90\\cpp_Sr90_genMeasFromNormRespMat_Nsamples200000_Ncopies1_Dose250_to_25000_Step10_36m_17h_02_07_2025_v4.csv"
    pure_DU_file = r"C:\\Users\\jboyd\\OneDrive - Fortive\\jboyd_luxel+\\generated_response_matrices\\DU\\cpp_DU_FrNormRMat_Nsamp25000_Ncop1_Dose250_to_25000_Step25_50m_13h_11_09_2025_v4.csv"

    MixedPhoton_Gen_file = r"C:\\Users\\jboyd\\OneDrive - Fortive\\jboyd_luxel+\\generated_response_matrices\\MixedPhotons\\MixedPhotons_General\\cpp_MixedPhotons_General_mixtures_Nsamples1500_Ncopies1_Dose50_to_5000_by_50_58m_10h_07_07_2025_v4.csv"

    MixedBL_file = r"C:\\Users\\jboyd\\OneDrive - Fortive\\jboyd_luxel+\\generated_response_matrices\\MixedBL\\cpp_MixedBL_mixtures_Samples100000_Ncopies1_mult0_25_Dose30_to_30000_by_25_20m_16h_01_07_2025_v4.csv"

    MixedBH_file = r"C:\\Users\\jboyd\\OneDrive - Fortive\\jboyd_luxel+\\generated_response_matrices\\MixedBH\\cpp_MixedBH_mixtures_Samp3000_Ncopy1_mult0_50_Dose30_to_30000_by_25_04m_14h_13_10_2025_v5.csv"
    MixedBH_file_2 = r"C:\\Users\\jboyd\\OneDrive - Fortive\\jboyd_luxel+\\generated_response_matrices\\MixedBH\\cpp_MixedBH_mixtures_Samp3000_Ncopy1_mult0_50_Dose30_to_30000_by_25_09m_14h_13_10_2025_v5.csv"
    MixedBH_file_3 = r"C:\\Users\\jboyd\\OneDrive - Fortive\\jboyd_luxel+\\generated_response_matrices\\MixedBH\\cpp_MixedBH_mixtures_Samp3000_Ncopy1_mult0_50_Dose30_to_30000_by_25_10m_14h_13_10_2025_v5.csv"
    MixedBH_file_4 = r"C:\\Users\\jboyd\\OneDrive - Fortive\\jboyd_luxel+\\generated_response_matrices\\MixedBH\\cpp_MixedBH_mixtures_Samp3000_Ncopy1_mult0_50_Dose30_to_30000_by_25_11m_14h_13_10_2025_v5.csv"

    MixedBH_file_5 = r"C:\\Users\\jboyd\\OneDrive - Fortive\\jboyd_luxel+\\generated_response_matrices\\MixedBH\\cpp_MixedBH_mixtures_Samp3000_Ncopy1_mult0_50_Dose30_to_30000_by_25_12m_14h_13_10_2025_v5.csv"
    MixedBH_file_6 = r"C:\\Users\\jboyd\\OneDrive - Fortive\\jboyd_luxel+\\generated_response_matrices\\MixedBH\\cpp_MixedBH_mixtures_Samp3000_Ncopy1_mult0_50_Dose30_to_30000_by_25_14m_14h_13_10_2025_v5.csv"
    MixedBH_file_7 = r"C:\\Users\\jboyd\\OneDrive - Fortive\\jboyd_luxel+\\generated_response_matrices\\MixedBH\\cpp_MixedBH_mixtures_Samp3000_Ncopy1_mult0_50_Dose30_to_30000_by_25_15m_14h_13_10_2025_v5.csv"
    MixedBH_file_8 = r"C:\\Users\\jboyd\\OneDrive - Fortive\\jboyd_luxel+\\generated_response_matrices\\MixedBH\\cpp_MixedBH_mixtures_Samp3000_Ncopy1_mult0_50_Dose30_to_30000_by_25_17m_14h_13_10_2025_v5.csv"
 

    files = [
        pure_photon_file,
        pure_photon_file1,
        pure_Kr85_file,
        pure_Sr90_file,
        # pure_DU_file,
        MixedPhoton_Gen_file,
        MixedBL_file,
        MixedBH_file,
        MixedBH_file_2,
        MixedBH_file_3,
        MixedBH_file_4,
        MixedBH_file_5,
        MixedBH_file_6,
        MixedBH_file_7,
        MixedBH_file_8,
    ]



    # --------------------------------------------------------
    # Fit polynomial model (Nominal)
    # --------------------------------------------------------

    # known_source_type = "BL"
    # known_source_type = "BH"

    # known_source_type = "Mixed_BH_SDE"
    # known_source_type = "Mixed_BH_DDE"

    # known_source_type = "Mixed_BL_SDE"
    known_source_type = "Mixed_BL_DDE"   

    # known_source_type = "PH_SDE"
    # known_source_type = "PH_DDE"

    # known_source_type = "PL_PM_SDE"
    # known_source_type = "PL_PM_DDE"

    # --------------------------------------------------------
    # Load and combine all data
    # --------------------------------------------------------

    print("---------------------------------------")
    print(f"Loading data for {known_source_type}")
    print("---------------------------------------")

    df = load_multiple(files, known_source_type)
    print(f"✅ Combined dataset: {len(df)} total entries from {len(files)} files.")

    # Convert numeric fields
    OW, PL, Al, Cu = [df[k].astype(float).values for k in ["ow", "pl", "al", "cu"]]
    DDE = df["dde"].astype(float).values
    SDE = df["sde"].astype(float).values


    if (known_source_type == "BL") or (known_source_type == "BH"):
        dose_type_to_fit = "SDE"
    else:
        if known_source_type.endswith("SDE"):
            dose_type_to_fit = "SDE"
        if known_source_type.endswith("DDE"):
            dose_type_to_fit = "DDE"

    if dose_type_to_fit == "SDE":
        dose_type = SDE
    if dose_type_to_fit == "DDE":
        dose_type = DDE

    calc_nominal = dose_type

    # # --- Sort everything by SDE before fitting ---
    # OW, PL, Al, Cu, SDE = sort_by_SDE(OW, PL, Al, Cu, dose_type)

    # print(f"Number of data loaded: {len(df)}")
    # print(f"Number of data in SDE: {len(SDE)}")

    # --------------------------------------------------------
    # Fit linear model (Refined)
    # --------------------------------------------------------
    print("---------------------------------------")
    print("\n    Starting fit on input data\n")
    print("---------------------------------------")

    if np.any(np.isnan([OW, PL, Al, Cu, dose_type])):
        print("⚠️ NaNs detected in input data")
    if np.any(np.isinf([OW, PL, Al, Cu, dose_type])):
        print("⚠️ Infs detected in input data")

    refined_params = fit_parameters(OW, PL, Al, Cu, dose_type, known_source_type)
    calc_refined = model_function(OW, PL, Al, Cu, refined_params)
    initial_error = 0
    refined_error = 0 

    for calc_OW, nom_OW in zip(calc_refined, dose_type):
        refined_error += math.pow(calc_OW - nom_OW, 2)

    refined_error /= len(dose_type)

    print(f"Lenght of dose_type: {len(dose_type)}")
    print("Refined parameters (A-D):", refined_params)

    print("")
    print("c1 = %0.6f" % (refined_params[0]))
    print("c2 = %0.6f" % (refined_params[1]))
    print("c3 = %0.6f" % (refined_params[2]))
    print("c4 = %0.6f" % (refined_params[3]))
    print("# Refined mean squared error: %0.1f (%0.3f)" % (refined_error, math.sqrt(refined_error)))


    # --------------------------------------------------------
    # Plot comparison
    # --------------------------------------------------------
    plot_fit_comparison_by_index_sorted(SDE, calc_nominal, calc_refined, title="SDE Fit Comparison")


# ============================================================
# 7. Entrypoint
# ============================================================

if __name__ == "__main__":
    main()

