# Climate–Smart Agriculture Yield Sensitivity Analysis

### **Overview**
This repository contains MATLAB R2023b source code, scripts, and example data for analyzing the effects of **climate-smart agricultural (CSA) practices**—including no-till, cover cropping, tile drainage, and irrigation—on U.S. crop yield sensitivity and drought resilience from 1980 to 2020.  
The workflow integrates multi-decadal yield, management, climate, and soil datasets using **Random Forest (RF)** models and statistical analyses to derive the **Agricultural Drought Sensitivity Index (ADSI)** and related indicators.

---

##  Repository Structure

```
├── Data/                                 # Input datasets (.mat files)
├── Fig/                                  # Output figures generated automatically
│
├── CSA_ADSI_proc.m                       # Main ADSI computation and mapping script
├── Decadal_changes_in_yield_sensitivity.m # Decadal drought sensitivity mapping
├── Field_crop_yield_sensitivity_analysis.m# Site-level dry vs. wet year analysis
├── Regional_Conservation_Practice_Mapping.m # Maps of no-till, drainage, cover crops
│
├── RF_Model.m                            # Random Forest model training function
├── loadAllVars.m                         # Helper for loading all required variables
├── addVariableArrows.m                   # Adds arrow annotations for map visualization
├── plotStackedBarWithTruePies.m          # Combined stacked-bar + pie chart plotter
├── daviolinplot.m                        # Custom violin plot for management comparisons
├── brewmermap.m                          # ColorBrewer colormap function
├── extract_fixed_effect_Res.m            # Utility for extracting fixed-effect results
├── subplot.m                             # Extended subplot layout function
├── Simulation.m                                # Unified entry point running full workflow
├── README.md                             # This file
```

---

##  System Requirements

- **Operating System:** Windows 11
- **Processor:** 12th Gen Intel® Core™ i9-12900 (2.40 GHz, 16 cores / 24 threads)
- **Installed RAM:** 64 GB
- **System Type:** x64-based processor
- **MATLAB Version:** R2023b or later

### **Required MATLAB Toolboxes**
- Statistics and Machine Learning Toolbox
- Mapping Toolbox
- Curve Fitting Toolbox
- Parallel Computing Toolbox (optional for faster training)

### **Performance Notes**
All analyses were tested on this configuration and completed efficiently without high-performance cluster resources.

---

##  Installation Guide

1. Clone or download this repository:
   ```bash
   git clone https://github.com/YakaiWang/ClimateSmart-Ag-Yield-Sensitivity.git
   ```
2. Open MATLAB and set the project folder as your working directory.
3. Ensure the required `.mat` data files are in the `Data/` directory.
4. Run the startup script:
   ```matlab
   Simulation
   ```


---

**Expected output:**
- Spatial maps of ADSI and management effects saved to `/Fig`.

**Runtime:** ~1–5 minutes on a standard workstation.

---

## Instructions for Use

### **1. Run the full CSA workflow**
```matlab
main
```
This runs the sequential workflow:
1. `Regional_Conservation_Practice_Mapping.m`
2. `Decadal_changes_in_yield_sensitivity.m`
3. `CSA_ADSI_proc.m`

### **2. Run site-level drought analysis**
```matlab
Field_crop_yield_sensitivity_analysis
```

---

##  Code Description and Workflow Summary

### **Workflow Logic (Pseudocode)**

```matlab
# 1. Load and preprocess data
loadAllVars
Preprocess climate (gridMET, SPEI)
Preprocess management (NT, CC, DRA, IRR)
Aggregate yield and soil variables by county

# 2. Train Random Forest models for each crop
for each crop
    for each decade
        Train RF_Model with climate, soil, and management predictors
        Evaluate performance (R², RMSE)
        Store predictor importance
    end
end

# 3. Compute Agricultural Drought Sensitivity Index (ADSI)
For each crop:
    Calculate yield–SPEI regression slopes under CSA and baseline
    Aggregate results using crop area weighting

# 4. Visualization
Generate national maps of FEU, Energy, Protein, and ADSI
Add variable arrows and inset regional bar charts
Export figures to /Fig
```

---

##  License

This project is released under the **MIT License**.
You are free to use, modify, and distribute the code with proper citation.

---


