# Climate–Smart Agriculture Yield Sensitivity Analysis

### Overview  
This repository contains MATLAB R2023b scripts and datasets for analyzing how **climate, soil, and management practices** (no-till, cover cropping, drainage, irrigation) affect crop yield and drought sensitivity across U.S. agricultural regions (1980s – 2010s).  

The workflow integrates multi-decadal datasets on yield, management, soil, and climate with machine-learning (Random Forest) and statistical analyses to quantify the effects of **climate-smart agriculture (CSA)** practices on the **Agricultural Drought Sensitivity Index (ADSI)** and food production components (Food-Equivalent Units, Energy, Protein).

---

## Folder Structure  

```
├── Fig/                        # Output figures (automatically generated)
├── Shapefile/                  # U.S. county boundary shapefiles (.shp, .shx, .dbf)
├── Data                 		# Data│
├── addVariableArrows.m         # Adds climate / management / soil arrows on maps
├── Analysis_of_Food_Production_Components.m   # FEU, Energy, and Protein analysis (Fig. 1)
├── CSA_ADSI_proc.m             # Main ADSI computation and CSA impact mapping  (Fig. 5)
├── Decadal_changes_in_yield_sensitivity.m     # Decadal yield–SPEI trend mapping (Fig. 2)
├── extract_fixed_effect_Res.m  # Helper for extracting fixed-effect residuals
├── Field_crop_yield_sensitivity_analysis.m    # Site-level analysis (dry vs. wet years) (Fig. 4)
├── plotStackedBarWithTruePies.m # Function for stacked-bar + pie-chart visualization
├── Regional_Conservation_Practice_Mapping.m   # Mapping of No-Till / Drainage / Cover Crop adoption (Fig. 3)
├── RF_Model.m                  # Random Forest model training and validation function
```

---

##  Main Scripts and Their Functions  

| Script | Purpose |
|--------|----------|
| **`Analysis_of_Food_Production_Components.m`** | Quantifies county-level FEU, energy, and protein production and analyzes drivers via Random Forest models. |
| **`Regional_Conservation_Practice_Mapping.m`** | Maps adoption rates and decadal changes in No-Till, Tile Drainage, and Cover Cropping. |
| **`Decadal_changes_in_yield_sensitivity.m`** | Assesses decadal changes in crop yield sensitivity to drought using SPEI correlations. |
| **`Field_crop_yield_sensitivity_analysis.m`** | Evaluates yield responses under combined management systems during dry and wet years. |
| **`CSA_ADSI_proc.m`** | Computes ADSI, integrates CSA management impacts, and visualizes regional/national patterns. |
| **`RF_Model.m`** | Wrapper function for Random Forest training, prediction, and variable-importance extraction. |

---

##  Core Data Files  

| File | Description |
|------|--------------|
| **`CSAdata.mat`** | Master dataset combining yield, climate, soil, and management information. |
| **`Inputdata_All.mat`** | Consolidated inputs for Random Forest modeling. |
| **`speidata.mat`** | County- and month-level SPEI (drought index). |
| **`crop_yield_spei.mat`** | Yield–SPEI sensitivity results by crop. |
| **`Siteyield.mat`** | Site-level yield, management, and climate data for violin-plot analysis. |
| **`UStable.mat`** | County–region lookup table for aggregations. |
| **`USregion.mat`** | Regional polygons for figure overlays. |

---

##  Workflow Summary  

1. **Map Conservation Practices**  
   → `Regional_Conservation_Practice_Mapping.m`  
   Visualize decadal adoption of No-Till, Tile Drainage, and Cover Cropping (1980s–2010s).  

2. **Quantify Food Production Components**  
   → `Analysis_of_Food_Production_Components.m`  
   Compute FEU, energy, and protein yields and their climatic/management drivers.  

3. **Evaluate Yield–Drought Sensitivity**  
   → `Decadal_changes_in_yield_sensitivity.m`  
   Analyze multi-decadal yield sensitivity to drought using SPEI.  

4. **Assess Site-Level Management Effects**  
   → `Field_crop_yield_sensitivity_analysis.m`  
   Compare dry vs. wet-year yield responses under different management combinations.  

5. **Compute and Map ADSI**  
   → `CSA_ADSI_proc.m`  
   Train Random Forest models, simulate CSA effects, and calculate Agricultural Drought Sensitivity Index.  

---

##  System Requirements  

- **MATLAB R2023b**  
- **Required Toolboxes:**  
  - Statistics and Machine Learning Toolbox  
  - Mapping Toolbox  
  - Curve Fitting Toolbox  
- **Optional Utilities:**  
  - `brewermap()` for advanced colormaps (ColorBrewer)  
  - Custom helper functions included:  
    - `addVariableArrows.m`  
    - `plotStackedBarWithTruePies.m`  

---

## Outputs  

All figures are automatically saved to the `Fig/` directory, including:  
- FEU / Energy / Protein production maps  
- Decadal yield sensitivity (SPEI) maps  
- Conservation practice adoption trends  
- ADSI spatial distribution and regional bar charts  



