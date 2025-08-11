import numpy as np

# Path configurations
geotiff_dir = "/home/vaithisa/Downscale_Emissions/Downscale_Winter/"
static_path = "/home/vaithisa/GEO4PALM-main/JOBS/Augsburg_3/OUTPUT/O3_small_static"
output_dir = "/home/vaithisa/GEO4PALM-main/JOBS/Augsburg_3/OUTPUT/"

# Emission categories and species
selected_band_prefix = "SumAllSectors"
active_categories = [
   # 'A_PublicPower', 
   # 'B_Industry', 
   # 'C_OtherStationaryComb', 
    #'D_Fugitives',
    #'E_Solvents', 
    'F_RoadTransport', 
    #'G_Shipping', 
    #'H_Aviation',
   # 'I_OffRoad', 
    #'J_Waste', 
   # 'K_AgriLivestock', 
   # 'L_AgriOther',
    #'SumAllSectors'
]
spec_name_str = ('pm10', 'no', 'no2', 'o3')
#spec_name_str = ('n2o', 'nox', 'nmvoc', 'so2', 'co', 'pm10', 'pm2_5', 'nh3', 'pb', 'cd', 'hg', 'as', 'ni', 'bc', 'co2', 'ch4', 'no', 'no2', 'ec', 'oc', 'na', 'so4', 'othmin', 'o3')
spec_name = np.array(spec_name_str, dtype='S64')

# Molar masses (g/mol) for unit conversion; None for particulate matter
molar_mass = {
    'n2o': 44.01,      # N₂O (gas-phase)
    'nox': 46.01,      # Approximated as NO₂ (gas-phase)
    'nmvoc': 44.10,    # Approximated as propane, C₃H₈ (gas-phase)
    'so2': 64.06,      # SO₂ (gas-phase)
    'co': 28.01,       # CO (gas-phase)
    'pm10': None,      # Particulate matter
    'pm2_5': None,     # Particulate matter
    'nh3': 17.03,      # NH₃ (gas-phase)
    'pb': None,        # Lead, treated as particulate matter
    'cd': None,        # Cadmium, treated as particulate matter
    'hg': None,        # Mercury, treated as particulate matter
    'as': None,        # Arsenic, treated as particulate matter
    'ni': None,        # Nickel, treated as particulate matter
    'bc': None,        # Black carbon, treated as particulate matter
    'co2': 44.01,      # CO₂ (gas-phase)
    'ch4': 16.04,      # CH₄ (gas-phase)
    'no': 30.01,       # NO (gas-phase)
    'no2': 46.01,      # NO₂ (gas-phase)
    'ec': None,        # Elemental carbon, treated as particulate matter
    'oc': None,        # Organic carbon, treated as particulate matter
    'na': None,        # Sodium, treated as particulate matter
    'so4': None,       # Sulfate, treated as particulate matter
    'othmin': None,    # Other minerals, treated as particulate matter
    'o3': 48.00        # O₃ (gas-phase)
}

# Model parameters
layer_height = 3  # Thickness of each vertical grid cell in (m) as specified in the PALM configuration file
non_zero_threshold = 1e-20
model_name = "O3_small_3"  #PALM config name
emission_mode = "traffic"
field_length = 64