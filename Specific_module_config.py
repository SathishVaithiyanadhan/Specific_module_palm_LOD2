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
    'n2o': 44.01,      
    'nox': 46.01,      
    'nmvoc': 44.10,    # Approximated as propane C3H8
    'so2': 64.06,      
    'co': 28.01,       
    'pm10': None,      
    'pm2_5': None,     
    'nh3': 17.03,     
    'pb': None,        
    'cd': None,        
    'hg': None,        
    'as': None,        
    'ni': None,        
    'bc': None,        
    'co2': 44.01,      
    'ch4': 16.04,      
    'no': 30.01,      
    'no2': 46.01,      
    'ec': None,        
    'oc': None,       
    'na': None,       
    'so4': None,       
    'othmin': None,    
    'o3': 48.00       
}

# Model parameters
layer_height = 3  # Thickness of each vertical grid cell in (m) as specified in the PALM configuration file
non_zero_threshold = 1e-20
model_name = "O3_small_3"  #PALM config name
emission_mode = "traffic"
field_length = 64