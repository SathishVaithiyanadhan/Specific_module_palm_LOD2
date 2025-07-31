import numpy as np
from pyproj import Proj, Transformer
import warnings
import os

# Suppress warnings
warnings.filterwarnings("ignore")
warnings.simplefilter("ignore")

print('Reading PALM chemistry configuration')
# Projection configurations
config_proj = "EPSG:25832"  # UTM Zone 32N
default_proj = "EPSG:4326"  # WGS84

# Coordinate transformers
transformer_to_utm = Transformer.from_crs(default_proj, config_proj, always_xy=True)
transformer_to_wgs = Transformer.from_crs(config_proj, default_proj, always_xy=True)

# Path configurations
emis_geotiff_pth = '/home/vaithisa/Downscale_Emissions/Downscale_Winter/'
static_pth = '/home/vaithisa/GEO4PALM-main/JOBS/Augsburg_3/OUTPUT/'
static = 'Augsburg_H'

# Active emission categories (edit these to select sectors)
active_categories = [
    #'A_PublicPower', 
    #'B_Industry', 
    #'C_OtherStationaryComb', 
    #'D_Fugitives',
    #'E_Solvents', 
    'F_RoadTransport', 
    #'G_Shipping', 
    #'H_Aviation',
    #'I_OffRoad', 
    #'J_Waste', 
    #'K_AgriLivestock', 
    #'L_AgriOther',
    #'SumAllSectors'
]
cat_name_str = tuple(active_categories)
cat_name = np.array(cat_name_str, dtype='S64')

# Chemical species configuration (customize as needed)
spec_name_str = ('no2', 'pm10', 'co')  # Example subset; can be modified to any subset of species
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