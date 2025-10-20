library(sf)

# URL of the ArcGIS service
url <- "https://services1.arcgis.com/fBc8EJBxQRMcHlei/arcgis/rest/services/GRSM_HYDROLOGY/FeatureServer/0/query?where=1=1&outFields=*&f=geojson"

# Where to save (change to your preferred path)
save_path <- "/Users/jgradym/Library/CloudStorage/GoogleDrive-jgradym@gmail.com/Shared drives/GRSM_CESU/Maine/Locations/Streams/GRSM_HYDROLOGY.geojson"

# Download and read once
grsm_hydro <- st_read(url)

# Save locally
st_write(grsm_hydro, save_path, delete_dsn = TRUE)

# ---- future sessions ----
# Load directly from disk (no re-download)
grsm_hydro <- st_read(save_path)
