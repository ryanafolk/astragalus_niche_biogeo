source("CropRaster.r")
source("CropRaster.r")

CropRaster(filelist=c("/Volumes/monteverdi/Saxifragales_all_layers_30s/BIOCLIM_2.asc"), ShapeFile="saxifragales_simple_seven_regions_newworld.shp", sufix = '_clipped')
CropRaster(filelist=c("/Volumes/monteverdi/Saxifragales_all_layers_30s/BIOCLIM_3.asc"), ShapeFile="saxifragales_simple_seven_regions_newworld.shp", sufix = '_clipped')

CropRaster(filelist=c("/Volumes/monteverdi/Saxifragales_all_layers_30s/BIOCLIM_2.asc"), ShapeFile="saxifragales_simple_seven_regions_oldworld.shp", sufix = '_clipped')
CropRaster(filelist=c("/Volumes/monteverdi/Saxifragales_all_layers_30s/BIOCLIM_3.asc"), ShapeFile="saxifragales_simple_seven_regions_oldworld.shp", sufix = '_clipped')