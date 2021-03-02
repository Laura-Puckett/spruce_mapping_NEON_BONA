library(caTools); library(raster)
setwd('//nau.froot.nau.edu/cirrus/scratch/lkp58/datasets/AVIRIS/sites/BONA/')

# filename='./ang20180723t193341_rfl_v2r2/ang20180723t193341_corr_v2r2_img'
filename='./2_gdalwarp_to_tiff/ang20180723t193341_corr_v2r2_img.tif'
img = brick(filename) 

