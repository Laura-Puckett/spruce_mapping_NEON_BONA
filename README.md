# Description
exploratory spruce fractional cover mapping using LVIS, AVIRIS, and NEON data


This goal of this project was to assess the utility of airborne data collected as part of NASA's ABoVE airborne campaign (AVIRIS-NG hyperspectral and LVIS waveform LiDAR data) for estimating percent cover of black spruce forest. 
While this workflow could be expanded over a larger area, this small project was focused on an area near Fairbanks, ALaska, where NEON field measurements of tree inventory are available across a gradient of environmental conditions. 
The NEON field data was manipulated and used as the "ground-truth" values for the proprtional cover of black spruce. 
AVIRIS-NG spectral bands and LVIS relative height metrics were used as predictors in a simple model. 
The main bottleneck in this project was loading and intersecting these large geospatial datasets to produce a training dataset containing plot-level values for all datasets. 

In order to continue this analysis, more sites should be incorporated across a larger area. Hundreds of sites with AVIRIS, LVIS, and ground-truth would be required to train a model for a real mapping effort. Because this didn't seem feasible with field sites that I currently have access to, this project was discontinued and I've shifted to a similar project using Landsat and GEDI instead. 

Dataset References:

NEON Vegetation Structure Data:
National Ecological Observatory Network. 2020. Data Product DP1.10098.001, Woody plant vegetation structure. Provisional data
do wnloaded from http://data.neonscience.org
on November 15, 2020. Battelle, Boulder, CO, USA NEON. 2020.

LVIS LiDAR Data:
Blair, J. B. and M.
Hofton 2020. LVIS Facility L2 Geolocated Surface Elevation and Canopy Height Product, Version 1 . [Indicate subset used]. Boulder, Colorado USA. NASA
National Snow and Ice Data Center Distributed Active Archive Center. doi https://doi.org/10.5067/VP7J20HJQISD . [10/01/

AVIRIS NG Hyperspectral Data:
Miller, C.E., R.O. Green, D.R. Thompson, A.K. Thorpe, M. Eastwood, I.B.
Mccubbin , W. Olson duvall , M. Bernas , C.M. Sarture , S. Nolte, L.M. Rios, M.A. Hernandez, B.D. Bue , and
S.R. Lundeen. 2019. ABoVE: Hyperspectral Imagery from AVIRIS NG, Alaskan and Canadian Arctic, 2017 2018. ORNL DAAC, Oak Ridge, T ennessee,
USA. https://doi.org/10.3334/ORNLDAAC/1569
