# RZPFA_matlab
A demo of Gotcha Large Scene data orthorectified imaging with DEM via Orthorectified Polar Format Algorithm (OPFA)

## When using this code please cite the relevant papers below: 

    [1] Hu, Ruizhi, et al. "Orthorectified Polar Format Algorithm for Generalized Spotlight SAR Imaging With DEM." IEEE Transactions on Geoscience and Remote Sensing (2020).
    [2] Hu, Ruizhi, et al. "Refocusing and Zoom-In Polar Format Algorithm for Curvilinear Spotlight SAR Imaging on Arbitrary Region of Interest." IEEE Transactions on Geoscience and Remote Sensing 57.10 (2019): 7995-8012
    [3] Hu, Ruizhi, et al. "Curvilinear Video-SAR Persistent Imaging with Distortion Correction Based on Nufft-3." IGARSS 2019-2019 IEEE International Geoscience and Remote Sensing Symposium. IEEE, 2019

## The papers are available on:
https://www.researchgate.net/publication/343356816_Orthorectified_Polar_Format_Algorithm_for_Generalized_Spotlight_SAR_Imaging_With_DEM

https://www.researchgate.net/publication/333733708_Refocusing_and_Zoom-In_Polar_Format_Algorithm_for_Curvilinear_Spotlight_SAR_Imaging_on_Arbitrary_Region_of_Interest

https://www.researchgate.net/publication/332462307_Curvilinear_Video-SAR_Persistent_Imaging_with_Distortion_Correction_Based_on_NuFFT-3

## Contents:
1) DEM_imaging_OPFA_demo.m : main program for DEM imaging
2) getdem.m ： get the DEM for any ROI 
3) geotiffread_modified.m : a modified function to read DEM data to avoid error
4) opfaw.m : Imgaing function with DEM (OPFA in [1])
5) opfawo.m : Imgaing function without DEM (RZPFA-3 in [2])
6) db20.m : a simple function to show image in Decibel
## Results
![image](https://github.com/fred9076/OPFA_matlab/master/result1.png)
![image](https://github.com/fred9076/OPFA_matlab/master/result2.png)
## Dependencies: 
1) Gotcha Large Scene data (Disc1.zip and Disc2.zip) are avaiable on
https://www.sdms.afrl.af.mil./content/public-data/s3_scripts/index.php?file=GotchaLargeSceneData-Disc1.zip
with registration

2) The DEM data "USGS_13_n40w085.tif", avaiable on https://www.sciencebase.gov/catalog/item/5deb329ae4b02caea0f0ea8f

2) NuFFT-3 code in this program is provided on: https://cims.nyu.edu/cmcl/nufft/nufftall-1.3.3.tar.gz 
   
   Mexfile should be built for matlab to call.

    Please also cite their papers below if you are using their codes. 

    * Accelerating the Nonuniform Fast Fourier Transform: (L. Greengard and J.-Y. Lee) SIAM Review 46, 443 (2004).

    * The type 3 nonuniform FFT and its applications: (J.-Y. Lee and L. Greengard) J. Comput. Phys. 206, 1 (2005).

    Other NuFFT-3 schemes is be applicable, such as FINUFFT on: https://finufft.readthedocs.io/en/latest/
    Quite fast, but tend to crush when data or image is large.
    * Barnett, Alexander H., Jeremy Magland, and Ludvig af Klinteberg. "A Parallel Nonuniform Fast Fourier Transform Library Based on an “Exponential of Semicircle" Kernel." SIAM Journal on Scientific Computing 41.5 (2019): C479-C504.

