# RZPFA_matlab
A demo of Gotcha Large Scene data imaging via RZPFA3

When using this code please cite the relevant papers below: 

    [1] Hu, Ruizhi, et al. "Refocusing and Zoom-In Polar Format Algorithm for Curvilinear Spotlight SAR Imaging on Arbitrary Region of Interest." IEEE Transactions on Geoscience and Remote Sensing 57.10 (2019): 7995-8012

    [2] Hu, Ruizhi, et al. "Curvilinear Video-SAR Persistent Imaging with Distortion Correction Based on Nufft-3." IGARSS 2019-2019 IEEE International Geoscience and Remote Sensing Symposium. IEEE, 2019

Dependencies: 
1) Gotcha Large Scene data (Disc1.zip and Disc2.zip) are avaiable on
https://www.sdms.afrl.af.mil./content/public-data/s3_scripts/index.php?file=GotchaLargeSceneData-Disc1.zip
with registration

2) NuFFT-3 code in this program is provided on: https://cims.nyu.edu/cmcl/nufft/nufftall-1.3.3.tar.gz 

    Please also cite their papers below if you are using their codes. 

    [1] Accelerating the Nonuniform Fast Fourier Transform: (L. Greengard and J.-Y. Lee) SIAM Review 46, 443 (2004).

    [2] The type 3 nonuniform FFT and its applications: (J.-Y. Lee and L. Greengard) J. Comput. Phys. 206, 1 (2005).

    Other NuFFT-3 schemes may also be applicable, such as FINUFFT on: https://finufft.readthedocs.io/en/latest/

3) db20.m is just a simple function to show image in Decibel
