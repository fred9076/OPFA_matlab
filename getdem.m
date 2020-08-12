function [XX,YY,ZZ,Zi] = getdem(xr,yr,rsize,Nr)
% Interpolation to get DEM for any ROI
% (xr,yr) is the center of the ROI
% rsize is the size of ROI rsize(m) x rsize(m) 
% Nr is the number of pixels Nr x Nr
% Zi is the elevation of the ROI center
% XX, YY are imaging grid
% ZZ is the output DEM of ROI

[DEMa,~] = geotiffread_modified('USGS_13_n40w085.tif'); % read the whole DEM data, modified from the geotiffread function in the map toolbox
DEMa = double(flipud(DEMa)); % latitude value vary from large to small.
Z0 = 244;

NI = size(DEMa,1);

LonMin = -85.0005559296; % from the DEM data 
LonMax = -83.9993522214;
LatMin = 38.9993520234;
LatMax = 40.0005557316;

Lo = linspace(LonMin,LonMax,NI);
La = linspace(LatMin,LatMax,NI);

xxr = linspace(xr - rsize, xr + rsize, Nr);
yyr = linspace(yr - rsize, yr + rsize, Nr);

LonC = -84.098365; % calibrated value
LatC = 39.779721;

OneLat = 111132.92 - 559.82*cosd(2*LatC)+ 1.175*cosd(4*LatC)-0.0023*cosd(6*LatC);
OneLon = 111412.84*cosd(LatC) - 93.5*cosd(3*LatC) + 0.118*cosd(5*LatC);

[XX, YY] = meshgrid(xxr,yyr);

RLoM = XX/OneLon + LonC;
RLaM = YY/OneLat + LatC;


[LoM, LaM] = meshgrid(Lo,La);

DEMs = interp2(LoM,LaM,DEMa,RLoM,RLaM,'spline');

DEMs = DEMs - Z0;
Zi= sum(sum(DEMs(Nr/2:Nr/2+1,Nr/2:Nr/2+1)))/4;

ZZ = DEMs;

return



