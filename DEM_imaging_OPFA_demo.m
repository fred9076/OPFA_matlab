close all;clear all;clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program performs reading Gotcha large scene data and imaging on DEM via OPFA & RZPFA%
% following fields need to be populated: %
% %
% data.Nfft: Size of the FFT to form the range profile in bp %
% data.deltaF: Step size of frequency data (Hz) %
% data.minF: Vector containing the start frequency of each pulse (Hz) %
 % data.x mat: The x position of each pixel (m) %
 % data.y mat: The y position of each pixel (m) %
 % data.z mat (DEM): The z position of each pixel (m) %
 % data.AntX: The x_a position of the sensor at each pulse (m) %
 % data.AntY: The y_a position of the sensor at each pulse (m) %
 % data.AntZ: The z_a position of the sensor at each pulse (m) %
 % data.R0: The range to scene center (m) %
 % data.phdata: Phase history data (frequency domain) %
 % Fast time in rows, slow time in columns %
 % %
 % The output is: %
 % data.im_opfa: The complex OPFA image value at each pixel
 % data.im_opfawo : The complex RZPFA image value at each pixel

 % Written by Ruizhi Hu, Interdisciplinary Centre for Security, Reliability and Trust, University of Luxembourg %
 % Email: fred9076@gmail.com %
 % Date Released: 12 Aug 2020 %
 % %
 % Part of this code is adapted from the code in .. 
 % Gorham, L.A. and Moore, L.J., "SAR image formation toolbox for %
 % MATLAB," Algorithms for Synthetic Aperture Radar Imagery XVII %
 % 7669, SPIE (2010). %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 %%%%%%%%%%%%%%%%%%%%  ROI parameters   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% location = 'Buildings'; % Imaging area to choose
% location = 'Center';
location = 'Uhill';
switch location
    case 'Flyover' 
        Xi = 667; % Flyover
        Yi = -787;
    case 'Center'
        Xi = 0;
        Yi = 0;
    case 'Uhill'
        Xi = 713;
        Yi = 431;
end

%Choose Imaging Methods%
Opw = 1; % OPFA
Opwo = 1; % RZPFA
 
rsize = 512;% The extent of the image 2rsize m x 2rsize m
Nr = 512;% output pixel size in a ROI Nr x Nr
  
 %%%%%%%%%%%%%%%%%%%   Read Data  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C = 299792458;
fileName = 'subData0'; % Make sure file names are: subData01, subData02,...,subData010
fileSt = 1;% DataStart
fileEn = 2; %DataEnd

%Read data and parameters
data.phdata= [];
% Airplane Trajetory
data.R0 = [];
data.AntX = [];
data.AntY = [];
data.AntZ = [];
data.Np = 0;% Samples in Slow time
for ii = (fileSt:fileEn) 

    load([fileName,num2str(ii),'.mat']);
    data.phdata = [data.phdata subData.phdata];
    data.R0 = [data.R0;subData.R0];
    data.AntX = [data.AntX;subData.AntX];
    data.AntY = [data.AntY;subData.AntY];
    data.AntZ = [data.AntZ;subData.AntZ];
    if ii == fileEn
        data.K = subData.K; % Number of Samples in range
        data.deltaF = subData.deltaF; 
        data.minF = subData.minF;
        data.freq = (0:data.K-1)*data.deltaF + data.minF;
        data.Fc = mean(data.freq);
    end
    data.Np = data.Np + length(subData.Np);% Total number in slow time
    clear subData;
end
data.freq = data.freq.';


    [data.x_mat,data.y_mat,data.z_mat,Zi] = getdem(Xi,Yi,rsize,Nr); % get the DEM of ROI, already subtracted Z0
    data.xaxis = data.x_mat(1,:);
    data.yaxis = data.y_mat(:,1);
    
    figure 
    colormap(jet)
    imagesc(data.xaxis,data.yaxis,data.z_mat);
    axis image
    title('DEM','fontsize',16);
    xlabel('X (m)','fontsize',16);
    ylabel('Y (m)','fontsize',16);
    set(gca,'ydir','normal');
    c=colorbar;
    set(get(c,'label'),'string','Elevation (m)','fontsize',16)     
    saveas(gcf,[location,'_DEM.jpeg']) 

    
    %Imaging via different algorithms%%
    tol = 1e-6;
    
   
   if  Opw == 1
       
       data = opfaw(data,Xi,Yi,Zi,tol,'fgg');
       
       figure; 
        colormap(jet)
        imagesc(data.xaxis,data.yaxis,db20(data.im_opfaw),[-70,0]);
        axis image
         title('OPFA (w/ DEM)','fontsize',16);
        xlabel('X (m)','fontsize',16);
        ylabel('Y (m)','fontsize',16);
        set(gca,'ydir','normal');
        saveas(gcf,[location,'_OPFA_w_DEM.jpeg']) 
   end
     
   if Opwo == 1
       data.z_mat = zeros(size(data.z_mat));
          Zi = 0;
       data = opfawo(data,Xi,Yi,tol,'fgg');
              
        figure; 
        colormap(jet)
        imagesc(data.xaxis,data.yaxis,db20(data.im_opfawo),[-70,0]);
        axis image
         title('OPFA (w/o DEM)','fontsize',16);
        xlabel('X (m)','fontsize',16);
        ylabel('Y (m)','fontsize',16);
        set(gca,'ydir','normal');
        saveas(gcf,[location,'_OPFA_wo_DEM.jpeg']) 
   end
   
    
