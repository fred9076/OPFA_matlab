%%Gotcha LargeScene data imaging using RZPFA-3

% When using this code please cite the relevant paper below: 
% [1] Hu, Ruizhi, et al. "Refocusing and Zoom-In Polar Format Algorithm for Curvilinear Spotlight SAR Imaging on Arbitrary Region of Interest." 
% IEEE Transactions on Geoscience and Remote Sensing 57.10 (2019): 7995-8012
% [2] Hu, Ruizhi, et al. "Curvilinear Video-SAR Persistent Imaging with Distortion Correction Based on Nufft-3." 2019 IGARSS, IEEE International Geoscience and Remote Sensing Symposium. IEEE.


% Dependencies: 
% 1) Gotcha Large Scene data is avaiable on
% https://www.sdms.afrl.af.mil./content/public-data/s3_scripts/index.php?file=GotchaLargeSceneData-Disc1.zip
% and Disc2.zip
% with registration

% 2) NuFFT-3 code in this program is provided on: https://cims.nyu.edu/cmcl/nufft/nufftall-1.3.3.tar.gz 
% Please also cite their papers below if you are using their codes. 
% [1] Accelerating the Nonuniform Fast Fourier Transform: (L. Greengard and J.-Y. Lee) SIAM Review 46, 443 (2004).
% [2] The type 3 nonuniform FFT and its applications: (J.-Y. Lee and L. Greengard) J. Comput. Phys. 206, 1 (2005).
% Other NuFFT-3 schemes may also applicable, such as FINUFFT on: https://finufft.readthedocs.io/en/latest/

% 3) db20.m is just a simple function to show image in Decibel

close all;clear all;clc;
C = 3e8;
fileName = 'subData0'; % Make sure file names are: subData01, subData02,...,subData010
fileSt = 1;% DataStart
fileEn = 1; %DataEnd, 5 in the paper, 10 could be too large for some computers.

Dt = 4096;% The extent of the whole image 4096m*4096m 
Ni = 2048;% pixel number in a sub-image 2048*2048

xc = linspace(-Dt*3/8,Dt*3/8,4);%The focusing centers of 16 patches
yc = linspace(-Dt*3/8,Dt*3/8,4);
[XXC, YYC] = meshgrid(xc,yc);

xx0 = linspace(-Dt/8,Dt/8,Ni); % Range and Azimuth axis in a sub-image.
yy0 = linspace(-Dt/8,Dt/8,Ni);

xwh = linspace(-Dt/2,Dt/2,Ni*4);% Range and Azimuth axis in the whole image.
ywh = linspace(-Dt/2,Dt/2,Ni*4);
[XX0, YY0] = meshgrid(xx0,yy0);
Ind = [ 4 8 12 16 3 7 11 15 2 6 10 14 1 5 9 13]; % Index for each sub-image center, order based on a ¡®flipud¡¯ relation, maybe caused by NuFFT-3 or the reshape process

%Read data and parameters
for kf = 1 : 1   % Number of frames, 2 in the paper
dataA = [];
% Airplane Trajetory
R0A = [];
Xa = [];
Ya = [];
Za = [];
Nt = 0;% Samples in Slow time
for ii = (fileSt:fileEn) + (kf-1)*5 % Read data

    load([fileName,num2str(ii),'.mat']);
    dataA = [dataA subData.phdata];
    R0A = [R0A;subData.R0];
    Xa = [Xa;subData.AntX];
    Ya = [Ya;subData.AntY];
    Za = [Za;subData.AntZ];
    if ii == (fileEn+(kf-1)*5)
        dataK = subData.K; % Number of Samples in range
        deltaF = subData.deltaF; 
        minF = subData.minF;
    end
    Nt = Nt + length(subData.Np);% Total number in slow time
    clear subData;
end

ic = round(Nt/2); % center of Azimuth
data.freq = minF + (0:dataK-1)*deltaF; % k vector
data.freq = data.freq.';
fc = mean(data.freq); 
lambda = C/fc;

Xa = Xa';
Ya = Ya';
Za = Za';
R0A = R0A';

%Apertuer center
Xc = Xa(ic);
Yc = Ya(ic);
Zc = Za(ic);

%First-order derivatives
dt = 2/Nt;
Xdk = diff(Xa)/dt;%
Ydk = diff(Ya)/dt;
Zdk = diff(Za)/dt;
Xdc = Xdk(ic); % Center derivatives
Ydc = Ydk(ic);
Zdc = Zdk(ic);
%Second-order derivative
Xddk = diff(Xdk)/dt;
Yddk = diff(Ydk)/dt;
Zddk = diff(Zdk)/dt;
Xddc = Xddk(ic); % center second-derivatives
Yddc = Yddk(ic);
Zddc = Zddk(ic);
clear Xdk Ydk Zdk Xddk Yddk Zddk; 

ImWhole = zeros(4*Ni,4*Ni); % Image of the whole scene
for np = 1:16 % Number of sub-images, 16 in the paper
    np
    niy = [0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3];% Index for stitch sub-images to Whole image  
    nix = [0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3];
    Xi = XXC(Ind(np));% Refocusing centers
    Yi = YYC(Ind(np));
    
%%%%%%%%%%%%%%%%%%%  Non-uniform Input Calculation   %%%%%%%%%%%%%

Rfi = sqrt((Xa-Xi).^2 + (Ya-Yi).^2 + Za.^2);% w.r.t Refocusing point
El = asind((Za)./Rfi);
Az = atan2d((Ya-Yi),(Xa-Xi));

Kx = 4*pi/C*data.freq*(cosd(Az).*cosd(El));% Kxi Kyi Non-uniform input wavenumber locations after refocusing
Ky = 4*pi/C*data.freq*(sind(Az).*cosd(El));

Ac = Az(ic);%Center angles
Ec = El(ic);

Kxc = 4*pi/C*fc*(cosd(Ac)*cosd(Ec));%center wavenumber
Kyc = 4*pi/C*fc*(sind(Ac)*cosd(Ec));
disp('Input Calculation finished')
   
%%%%%%%%%%%%%%%%%%%  Refocusing   %%%%%%%%%%%%%

ddi = R0A - Rfi; % differential range for refocusing
sr = dataA.*exp(-1i*4*pi/C*data.freq*ddi); % Refocusing
datakk = sr;
disp('Refocusing finished')
clear sr    
%%%%%%%%%%%%%%%%%%%  Non-uniform Output Calculation   %%%%%%%%%%%%%
 
Xz = Xi; % Zoom-in centers, in the paper, Xz is not equal to Xi to show the effects of distortion and defocus.
Yz = Yi;
xx1 = xx0+Xz ; % Axis vector centered on (Xz Yz)
yy1 = yy0+Yz ;
[XX1, YY1] = meshgrid(xx1,yy1);% Ideal pixel grids for a sub-image

%Calculation of the mapping
A0 = (Xc-XX1)*Xdc + (Yc-YY1)*Ydc + Zc*Zdc;%Mat for image reconstrution
Ai0 = (Xc-Xi)*Xdc + (Yc-Yi)*Ydc + Zc*Zdc;% For calculating the mapping
Ri0 = sqrt((Xc-Xi)^2+(Yc-Yi)^2 + Zc^2 );
Rt0 = sqrt((Xc-XX1).^2+(Yc-YY1).^2 + Zc^2 );%Mat
Di0 =  Ri0.*(Ri0-Rt0)  ;%Mat
Ei0 = 2*Ai0 - A0./Rt0*Ri0 - Ai0*Rt0/Ri0 ;
Fi0 = (Xc-Xi)*Ydc - (Yc- Yi)*Xdc;
%Non-uniform output grids after mapping
Xh0 = (Ydc*Di0-(Yc-Yi)*Ei0)/Fi0  ; 
Yh0 = (-Xdc*Di0+(Xc-Xi)*Ei0)/Fi0 ; 

%%%%%%%%%%% Calculating effect scene due to defocus  %%%%%%%%%%%%%%%
% Bi0 = (Xc - Xi)* Xh0 +(Yc - Yi)*Yh0; %Mat
% Ci0 = Xdc* Xh0 +Ydc*Yh0; %Mat
% Hi0 = Xddc* Xh0 +Yddc*Yh0; %Mat
% G0 = Xdc^2 + Ydc^2 + Zdc^2 + (Xc - XX1)*Xddc + (Yc - YY1)*Yddc + Zc*Zddc;%Mat
% Gi = Xdc^2 + Ydc^2 + Zdc^2 + (Xc - Xi)*Xddc + (Yc - Yi)*Yddc + Zc*Zddc;
% Rdt0 = A0./Rt0;%Mat
% Rdi0 = Ai0/Ri0;
% Rddt0 = (G0.*(Rt0.^2) - A0.^2)./(Rt0.^3);%Mat
% Rddi0 = (Gi*(Ri0^2) - Ai0^2)/(Ri0^3);
% Phqi = 2*pi/lambda*abs(Rddt0 - Rddi0 - Bi0*(Ri0*Rddi0 - 2*Rdi0^2)/(Ri0^3)-(2*Ci0*Rdi0 - Hi0*Ri0 )/(Ri0^2));%Mat

% clear A0 Ai0 Ri0 Rt0 Di0 Ei0 Fi0 Bi0 Ci0 Hi0 G0 Gi Rdt0 Rdi0 Rddt0 Rddi0;
disp('Output Calculation finished')

% 
% %%%%%%%%%%___NuFFt3_Imaging____%%%%%%%%%%%%%%
nj = dataK*Nt;% Number of Input Samples
KKKx = Kx(:)-Kxc;
KKKy = Ky(:)-Kyc;

xj = KKKx;% Input locations
yj = KKKy;
iflag = 0;
 sk = -Xh0(:);% Output locations, minus sign for compensation
 tk = -Yh0(:);
 nk = length(sk); % Number of output samples
eps=1e-3; % error tolerance

[fk,ier]=nufft2d3(nj,xj,yj,datakk(:),iflag,eps,nk,sk,tk); % Call Nufft-3 
disp('NuFFT-3 finished')

imnu3 = reshape(fk,[length(Xh0(1,:)),length(Xh0(1,:))]); % Reshape to image

% figure
% colormap('jet')
% imagesc(xx1,yy1,db20(imnu3),[-70,0]); %  Show each sub-image
ImWhole(nix(np)*Ni+1:(nix(np)+1)*Ni,niy(np)*Ni+1:(niy(np)+1)*Ni)=flipud(imnu3); % Store sub-image in the whole image
end
figure % Show whole image
colormap('jet')
imagesc(xwh,ywh,db20(ImWhole),[-70,0]);
end
