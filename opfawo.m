function data = opfawo(data,Xi,Yi,tol,meth)
% OPFA Algorithm without DEM, ie Zi = 0, data.z_mat = 0
% DEM calculated by function getdem
% Xi,Yi,Zi is the location of refocusing point
% tol is the tolerence for Nufft
% (Xi,Yi,Zi) is the refocusing point
% meth can be 'fgg' or 'finufft'
C = 299792458;
ic = round(data.Np/2); % center of Azimuth
Xa = data.AntX';
Ya = data.AntY';
Za = data.AntZ';
R0A = data.R0';
lambda = C/data.Fc;


% 
%Apertuer center
Xc = Xa(ic);
Yc = Ya(ic);
Zc = Za(ic);
% 
%First-order derivatives
dt = 2/data.Np;
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
% clear Xdk Ydk Zdk Xddk Yddk Zddk; 


%%%%%%%%%%%%%%%%%%%  Non-uniform Input Calculation   %%%%%%%%%%%%%
% 
Rfi0 = sqrt((Xa-Xi).^2 + (Ya-Yi).^2 + Za.^2);% w.r.t Refocusing point
El = asind(Za./Rfi0);
Az = atan2d((Ya-Yi),(Xa-Xi));

Kx0 = 4*pi/C*data.freq*(cosd(Az).*cosd(El));% Kxi Kyi Non-uniform input wavenumber locations after refocusing
Ky0 = 4*pi/C*data.freq*(sind(Az).*cosd(El));

Ac = Az(ic);%Center angles
Ec = El(ic);

Kxc0 = 4*pi/C*data.Fc*(cosd(Ac)*cosd(Ec));%center wavenumber
Kyc0 = 4*pi/C*data.Fc*(sind(Ac)*cosd(Ec));
disp('Input Calculation finished')

% %%%%%%%%%%%%%%%%%%%  Refocusing   %%%%%%%%%%%%%
% 
ddi0 = R0A - Rfi0; % differential range for refocusing
srwo = double(data.phdata.*exp(-1i*4*pi/C*data.freq*ddi0)); % Refocusing
disp('Refocusing finished')
%%%%%%%%%%%%%%%%%%%  Non-uniform Output Calculation   %%%%%%%%%%%%%
%  

% %Calculation of the mapping
A0 = (Xc-data.x_mat)*Xdc + (Yc-data.y_mat)*Ydc + Zc*Zdc;%Mat for image reconstrution
Ai0 = (Xc-Xi)*Xdc + (Yc-Yi)*Ydc + Zc*Zdc;% For calculating the mapping
Ri0 = sqrt((Xc-Xi)^2+(Yc-Yi)^2 + Zc^2 );
Rt0 = sqrt((Xc-data.x_mat).^2+(Yc-data.y_mat).^2 + Zc^2 );%Mat
Di0 =  Ri0.*(Ri0-Rt0)  ;%Mat
Ei0 = 2*Ai0 - A0./Rt0*Ri0 - Ai0*Rt0/Ri0 ;% Mat
Fi0 = (Xc-Xi)*Ydc - (Yc- Yi)*Xdc;

% %Non-uniform output grids after mapping
data.Xh0 = (Ydc*Di0-(Yc-Yi)*Ei0)/Fi0  ; 
data.Yh0 = (-Xdc*Di0+(Xc-Xi)*Ei0)/Fi0 ; 
disp('Output calculation finished')
%%%%%%%%%%___NuFFt3_Imaging____%%%%%%%%%%%%%%
nj = data.K*data.Np;% Number of Input Samples
Kx0 = Kx0(:)-Kxc0;
Ky0 = Ky0(:)-Kyc0;

xj = Kx0;% Input locations
yj = Ky0;
iflag = 0;
 sk = -data.Xh0(:);% Output locations, minus sign for distortioncompensation
 tk = -data.Yh0(:);
 nk = length(sk); % Number of output samples
if strcmp(meth, 'fin')
[fkwo,~]=finufft2d3(xj,yj,srwo(:),iflag,tol,sk,tk); % Call Nufft-3 
elseif strcmp(meth, 'fgg')
[fkwo,~]=nufft2d3(nj,xj,yj,srwo(:),iflag,tol,nk,sk,tk);

end
disp('NuFFT-3 finished')

data.im_opfawo = reshape(fkwo,[length(data.x_mat(1,:)),length(data.x_mat(1,:))]); % Reshape to image
   %%%%%%%%%%% Calculating effect scene due to defocus without DEM %%%%%%%%%%%%%%%
% data.rqpewo = zeros(size(data.x_mat));
% Bi0 = (Xc - Xi)* data.Xh0 +(Yc - Yi)*data.Yh0; %Mat
% Ci0 = Xdc* data.Xh0 +Ydc*data.Yh0; %Mat
% Hi0 = Xddc* data.Xh0 +Yddc*data.Yh0; %Mat
% G0 = Xdc^2 + Ydc^2 + Zdc^2 + (Xc - data.x_mat)*Xddc + (Yc - data.y_mat)*Yddc + Zc*Zddc;%Mat
% Gi0 = Xdc^2 + Ydc^2 + Zdc^2 + (Xc - Xi)*Xddc + (Yc - Yi)*Yddc + Zc*Zddc;
% 
% Rdi0 = Ai0/Ri0;
% Rddt0 = (G0.*(Rt0.^2) - A0.^2)./(Rt0.^3);%Mat
% Rddi0 = (Gi0*(Ri0^2) - Ai0^2)/(Ri0^3);
% data.rqpewo = 2*pi/lambda*abs(Rddt0 - Rddi0 - Bi0*(Ri0*Rddi0 - 2*Rdi0^2)/(Ri0^3)-(2*Ci0*Rdi0 - Hi0*Ri0 )/(Ri0^2));%Mat
     
   
return