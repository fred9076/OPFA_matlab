function data = opfaw(data,Xi,Yi,Zi,tol,meth)
% OPFA Algorithm with DEM
% DEM calculated by function getdem
% Xi,Yi,Zi is the location of refocusing point
% tol is the tolerence for nufft;
% (Xi,Yi,Zi) is the refocusing point
% meth can be 'fgg' or 'finufft'
% 
C = 299792458;
ic = round(data.Np/2); % center of Azimuth

Xa = data.AntX'; % curve smoothing may be needed for some trajectory
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
Rfi = sqrt((Xa-Xi).^2 + (Ya-Yi).^2 + (Za-Zi).^2);% w.r.t Refocusing point
El = asind((Za-Zi)./Rfi);
Az = atan2d((Ya-Yi),(Xa-Xi));

Kx = 4*pi/C*data.freq*(cosd(Az).*cosd(El));% Kxi Kyi Non-uniform input wavenumber locations after refocusing
Ky = 4*pi/C*data.freq*(sind(Az).*cosd(El));

Ac = Az(ic);%Center angles
Ec = El(ic);

Kxc = 4*pi/C*data.Fc*(cosd(Ac)*cosd(Ec));%center wavenumber
Kyc = 4*pi/C*data.Fc*(sind(Ac)*cosd(Ec));
disp('Input Calculation finished')

% %%%%%%%%%%%%%%%%%%%  Refocusing   %%%%%%%%%%%%%
% 
ddi = R0A - Rfi; % differential range for refocusing
sr = double(data.phdata.*exp(-1i*4*pi/C*data.freq*ddi)); % Refocusing
disp('Refocusing finished')
%%%%%%%%%%%%%%%%%%%  Non-uniform Output Calculation   %%%%%%%%%%%%%
%  

% %Calculation of the mapping
A = (Xc-data.x_mat)*Xdc + (Yc-data.y_mat)*Ydc + (Zc-data.z_mat)*Zdc;%Mat for image reconstrution
Ai = (Xc-Xi)*Xdc + (Yc-Yi)*Ydc + (Zc-Zi)*Zdc;% For calculating the mapping
Ri = sqrt((Xc-Xi)^2+(Yc-Yi)^2 + (Zc-Zi)^2 );
Rt = sqrt((Xc-data.x_mat).^2+(Yc-data.y_mat).^2 + (Zc-data.z_mat).^2 );%Mat
Di =  Ri.*(Ri-Rt)  ;%Mat
Ei = 2*Ai - A./Rt*Ri - Ai*Rt/Ri ;% Mat
Fi = (Xc-Xi)*Ydc - (Yc- Yi)*Xdc;

% %Non-uniform output grids after mapping
data.Xh = (Ydc*Di-(Yc-Yi)*Ei)/Fi  ; 
data.Yh = (-Xdc*Di+(Xc-Xi)*Ei)/Fi ; 
disp('Output calculation finished')
%%%%%%%%%%___NuFFt3_Imaging____%%%%%%%%%%%%%%
nj = data.K*data.Np;% Number of Input Samples
Kx = Kx(:)-Kxc;
Ky = Ky(:)-Kyc;

xj = Kx;% Input locations
yj = Ky;
iflag = 0;
 sk = -data.Xh(:);% Output locations
 tk = -data.Yh(:);
 nk = length(sk); % Number of output samples
% tol=1e-15; % error tolerance

if strcmp(meth, 'fin')  % Finufft may crush for large data and large image
[fk,~]=finufft2d3(xj,yj,sr(:),iflag,tol,sk,tk); % Call Nufft-3 

elseif strcmp(meth, 'fgg')
[fk,~]=nufft2d3(nj,xj,yj,sr(:),iflag,tol,nk,sk,tk);

end
disp('NuFFT-3 finished')

data.im_opfaw = reshape(fk,[length(data.x_mat(1,:)),length(data.x_mat(1,:))]); % Reshape to image



    
%%%%%%%%%%% Calculating effect scene due to defocus  %%%%%%%%%%%%%%%
% data.rqpew = zeros(size(data.x_mat));
% Bi = (Xc - Xi)* data.Xh +(Yc - Yi)*data.Yh; %Mat
% Ci = Xdc* data.Xh +Ydc*data.Yh; %Mat
% Hi = Xddc* data.Xh +Yddc*data.Yh; %Mat
% G = Xdc^2 + Ydc^2 + Zdc^2 + (Xc - data.x_mat)*Xddc + (Yc - data.y_mat)*Yddc + (Zc-data.z_mat)*Zddc;%Mat
% Gi = Xdc^2 + Ydc^2 + Zdc^2 + (Xc - Xi)*Xddc + (Yc - Yi)*Yddc + (Zc-Zi)*Zddc;
% % Rdt0 = A0./Rt0;%Mat
% Rdi = Ai/Ri;
% Rddt = (G.*(Rt.^2) - A.^2)./(Rt.^3);%Mat
% Rddi = (Gi*(Ri^2) - Ai^2)/(Ri^3);
% data.rqpew = 2*pi/lambda*abs(Rddt - Rddi - Bi*(Ri*Rddi - 2*Rdi^2)/(Ri^3)-(2*Ci*Rdi - Hi*Ri )/(Ri^2));%Mat
   
return