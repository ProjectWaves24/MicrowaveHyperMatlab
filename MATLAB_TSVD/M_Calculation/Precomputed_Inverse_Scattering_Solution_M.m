clear all
close all
clc

addpath('./EfielsPlusSMatrix_Baseline')
addpath('./EfielsPlusSMatrix_Baseline/E_Fields')

%%
  % location of every antenna
Ports = importdata('Monopole16 (Discrete Port Locations).txt');

% geometry of the computational domain
Geometrics = load('Geo_04_Monopole16_for_Tomography_Baseline_t0_temp55degc_LowerContrast2MatchingFluid.mat'); %A struct containing {Geo, XX, YY, ZZ}

%Spatial Resolutions [m]
dX = 0.002;
dY = 0.002;
dZ = 0.002;

%Number of Antennas
NAnt = 16;

%Frequency
Freq = 0.915;   %GHz
Fnum = length(Freq);

%Regularization settings
CSTAccdB = -40;     %dB
ForwardAccuracy = 10^(CSTAccdB/20);
NoiseLevel = ForwardAccuracy;       %Noise-free simulation input data

%%
eps0 = 1e-9/36/pi;
mu0 = 4*pi*1e-7;
epsr_b = 24.5;

PortLocations = [Ports(:,2),Ports(:,3),(Ports(:,4)+Ports(:,7))/2];     %for rm and rt (They are coming from the forward solver, so I don't have to get involved in this.)

%Corner points of the imaging area [m]
Xmin = min(Geometrics.XX);                
Xmax = max(Geometrics.XX); 
Ymin = min(Geometrics.YY); 
Ymax = max(Geometrics.YY); 
Zmin = min(Geometrics.ZZ); 
Zmax = max(Geometrics.ZZ); 

Pmin = [Xmin Ymin Zmin];     %[m]
Pmax = [Xmax Ymax Zmax];

X = Geometrics.XX;
Y = Geometrics.YY;
Z = Geometrics.ZZ;

[XM, YM, ZM] = meshgrid(X, Y, Z);
NDeltaO = size(XM,1)*size(XM,2)*size(XM,3);     

%import E-Fields
E_Fields = cell(Fnum,NAnt);      %4-D electric field distributions for each antenna radiation (A*B*C*3).The fourth dimension is for x, y, and z direction. 
E_FileName = @(f,a)['e-field (f=' num2str(f) ') [' num2str(a) ']'];

for Fr=1:Fnum
    for Ant=1:NAnt
        DataTemp = load([E_FileName(Freq(1,Fr),Ant) '.mat']);  %This is a struct. 
        ConHandle = fieldnames(DataTemp);
        E_Fields{Fr,Ant} = DataTemp.(ConHandle{1});         %Conversion from struct to simple matrix.
    end
end

%% Calculating the discritized linear forward Operator (the kernel G)
Gbig = zeros(NAnt^2,NDeltaO);


count = 0;
C = 0.5*1j*pi*Freq*1e9*eps0*epsr_b;         % eq(7)
DeltaV = dX*dY*dZ;

for i=1:NAnt
    Ei = E_Fields{1,i};
    for j=1:NAnt
        count = count + 1;
        Ej = E_Fields{1,j};
        Ebij = sum(Ei.*Ej,4);
        Gbig(count,:) = C*DeltaV*reshape(Ebij, 1, NDeltaO);
    end
end

%% Calculating the precomputed inverse scattering solution 
G = tall(single(Gbig));

[U, S, V] = svd(G,'econ');    %Memory limit and limitation in using tall vectors for svd.
%[U, S, V] = svd(single(G));

UU = gather(U);
SS = gather(S);
VV = gather(V);


save('../TimeSteps_Images/M_16Ant_40x42x26GeoMatrix_Resolution2mm_SVD_Economy.mat','M') % To save the data to be later used in the ImageWraper routine for time-step images.



