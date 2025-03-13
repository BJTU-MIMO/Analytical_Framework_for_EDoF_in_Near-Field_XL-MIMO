%---This file is used to compare the EDoF performance for the UPA-based
%system with patch antennas and the UPA-based system with infinitely thin
%dipoles

%%=============================================================
%This file generates the data applied in Figure 5 of the paper:
%
%Zhe Wang, Jiayi Zhang, Wenhui Yi, Huahua Xiao, Hongyang Du, Dusit Niyato,
%Bo Ai, and Derrick Wing Kwan Ng, "Analytical Framework for Effective Degrees of Freedom in Near-Field XL-MIMO,"
%IEEE Transactions on Wireless Communications, to appear, 2025, %doi: 10.1109/TWC.2025.3531418.
%
%Download article: https://arxiv.org/abs/2401.15280 or https://ieeexplore.ieee.org/document/10856805
%
%License: This code is licensed under the GPLv2 license. If you in any way
%use this code for research that results in publications, please cite our
%paper as described above.
%%=============================================================

clc
clear all
close all

tic

%----System setup----%
cc = 3*10^8; %Light speed
fc = 30*10^9; %Carrier frequency
wavelength = cc/fc; %Wavelength (meter)
d = 10*wavelength;
wavenumber = 2*pi/wavelength;
alpha = 1/2;


N_antenna = 2:2:14;

%--Transmitter
L_tv = 10*wavelength; %Vertical sige-length (meter)
L_th = 10*wavelength; %Horizontal sige-length (meter)


Aperture_t = sqrt(L_tv^2+L_th^2); %Array aperture



Atv = alpha*wavelength; %Physical size of each patch antenna element
Ath = alpha*wavelength;

%--Receiver
L_rv = 10*wavelength; %Vertical sige-length (meter)
L_rh = 10*wavelength; %Horizontal sige-length (meter)

Aperture_r = sqrt(L_rv^2+L_rh^2);


Arv = alpha*wavelength;
Arh = alpha*wavelength;

EDoF_Point_Scalar_Total = zeros(length(N_antenna),1);
EDoF_Point_Dyadic_Total = zeros(length(N_antenna),1);
EDoF_Patch_Scalar_Total = zeros(length(N_antenna),1);
EDoF_Patch_Dyadic_Total = zeros(length(N_antenna),1);


% % % % % % % ----Parallel computing
core_number = 2;           
parpool('local',core_number);
% % % Starting parallel pool (parpool) using the 'local' profile ...


parfor cc = 1:length(N_antenna)

    M_v = N_antenna(cc); %Vertical number of antennas
    M_h = N_antenna(cc); %Horizontal number of antennas
    M = M_v*M_h; %Number of antennas

    N_v = N_antenna(cc); %Vertical number of antennas
    N_h = N_antenna(cc); %Horizontal number of antennas
    N = N_v*N_h; %Number of antennas


    delta_rv = L_rv/N_v; %Vertical antenna spacing (meter)
    delta_rh = L_rh/N_h; %Horizontal antenna spacing (meter)

    delta_tv = L_tv/M_v; %Vertical antenna spacing (meter)
    delta_th = L_th/M_h; %Horizontal antenna spacing (meter)



    % Positions of transmitting antennas (the left-bottom is the origin)
    zt = zeros(M,1);
    yt = zeros(M,1);
    xt = zeros(M,1);

    for m = 1:M

        xt(m) = -L_th/2 + mod(m-1,M_h)*delta_th;
        yt(m) = -L_tv/2 + floor((m-1)/M_h)*delta_tv;

    end

    rt = [xt,yt,zt]';


    % Positions of receiving antennas
    xr = zeros(N,1);
    yr = zeros(N,1);
    zr = d*ones(N,1);

    for n = 1:N

        xr(n) = -L_rh/2 + mod(n-1,N_h)*delta_rh;
        yr(n) = -L_rv/2 + floor((n-1)/N_h)*delta_rv;

    end
    rr = [xr,yr,zr]';


    [EDoF_UPA_Patch_Scalar,~,~] = functionComputerEDoF_UPA_Patch_Scalar_average(rr,rt,M,N,Ath,Atv,Arh,Arv,wavenumber);
    [EDoF_UPA_Patch_Dyadic,~,~] = functionComputerEDoF_UPA_Patch_Dyadic_average(rr,rt,M,N,Ath,Atv,Arh,Arv,wavenumber);


    H_Point_Dyadic = functionChannelGeneration(rr,rt,wavenumber,M,N,3);
    H_Point_Scalar = functionChannelGeneration_Scalar(rr,rt,wavenumber,M,N);

    %--Compute EDoF

    R_Point_Dyadic = H_Point_Dyadic*H_Point_Dyadic';
    EDoF_Point_Dyadic = (trace(R_Point_Dyadic)/norm(R_Point_Dyadic,'fro'))^2;



    R_Point_Scalar = H_Point_Scalar*H_Point_Scalar';
    EDoF_Point_Scalar = (trace(R_Point_Scalar)/norm(R_Point_Scalar,'fro'))^2;


    EDoF_Point_Scalar_Total(cc) = EDoF_Point_Scalar;
    EDoF_Point_Dyadic_Total(cc) = EDoF_Point_Dyadic;
    EDoF_Patch_Scalar_Total(cc) = EDoF_UPA_Patch_Scalar;
    EDoF_Patch_Dyadic_Total(cc) = EDoF_UPA_Patch_Dyadic;
1
end
toc

