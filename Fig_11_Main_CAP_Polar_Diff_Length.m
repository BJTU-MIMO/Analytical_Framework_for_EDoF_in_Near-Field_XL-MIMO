%---This file is used to compare the EDoF performance for the 2D CAP
%plane-based system over different numbers of channel polarizations

%%=============================================================
%This file generates the data applied in Figure 11 of the paper:
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
d = 6*wavelength;
wavenumber = 2*pi/wavelength;
Polarnumber = 1:1:3;

Sidelength = (2:2:12)*wavelength;
N_Sample_t = 150;
N_Sample_r = 150;

NumberofSampling = 40;

EDoF_CAP_Total = zeros(length(Polarnumber),NumberofSampling,length(Sidelength));

% % % % % % % ----Parallel computing
core_number = 2;           
parpool('local',core_number);
% % % Starting parallel pool (parpool) using the 'local' profile ...

for ll = 1:length(Sidelength)

    L_tv = Sidelength(ll);
    L_th = Sidelength(ll);

    L_rv = Sidelength(ll);
    L_rh = Sidelength(ll);
    
    parfor i = 1:NumberofSampling
        
        [EDoF_CAP_All_Polar] = functionComputeEDoF_CAP_Different_Polar(L_tv,L_th,L_rv,L_rh,d,N_Sample_t,N_Sample_r,wavenumber,Polarnumber);
        EDoF_CAP_Total(:,i,ll) = EDoF_CAP_All_Polar;

    end


end


toc















