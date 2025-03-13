%---This file is used to compare the EDoF performance between 2D CAP
%plane-based XL-MIMO system and that of the 1D CAP plane-based XL-MIMO
%systems

%%=============================================================
%This file generates the data applied in Figure 6 of the paper:
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

%----System setup----%
cc = 3*10^8; %Light speed
fc = 30*10^9; %Carrier frequency
wavelength = cc/fc; %Wavelength (meter)
c = (1/(16*pi^2));
wavenumber = 2*pi/wavelength;

%--Transmitting distance
d_total = [20,30]; 

%--Transmitter
L_tt = 1:0.5:4;
L_rr = 1:0.5:4;

%--Simultaion sampling
N_Sample_t = 100;
N_Sample_r = 100;

%--Analytical sampling
M_s = 100;
N_s = 100;

NumberofSampling = 80;


EDoF_2D_CAP_Analytical_total = zeros(NumberofSampling,length(L_tt),length(d_total));
EDoF_2D_CAP_Simulation_total = zeros(NumberofSampling,length(L_tt),length(d_total));

EDoF_1D_CAP_Analytical_total = zeros(NumberofSampling,length(L_tt),length(d_total));
EDoF_1D_CAP_Simulation_total = zeros(NumberofSampling,length(L_tt),length(d_total));

% % % % % % % ----Parallel computing
core_number = 2;           
parpool('local',core_number);
% % % Starting parallel pool (parpool) using the 'local' profile ...

for dd = 1:length(d_total)

    d = d_total(dd);

    for ll = 1:length(L_tt)

        L_t = L_tt(ll); 
        L_r = L_rr(ll);
        
        L_tv = sqrt(2)/2*L_t; 
        L_th = sqrt(2)/2*L_t;
        
        L_rv = sqrt(2)/2*L_r;
        L_rh = sqrt(2)/2*L_r;


        parfor i = 1:NumberofSampling

            [EDoF_CAP_Simulation,Numerator,Denominator] = functionComputerEDoF_CAP_Scalar(L_tv,L_th,L_rv,L_rh,d,N_Sample_t,N_Sample_r,wavenumber);

            [Denominator_closed] = functionComputeAnalytical_EDoF_CAP_Denominator(L_th,L_rh,L_tv,L_rv,M_s,N_s,d,wavelength);

            Numerator_closed = functionComputeAnalytical_EDoF_CAP_Numerator(L_th,L_rh,L_tv,L_rv,d);

            EDoF_CAP_Analytical = Numerator_closed^2/Denominator_closed;

            EDoF_2D_CAP_Analytical_total(i,ll,dd) = EDoF_CAP_Analytical;
            EDoF_2D_CAP_Simulation_total(i,ll,dd) = EDoF_CAP_Simulation;


            [EDoF_1D_CAP_Scalar,~,~] = functionComputerEDoF_1D_CAP_Scalar(L_t,L_r,d,N_Sample_t,N_Sample_r,wavenumber);
            [EDoF_1D_CAP_closed] = functionComputeAnalytical_EDoF_1D_CAP(L_t,L_r,M_s,N_s,d,wavelength);

            EDoF_1D_CAP_Analytical_total(i,ll,dd) = EDoF_1D_CAP_closed;
            EDoF_1D_CAP_Simulation_total(i,ll,dd) = EDoF_1D_CAP_Scalar;

        end

    end

end


