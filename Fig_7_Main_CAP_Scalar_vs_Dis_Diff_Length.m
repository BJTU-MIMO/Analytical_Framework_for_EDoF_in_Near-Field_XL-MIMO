%---This file is used to show the EDoF performance for the 2D CAP
%plane-based system against the transmitting distance

%%=============================================================
%This file generates the data applied in Figure 7 of the paper:
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
wavenumber = 2*pi/wavelength;

%--Transmitting distance 
d_total = 5:5:35; 



%--Sidelength
Sidelength = 1:0.25:2;

%--Simultaion sampling
N_Sample_t = 150;
N_Sample_r = 150;

%--Analytical sampling
M_s = 150;
N_s = 150;


c = (1/(16*pi^2));

NumberofSampling = 40;


EDoF_CAP_Analytical_total = zeros(NumberofSampling,length(d_total),length(Sidelength));
EDoF_CAP_Simulation_total = zeros(NumberofSampling,length(d_total),length(Sidelength));

% % % % % % % ----Parallel computing
core_number = 2;           
parpool('local',core_number);
% % % Starting parallel pool (parpool) using the 'local' profile ...

for ll = 1:length(Sidelength)

    L_tv = Sidelength(ll); %Vertical sige-length (meter)
    L_th = Sidelength(ll); %Horizontal sige-length (meter)

    %--Receiver
    L_rv = Sidelength(ll); %Vertical sige-length (meter)
    L_rh = Sidelength(ll); %Horizontal sige-length (meter)

    for dd = 1:length(d_total)

        d = d_total(dd);

        parfor i = 1:NumberofSampling

            [EDoF_CAP_Simulation,~,~] = functionComputerEDoF_CAP_Scalar(L_tv,L_th,L_rv,L_rh,d,N_Sample_t,N_Sample_r,wavenumber);
            [Denominator_closed] = functionComputeAnalytical_EDoF_CAP_Denominator(L_th,L_rh,L_tv,L_rv,M_s,N_s,d,wavelength);

            Numerator_closed = functionComputeAnalytical_EDoF_CAP_Numerator(L_th,L_rh,L_tv,L_rv,d);

            EDoF_CAP_Analytical = Numerator_closed^2/Denominator_closed;

            EDoF_CAP_Analytical_total(i,dd,ll) = EDoF_CAP_Analytical;
            EDoF_CAP_Simulation_total(i,dd,ll) = EDoF_CAP_Simulation;
            
        end
    end


end

