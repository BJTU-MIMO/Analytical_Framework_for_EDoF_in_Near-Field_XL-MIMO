function [EDoF_CAP] = functionComputeEDoF_CAP_Different_Polar(L_tv,L_th,L_rv,L_rh,d,N_Sample_t,N_Sample_r,wavenumber,Polarnumber) 

%%=============================================================
%The file is used to compute the simulation results for the EDoF for the
%2D CAP-based system over channels of different numbers of polarizations of the paper:
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



%---Sampling in the transmitter and receiver
SamplingPoints_t = L_th*rand(N_Sample_t,1) + L_tv*rand(N_Sample_t,1)*1i;
SamplingPoints_tt = L_th*rand(N_Sample_t,1) + L_tv*rand(N_Sample_t,1)*1i;
SamplingPoints_r = L_rh*rand(N_Sample_r,1) + L_rv*rand(N_Sample_r,1)*1i;


% Positions of transmitting antennas (the left-bottom of the transmitter is the origin)
TransmitterPoints_Positions = [real(SamplingPoints_t),imag(SamplingPoints_t),zeros(N_Sample_t,1)]';
TransmitterPoints_Positions_2 = [real(SamplingPoints_tt),imag(SamplingPoints_tt),zeros(N_Sample_t,1)]';
ReceiverPoints_Positions = [real(SamplingPoints_r),imag(SamplingPoints_r),d*ones(N_Sample_r,1)]';


EDoF_CAP = zeros(length(Polarnumber),1);

for pp = 1:length(Polarnumber)

    polarnumber = Polarnumber(pp);
    EDoF_CAP(pp) = functionComputeEDoF_CAP_Particular_Polar(TransmitterPoints_Positions,TransmitterPoints_Positions_2,ReceiverPoints_Positions,N_Sample_t,N_Sample_r,wavenumber,polarnumber); 

end




clear TransmitterPoints_Positions TransmitterPoints_Positions_2 ReceiverPoints_Positions 





