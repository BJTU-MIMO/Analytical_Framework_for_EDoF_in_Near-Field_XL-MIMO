function [EDoF_CAP,Denominator,Numerator,TransmitterPoints_Positions,ReceiverPoints_Positions] = functionComputeEDoF_CAP(L_tv,L_th,L_rv,L_rh,d,N_Sample_t,N_Sample_r,wavenumber) 

%%=============================================================
%The file is used to compute the simulation results for the EDoF for the
%2D CAP-based system over the dyadic channel of the paper:
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


Numerator = 0;
Denominator = 0;
Denominator1 = 0;




for m = 1:N_Sample_t
    for mm = 1:N_Sample_t

        for p = 1:3
            for q = 1:3

                for n = 1:N_Sample_r

                    G_Dyadic_rt = functionDyadicGreen(ReceiverPoints_Positions(:,n),TransmitterPoints_Positions(:,m),wavenumber);
                    G_Dyadic_rtt = functionDyadicGreen(ReceiverPoints_Positions(:,n),TransmitterPoints_Positions_2(:,mm),wavenumber);

                    Numerator = Numerator + abs(G_Dyadic_rt(p,q))^2/(N_Sample_r*N_Sample_t*N_Sample_t); %Redundant iterations for mm

                    for i = 1:3

                        Denominator1 = Denominator1 + G_Dyadic_rt(i,p)'*G_Dyadic_rtt(i,q)/N_Sample_r;

                    end
                end

                Denominator = Denominator + abs(Denominator1)^2/(N_Sample_t*N_Sample_t);
                Denominator1 = 0;

                clear G_Dyadic_rt G_Dyadic_rtt

            end
        end
    end
end



EDoF_CAP = Numerator^2/Denominator;



