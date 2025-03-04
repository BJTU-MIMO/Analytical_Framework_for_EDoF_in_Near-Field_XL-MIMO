function [EDoF_CAP] = functionComputeEDoF_CAP_Particular_Polar(TransmitterPoints_Positions,TransmitterPoints_Positions_2,ReceiverPoints_Positions,N_Sample_t,N_Sample_r,wavenumber,polarnumber) 

%%=============================================================
%The file is used to compute the simulation results for the EDoF for the
%2D CAP-based system over channels of a particular number of polarizations of the paper:
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

Numerator = 0;
Denominator = 0;
Denominator1 = 0;



for m = 1:N_Sample_t
    for mm = 1:N_Sample_t

        for p = 1:polarnumber
            for q = 1:polarnumber

                for n = 1:N_Sample_r

                    G_Dyadic_rt = functionDyadicGreen_Arbitrary_Polar(ReceiverPoints_Positions(:,n),TransmitterPoints_Positions(:,m),wavenumber,polarnumber);
                    G_Dyadic_rtt = functionDyadicGreen_Arbitrary_Polar(ReceiverPoints_Positions(:,n),TransmitterPoints_Positions_2(:,mm),wavenumber,polarnumber);

                    Numerator = Numerator + abs(G_Dyadic_rt(p,q))^2/(N_Sample_r*N_Sample_t*N_Sample_t); %Redundant iterations for mm

                    for i = 1:polarnumber

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

clear TransmitterPoints_Positions TransmitterPoints_Positions_2 ReceiverPoints_Positions 





