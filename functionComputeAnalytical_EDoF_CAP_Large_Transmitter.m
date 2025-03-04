function [EDoF_CAP_Closed] = functionComputeAnalytical_EDoF_CAP_Large_Transmitter(L_th,L_rh,L_tv,L_rv,M_s,N_s,D,wavelength)

%%=============================================================
%The file is used to compute the approximate closed-form results in Remark 6 of the paper:
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

%---Sampling to compute exp
x_to_positions = -L_th/2 + rand(M_s,1)*L_th;
x_tu_positions = -L_th/2 + rand(M_s,1)*L_th;
y_to_positions = -L_tv/2 + rand(M_s,1)*L_tv;
y_tu_positions = -L_tv/2 + rand(M_s,1)*L_tv;
y_rkk_positions = -L_rv/2 + rand(N_s,1)*L_rv;
x_rkk_positions = -L_rh/2 + rand(N_s,1)*L_rh;



expp = 0;
exppp = 0;

for o = 1:M_s
    for u = 1:M_s
        for kk = 1:N_s

            d11 = sqrt(D^2 + (x_to_positions(o)-x_rkk_positions(kk))^2 + (y_to_positions(o)-y_rkk_positions(kk))^2);
            d22 = sqrt(D^2 + (x_tu_positions(u)-x_rkk_positions(kk))^2 + (y_tu_positions(u)-y_rkk_positions(kk))^2);

            exppp = exppp + exp(1j*2*pi*d11/wavelength)*exp(-1j*2*pi*d22/wavelength);

        end

        expp = expp + 1/(M_s*N_s)^2*abs(exppp)^2;
        exppp = 0;

    end
end



Numerator_closed = 2*T_function(L_tv,L_rv,L_th/2,D)^2;

Denominator_closed = 2*expp*L_th^2*L_tv^2*L_rv^2/(D^2*(4*D^2+L_th^2)) + expp*L_th*L_tv^2*L_rv^2/D^3*atan(L_th/(2*D))...
    + 8*expp*L_tv^2*L_rv^2/(4*D^2+L_th^2) - 2*expp*L_tv^2*L_rv^2/D^2;


EDoF_CAP_Closed = Numerator_closed/Denominator_closed;





