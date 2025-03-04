function Numerator_closed = functionComputeAnalytical_EDoF_CAP_Numerator(L_th,L_rh,L_tv,L_rv,D)

%%=============================================================
%The file is used to compute the closed-form results for the numerator of the EDoF for the
%2D CAP-based system of the paper:
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


mu0 = (1/(4*pi))^2;

L_H1 = abs(L_th-L_rh)/2;
L_H2 = (L_th+L_rh)/2;

Numerator_closed = 2*L_th*L_rh/max(L_th,L_rh)*T_function(L_tv,L_rv,L_H1,D) + (L_th+L_rh)*T_function(L_tv,L_rv,L_H2,D)...
    - (L_th+L_rh)*T_function(L_tv,L_rv,L_H1,D) - 2*Q_function(L_tv,L_rv,L_H2,D) + 2*Q_function(L_tv,L_rv,L_H1,D);

Numerator_closed = Numerator_closed*mu0;
