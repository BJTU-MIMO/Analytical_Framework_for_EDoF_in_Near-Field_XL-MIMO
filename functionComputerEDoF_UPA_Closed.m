function [EDoF_UPA_Closed] = functionComputerEDoF_UPA_Closed(L_rh,L_rv,L_th,L_tv,M,N,M_h,N_h,delta_th,delta_tv,delta_rh,delta_rv,D,wavenumber) 

%%=============================================================
%The file is used to compute the closed-form results for the EDoF for the
%UPA-based system of the paper:
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


%%Calculate the molecules of discrete EDOF
Numerator=0;
for m = 1:M
    i = mod(m-1,M_h);
    j = floor((m-1)/M_h);
    for n = 1:N
        t = mod(n-1,N_h);
        k = floor((n-1)/N_h);
        Numerator = Numerator + 1/(D^2+(-L_rh/2 + t*delta_rh + L_th/2 - i*delta_th)^2 + (-L_rv/2 + k*delta_rv + L_tv/2 - j*delta_tv)^2);
    end
end

Numerator = Numerator^2*D^4;

%%Calculate the denominator of EDOF
Denominator = 0;
for m1 = 1:M
    i1 = mod(m1-1,M_h);
    j1 = floor((m1-1)/M_h);
    sum_m2 = 0;
    for m2 = 1:M
        i2 = mod(m2-1,M_h);
        j2 = floor((m2-1)/M_h);
        sum_e = 0;
        for n = 1:N
            t1 = mod(n-1,N_h);
            k1 = floor((n-1)/N_h);
            x=-((wavenumber*delta_th)/D)*(i1-i2)*(-L_rh/2 + t1*delta_rh)-((wavenumber*delta_tv)/D)*(j1-j2)*(-L_rv/2 + k1*delta_rv);

            ex = complex(cos(x),sin(x));
            sum_e = sum_e + ex;
        end
        sum_m2 = sum_m2 + (abs(sum_e)^2);
    end
    Denominator = Denominator + sum_m2;   
end

EDoF_UPA_Closed = Numerator/Denominator;








