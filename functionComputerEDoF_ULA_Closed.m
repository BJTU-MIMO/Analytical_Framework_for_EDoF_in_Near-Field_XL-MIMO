function [EDoF_ULA_Closed] = functionComputerEDoF_ULA_Closed(L_r,L_t,M,N,delta_t,delta_r,D,wavenumber) 

%%=============================================================
%The file is used to compute the closed-form results for the EDoF for the
%ULA-based system of the paper:
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

    for n = 1:N

        Numerator = Numerator + 1/(D^2+(-L_r/2 + (n-1)*delta_r + L_t/2 - (m-1)*delta_t)^2);

    end
end

Numerator = Numerator^2*D^4;

%%Calculate the denominator of EDOF
Denominator = 0;
for m1 = 1:M

    sum_m2 = 0;
    for m2 = 1:M

        sum_e = 0;
        for n = 1:N


            x=-((wavenumber*delta_t)/D)*(m1 - m2)*(-L_r/2 + (n-1)*delta_r);

            ex = complex(cos(x),sin(x));
            sum_e = sum_e + ex;
        end
        sum_m2 = sum_m2 + (abs(sum_e)^2);
    end
    Denominator = Denominator + sum_m2;   
end

EDoF_ULA_Closed = Numerator/Denominator;








