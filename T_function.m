function output = T_function(L_tv,L_rv,x,D)

%%=============================================================
%The file generates the T_function of the paper:
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


mu1 = (L_tv-L_rv)^2 + 4*D^2;
mu2 = (L_tv+L_rv)^2 + 4*D^2;

output = 2*L_tv*L_rv/D*atan(x/D) + x*log(mu1+4*x^2) - x*log(mu2+4*x^2)...
    + sqrt(mu1)*atan(2*x/sqrt(mu1)) - sqrt(mu2)*atan(2*x/sqrt(mu2));


end

