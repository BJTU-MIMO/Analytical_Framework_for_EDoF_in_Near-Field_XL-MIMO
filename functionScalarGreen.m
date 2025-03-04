function G_Scalar = functionScalarGreen(rr,rt,wavenumber) 

%%=============================================================
%The file is used to generate scalar Green's function between particular 
% transmitting/receiving points of the paper:
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

d = sqrt((rr-rt)'*(rr-rt)); %Distance between particular transmitting/receiving points

G_Scalar = exp(-1j*wavenumber*d)/(4*pi*d);
end


