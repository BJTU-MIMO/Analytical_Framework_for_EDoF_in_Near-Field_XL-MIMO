function G = functionDyadicGreen_Arbitrary_Polar(rr,rt,wavenumber,polarnumber) 

%%=============================================================
%The file is used to generate Green's function with the particular number of polarizations between particular 
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
ar = (rr-rt)./d; %Unit direction vector
I = eye(3,3);

G_Scalar = exp(-1j*wavenumber*d)/(4*pi*d);
G1 = (1+1j/(wavenumber*d)-1/(wavenumber*d)^2);
G2 = (3/(wavenumber*d)^2-3j/(wavenumber*d)-1);


G_All_Polar = G1*I*G_Scalar + G2*(ar*ar')*G_Scalar; %Generate the dyadic Green's function channel


G = G_All_Polar(1:polarnumber,1:polarnumber);



