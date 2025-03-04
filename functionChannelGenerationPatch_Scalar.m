function H = functionChannelGenerationPatch_Scalar(rr,rt,wavenumber,M,N,Arv,Arh,Atv,Ath)


%%=============================================================
%The file generates the scalar channel matrices of the UPA-based system with patch antennas of the paper:
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

H = zeros(N,M); 

G = 0;

for n = 1:N
    for  m = 1:M
        
%-Transmitting
r_n = rr(:,n); %Position of n-th receiving patch antenna
r_m = rt(:,m); %Position of m-th transmitting patch antenna

%--Sampling positions of the patch antenna
r_n1 = r_n;

zn2 = r_n(1);
xn2 = r_n(2);
yn2 = r_n(3) + Arv;
r_n2 = [xn2,yn2,zn2]';

zn3 = r_n(1);
xn3 = r_n(2) + Arh;
yn3 = r_n(3) + Arv;
r_n3 = [xn3,yn3,zn3]';

zn4 = r_n(1);
xn4 = r_n(2) + Arh;
yn4 = r_n(3);
r_n4 = [xn4,yn4,zn4]';

zn5 = r_n(1);
xn5 = r_n(2) + Arh/2;
yn5 = r_n(3) + Arv/2;
r_n5 = [xn5,yn5,zn5]';

r_nn = [r_n1,r_n2,r_n3,r_n4,r_n5];


%--m transmitting antenna
r_m1 = r_m;

zm2 = r_m(1);
xm2 = r_m(2);
ym2 = r_m(3) + Atv;
r_m2 = [xm2,ym2,zm2]';

zm3 = r_m(1);
xm3 = r_m(2) + Ath;
ym3 = r_m(3) + Atv;
r_m3 = [xm3,ym3,zm3]';

zm4 = r_m(1);
xm4 = r_m(2) + Ath;
ym4 = r_m(3);
r_m4 = [xm4,ym4,zm4]';

zm5 = r_m(1);
xm5 = r_m(2) + Ath/2;
ym5 = r_m(3) + Atv/2;
r_m5 = [xm5,ym5,zm5]';

r_mm = [r_m1,r_m2,r_m3,r_m4,r_m5];



        for nn = 1:5
            for mm = 1:5

                G_nnmm = functionScalarGreen(r_nn(:,nn),r_mm(:,mm),wavenumber); 
                G = G + G_nnmm/25;

                clear G_nnmm

            end
        end

        H(n,m) = G;

        G = 0;
    end
end
