function [EDoF_UPA_Patch,Numerator,Denominator] = functionComputerEDoF_UPA_Patch_Scalar_average(rr,rt,M,N,Ath,Atv,Arh,Arv,wavenumber)

%%=============================================================
%The file is used to compute the simulation results for the EDoF for the
%UPA-based system with patch antennas over scalar channels of the paper:
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


for m = 1:M
    for mm = 1:M

        for m1 = 1:5
            for m2 = 1:5
                for n = 1:N
                    for nn = 1:5

                   

                        %--Regions of patch antennas

                        %-Transmitting
                        r_n = rr(:,n); %Position of n-th receiving patch antenna
                        r_m = rt(:,m); %Position of m-th transmitting patch antenna
                        r_mpie = rt(:,mm);

                        %--Sampling positions of the patch antenna
                        r_n1 = r_n;

                        zn2 = r_n(3);
                        xn2 = r_n(1) - 1/2*Arh;
                        yn2 = r_n(2) - 1/2*Arv;
                        r_n2 = [xn2,yn2,zn2]';

                        zn3 = r_n(3);
                        xn3 = r_n(1) - 1/2*Arh;
                        yn3 = r_n(2) + 1/2*Arv;
                        r_n3 = [xn3,yn3,zn3]';

                        zn4 = r_n(3);
                        xn4 = r_n(1) + 1/2*Arh;
                        yn4 = r_n(2) + 1/2*Arv;
                        r_n4 = [xn4,yn4,zn4]';

                        zn5 = r_n(3);
                        xn5 = r_n(1) + Arh/2;
                        yn5 = r_n(2) - Arv/2;
                        r_n5 = [xn5,yn5,zn5]';

                        r_nn = [r_n1,r_n2,r_n3,r_n4,r_n5];


                        %--m transmitting antenna
                        r_m1 = r_m;

                        zm2 = r_m(3);
                        xm2 = r_m(1) - 1/2*Ath;
                        ym2 = r_m(2) - 1/2*Atv;
                        r_m2 = [xm2,ym2,zm2]';

                        zm3 = r_m(3);
                        xm3 = r_m(1) - 1/2*Ath;
                        ym3 = r_m(2) + 1/2*Atv;
                        r_m3 = [xm3,ym3,zm3]';

                        zm4 = r_m(3);
                        xm4 = r_m(1) + 1/2*Ath;
                        ym4 = r_m(2) + 1/2*Atv;
                        r_m4 = [xm4,ym4,zm4]';

                        zm5 = r_m(3);
                        xm5 = r_m(1) + Ath/2;
                        ym5 = r_m(2) - Atv/2;
                        r_m5 = [xm5,ym5,zm5]';

                        r_mm = [r_m1,r_m2,r_m3,r_m4,r_m5];

                        %--mm transmitting antenna
                        r_mm1 = r_mpie;

                        zmm2 = r_mpie(3);
                        xmm2 = r_mpie(1) - 1/2*Ath;
                        ymm2 = r_mpie(2) - 1/2*Atv;
                        r_mm2 = [xmm2,ymm2,zmm2]';

                        zmm3 = r_mpie(3);
                        xmm3 = r_mpie(1) - 1/2*Ath;
                        ymm3 = r_mpie(2) + 1/2*Atv;
                        r_mm3 = [xmm3,ymm3,zmm3]';

                        zmm4 = r_mpie(3);
                        xmm4 = r_mpie(1) + 1/2*Ath;
                        ymm4 = r_mpie(2) + 1/2*Atv;
                        r_mm4 = [xmm4,ymm4,zmm4]';

                        zmm5 = r_mpie(3);
                        xmm5 = r_mpie(1) + Ath/2;
                        ymm5 = r_mpie(2) - Atv/2;
                        r_mm5 = [xmm5,ymm5,zmm5]';

                        r_mmm = [r_mm1,r_mm2,r_mm3,r_mm4,r_mm5];

                    

                        G_Scalar_rt = functionScalarGreen(r_nn(:,nn),r_mm(:,m1),wavenumber);
                        G_Scalar_rtt = functionScalarGreen(r_nn(:,nn),r_mmm(:,m2),wavenumber);

                        Numerator = Numerator + abs(G_Scalar_rt)^2/(5*M);

                        Denominator1 = Denominator1 + G_Scalar_rt'*G_Scalar_rtt;

                    end

                end
                Denominator = Denominator + abs(Denominator1)^2;
                Denominator1 = 0;
            end
        end

    end
end






% Numerator = Numerator*L_tv*L_th*L_rh*L_rv;
% Denominator = (L_tv*L_th*L_rh*L_rv)^2*Denominator;

EDoF_UPA_Patch = Numerator^2/Denominator;
end
