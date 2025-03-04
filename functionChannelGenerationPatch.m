function H = functionChannelGenerationPatch(rr,rt,wavenumber,M,N,I,l_rv,l_rh,l_tv,l_th)


%%=============================================================
%The file generates the dyadic channel matrices of the UPA-based system with patch antennas of the paper:
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

H = zeros(I*N,I*M); 

G = zeros(I,I);

for n = 1:N
    for  m = 1:M
        
        r_n = rr(:,n); %Position of n-th receiving patch antenna
        r_m = rt(:,m); %Position of m-th transmitting patch antenna


        %--Sampling positions of the patch antenna
        r_n1 = r_n;
        
        xn2 = r_n(1) - 1/2*l_th;
        yn2 = r_n(2) -1/2*l_tv;
        zn2 = r_n(3);
        r_n2 = [xn2,yn2,zn2]';

        xn3 = r_n(1) - 1/2*l_th;
        yn3 = r_n(2) + 1/2*l_tv;
        zn3 = r_n(3);
        r_n3 = [xn3,yn3,zn3]';

        xn4 = r_n(1) + 1/2*l_th;
        yn4 = r_n(2) + 1/2*l_tv;
        zn4 = r_n(3);
        r_n4 = [xn4,yn4,zn4]';

        xn5 = r_n(1) + 1/2*l_th;
        yn5 = r_n(2) - 1/2*l_tv;
        zn5 = r_n(3);
        r_n5 = [xn5,yn5,zn5]';

        r_nn = [r_n1,r_n2,r_n3,r_n4,r_n5];
        

        r_m1 = r_m;
       
        xm2 = r_m(1) - 1/2*l_rh;
        ym2 = r_m(2) - 1/2*l_rv;
        zm2 = r_m(3);
        r_m2 = [xm2,ym2,zm2]';

        xm3 = r_m(1) - 1/2*l_rh;
        ym3 = r_m(2) + 1/2*l_rv;
        zm3 = r_m(3);
        r_m3 = [xm3,ym3,zm3]';

        xm4 = r_m(1) + 1/2*l_rh;
        ym4 = r_m(2) + 1/2*l_rv;
        zm4 = r_m(3);
        r_m4 = [xm4,ym4,zm4]';

        xm5 = r_m(1) + 1/2*l_rh ;
        ym5 = r_m(2) - 1/2*l_rv;
        zm5 = r_m(3);
        r_m5 = [xm5,ym5,zm5]';

        r_mm = [r_m1,r_m2,r_m3,r_m4,r_m5];


        for nn = 1:5
            for mm = 1:5

                G_nnmm = functionDyadicGreen(r_nn(:,nn),r_mm(:,mm),wavenumber); 

                G = G + G_nnmm/25;

                clear G_nnmm

            end
        end

        for i = 1:I
            for ii = 1:I
                H((i-1)*N+n,(ii-1)*M+m) = G(i,ii);
            end
        end

        G = zeros(I,I);
    end
end
