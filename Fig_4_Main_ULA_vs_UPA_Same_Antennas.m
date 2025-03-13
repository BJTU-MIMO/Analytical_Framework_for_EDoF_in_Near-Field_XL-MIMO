%---This file is used to compare the EDoF performance for the UPA-based
%system with the ULA-based system with the same number of antennas

%%=============================================================
%This file generates Figure 4 of the paper:
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


clc
clear all
close all

tic

%----System setup----%
cc = 3*10^8; %Light speed
fc = 30*10^9; %Carrier frequency
wavelength = cc/fc; %Wavelength (meter)
d_total = 10;
wavenumber = 2*pi/wavelength;

alpha_t = [0.5,1,1.5];

%--ULA
M_vv = 2:2:20;
MM = M_vv.^2;



N_vv = 2:2:20;
NN = N_vv.^2;

M_hh = M_vv;
N_hh = N_vv; %Horizontal number of antennas


delta_t = alpha_t*wavelength;
delta_r = alpha_t*wavelength;

delta_rv = alpha_t*wavelength;
delta_rh = alpha_t*wavelength;

delta_tv = alpha_t*wavelength; %Vertical antenna spacing (meter)
delta_th = alpha_t*wavelength; %Horizontal antenna spacing (meter)


%---UPA
EDoF_UPA_total_Scalar = zeros(length(MM),length(alpha_t));

%---ULA
EDoF_ULA_total_Scalar = zeros(length(MM),length(alpha_t));





d = d_total;

for aa = 1:length(alpha_t)

    delta_t = alpha_t(aa)*wavelength;
    delta_r = alpha_t(aa)*wavelength;

    delta_rv = alpha_t(aa)*wavelength;
    delta_rh = alpha_t(aa)*wavelength;

    delta_tv = alpha_t(aa)*wavelength; %Vertical antenna spacing (meter)
    delta_th = alpha_t(aa)*wavelength; %Horizontal antenna spacing (meter)


    for mm = 1:length(M_vv)



        M = MM(mm);
        M_v = M_vv(mm);
        M_h = M_hh(mm);


        N = NN(mm);
        N_v = N_vv(mm);
        N_h = N_hh(mm);

        L_t = delta_t*M;
        L_r = delta_r*N;

        L_tv = delta_tv*M_v;
        L_th = delta_th*M_h;


        L_rv = delta_tv*N_v; %Vertical sige-length (meter)
        L_rh = delta_rh*N_h; %Horizontal sige-length (meter)



        %---UPA

        % Positions of transmitting antennas
        xt_UPA = zeros(M,1);
        yt_UPA = zeros(M,1);
        zt_UPA = zeros(M,1);


        for m = 1:M

            xt_UPA(m) = -L_th/2 + mod(m-1,M_h)*delta_th;
            yt_UPA(m) = -L_tv/2 + floor((m-1)/M_h)*delta_tv;

        end

        rt_UPA = [xt_UPA,yt_UPA,zt_UPA]';


        % Positions of receiving antennas
        xr_UPA = zeros(N,1);
        yr_UPA = zeros(N,1);
        zr_UPA = d*ones(N,1);

        for n = 1:N

            xr_UPA(n) = -L_rh/2 + mod(n-1,N_h)*delta_rh;
            yr_UPA(n) = -L_rv/2 + floor((n-1)/N_h)*delta_rv;

        end
        rr_UPA = [xr_UPA,yr_UPA,zr_UPA]';


        %---Channel generation
        H_UPA_Scalar = functionChannelGeneration_Scalar(rr_UPA,rt_UPA,wavenumber,M,N);

        
        %-Scalar
        R_UPA_Scalar = H_UPA_Scalar*H_UPA_Scalar';
        EDoF_UPA_total_Scalar(mm,aa) = (trace(R_UPA_Scalar)/norm(R_UPA_Scalar,'fro'))^2;



        %---ULA
        xt_ULA = zeros(M,1);
        yt_ULA = zeros(M,1);
        zt_ULA = zeros(M,1);

        for m = 1:M

            yt_ULA(m) = -1/2*L_t + (m-1)*delta_t;

        end
        
        rt_ULA = [xt_ULA,yt_ULA,zt_ULA]';

        xr_ULA = zeros(N,1);
        yr_ULA = zeros(N,1);
        zr_ULA = d*ones(N,1);

        for n = 1:N

            yr_ULA(n) = -L_r/2 + (n-1)*delta_r;

        end
        rr_ULA = [xr_ULA,yr_ULA,zr_ULA]';

        %---Channel generation
        H_ULA_Scalar = functionChannelGeneration_Scalar(rr_ULA,rt_ULA,wavenumber,M,N);


        %-Scalar
        R_ULA_Scalar = H_ULA_Scalar*H_ULA_Scalar';
        EDoF_ULA_total_Scalar(mm,aa) = (trace(R_ULA_Scalar)/norm(R_ULA_Scalar,'fro'))^2;

        clear R_UPA_Scalar H_UPA_Scalar R_ULA_Scalar H_ULA_Scalar

    end
end


figure;
subplot(2,1,1)
hold on; box on;
c1 = plot(N_vv,EDoF_ULA_total_Scalar(:,1),'rx-.','LineWidth',1.7);hold on
c2 = plot(N_vv,EDoF_ULA_total_Scalar(:,2),'b+-','LineWidth',1.7);hold on
c3 = plot(N_vv,EDoF_ULA_total_Scalar(:,3),'k+-','LineWidth',1.7);hold on
grid on
set(gca,'XLim',[2 20]);
set(gca,'xtick',N_vv)
xlabel('Square root of number of total antennas $\sqrt{M}$','Interpreter','Latex');
ylabel('EDoF','Interpreter','Latex');
legend([c1 c2 c3],{'$\Delta _t=0.5\lambda$ ','$\Delta _t=1\lambda$','$\Delta _t=1.5\lambda$'},'Interpreter','Latex','Location','Northwest')
set(gca,'FontSize',12);
title('(a). ULA','Interpreter','Latex','FontSize',12)

subplot(2,1,2)
hold on; box on;
c4 = plot(N_vv,EDoF_UPA_total_Scalar(:,1),'rx-.','LineWidth',1.7);hold on 
c5 = plot(N_vv,EDoF_UPA_total_Scalar(:,2),'b+-','LineWidth',1.7);hold on
c6 = plot(N_vv,EDoF_UPA_total_Scalar(:,3),'k+-','LineWidth',1.7);hold on
grid on
set(gca,'XLim',[2 20]);
set(gca,'YLim',[0.95 2.2]);
xlabel('Square root of number of total antennas $\sqrt{M}$','Interpreter','Latex');
ylabel('EDoF','Interpreter','Latex');
legend([c4 c5 c6],{'$\Delta _t=0.5\lambda$ ','$\Delta _t=1\lambda$','$\Delta _t=1.5\lambda$'},'Interpreter','Latex','Location','Northwest')
set(gca,'FontSize',12);
title('(b). UPA','Interpreter','Latex','FontSize',12)


