%---This file is used to show the EDoF performance in the near-field and
%the far-field regions

%%=============================================================
%This file generates Figure 10 of the paper:
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


d_total = 5:1:150;
wavenumber = 2*pi/wavelength;

%--Transmitter
L_tv = [40,60].*wavelength; %Vertical sige-length (meter)
L_th = [40,60].*wavelength; %Horizontal sige-length (meter)

delta_tv = 2*wavelength; %Vertical antenna spacing (meter)
delta_th = 2*wavelength; %Horizontal antenna spacing (meter)

M_vv = floor(L_tv./delta_tv);
M_hh = floor(L_th./delta_th);

MM = M_vv.*M_hh; %Number of antennas

%--Receiver
L_rv = [40,60].*wavelength; %Vertical sige-length (meter)
L_rh = [40,60].*wavelength; %Horizontal sige-length (meter)


delta_rv = 2*wavelength;
delta_rh = 2*wavelength;

N_vv = floor(L_rv./delta_rv);
N_hh = floor(L_rh./delta_rh);
NN = N_vv.*N_hh; %Number of receiving antennas


Aperture_t = sqrt(L_tv.^2+L_th.^2); %Array aperture

d_Rayleigh = 2*(Aperture_t).^2/wavelength; %Rayleigh distance (meter)

EDoF_DIS_Dyadic_total = zeros(length(d_total),length(L_tv));
EDoF_DIS_Scalar_total = zeros(length(d_total),length(L_tv));


for dd = 1:length(d_total)

    d = d_total(dd);

    for mm = 1:length(L_tv)

        M = MM(mm);
        M_v = M_vv(mm);
        M_h = M_hh(mm);


        N = NN(mm);
        N_v = N_vv(mm);
        N_h = N_hh(mm);

        % Positions of transmitting antennas
        xt = zeros(M,1);
        yt = zeros(M,1);
        zt = zeros(M,1);

        for m = 1:M

            yt(m) = mod(m-1,M_h)*delta_th;
            zt(m) = floor((m-1)/M_h)*delta_tv;

        end

        rt = [xt,yt,zt]';


        % Positions of receiving antennas
        xr = d*ones(N,1);
        yr = zeros(N,1);
        zr = zeros(N,1);

        for n = 1:N

            yr(n) = mod(n-1,N_h)*delta_rh;
            zr(n) = floor((n-1)/N_h)*delta_rv;

        end
        rr = [xr,yr,zr]';

        %---Channel generation


        H_Dyadic = functionChannelGeneration(rr,rt,wavenumber,M,N,3);

        H_Scalar = functionChannelGeneration_Scalar(rr,rt,wavenumber,M,N);


        %--Compute EDoF

        %-Dyadic

        R_Dyadic = H_Dyadic*H_Dyadic';

        EDoF_DIS_Dyadic_total(dd,mm) = (trace(R_Dyadic)/norm(R_Dyadic,'fro'))^2;

        %-Scalar

        R_Scalar = H_Scalar*H_Scalar';
        EDoF_DIS_Scalar_total(dd,mm) = (trace(R_Scalar)/norm(R_Scalar,'fro'))^2;

        disp([num2str(dd) '-th distance of ' num2str(length(d_total))]);
        disp([num2str(mm) '-th aperture of ' num2str(length(L_tv))]);




        clear H_Dyadic R_Dyadic R_Scalar H_Scalar

        1;

    end
end



figure;
hold on; box on;
c1 = plot(d_total,EDoF_DIS_Dyadic_total(:,1),'r-','LineWidth',1.7);hold on 
c2 = plot(d_total,EDoF_DIS_Dyadic_total(:,2),'b-','LineWidth',1.7);hold on 
c3 = plot(d_total,EDoF_DIS_Scalar_total(:,1),'r-.','LineWidth',1.7);hold on
c4 = plot(d_total,EDoF_DIS_Scalar_total(:,2),'b-.','LineWidth',1.7);hold on
line([d_Rayleigh(1) d_Rayleigh(1)],[0 35],'linestyle','--','Color','r','LineWidth',1.7);
c5 = line([d_Rayleigh(2) d_Rayleigh(2)],[0 35],'linestyle','--','Color','b','LineWidth',1.7);
grid on
set(gca,'XLim',[8 150]);
xlabel('Transmitting distance $D \ [m] $','Interpreter','Latex');
ylabel('EDoF','Interpreter','Latex');
legend([c2 c1 c4 c3 c5],{'Dyadic, $L=60\lambda$','Dyadic, $L=40\lambda$','Scalar, $L=60\lambda$','Scalar, $L=40\lambda$','Rayleigh distance'},'Interpreter','Latex','Location','Northwest')
set(gca,'FontSize',12);





