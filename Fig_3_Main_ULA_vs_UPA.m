%---This file is used to compare the EDoF performance for the UPA-based
%system with the ULA-based system with the same array aperture

%%=============================================================
%This file generates the data applied in Figure 3 of the paper:
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
d_total = [10,13]; %Distance between the left-bottom antenna elements of the transmitter and the receiver
wavenumber = 2*pi/wavelength;



%--ULA
L_t = 1;
L_r = 1;

M_vv = 1:1:20;
MM = M_vv.^2;



N_vv = 1:1:20;
NN = N_vv.^2;

L_tv = sqrt(2)/2*L_t;
L_th = sqrt(2)/2*L_t;


L_rv = sqrt(2)/2*L_r; %Vertical sige-length (meter)
L_rh = sqrt(2)/2*L_r; %Horizontal sige-length (meter)



M_hh = M_vv;


N_hh = N_vv; %Horizontal number of antennas


%---UPA
EDoF_UPA_total_Scalar = zeros(length(MM),length(d_total));
EDoF_UPA_total_Scalar_Closed = zeros(length(MM),length(d_total));

%---ULA
EDoF_ULA_total_Scalar = zeros(length(MM),length(d_total));
EDoF_ULA_total_Scalar_Closed = zeros(length(MM),length(d_total));


for dd = 1:length(d_total)

    d = d_total(dd);

    for mm = 1:length(M_vv)

        M = MM(mm);
        M_v = M_vv(mm);
        M_h = M_hh(mm);

        delta_tv = L_tv/M_v; %Vertical antenna spacing (meter)
        delta_th = L_th/M_h; %Horizontal antenna spacing (meter)

        N = NN(mm);
        N_v = N_vv(mm);
        N_h = N_hh(mm);

        delta_rv = L_rv/N_v;
        delta_rh = L_rh/N_h;


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
        H_ULA_Scalar = functionChannelGeneration_Scalar(rr_UPA,rt_UPA,wavenumber,M,N);


        %--Compute EDoF
        
        %-Scalar

        R_UPA_Scalar = H_ULA_Scalar*H_ULA_Scalar';
        EDoF_UPA_total_Scalar(mm,dd) = (trace(R_UPA_Scalar)/norm(R_UPA_Scalar,'fro'))^2;

        %-Closed
        [EDoF_UPA_Closed] = functionComputerEDoF_UPA_Closed(L_rh,L_rv,L_th,L_tv,M,N,M_h,N_h,delta_th,delta_tv,delta_rh,delta_rv,d,wavenumber); 
        EDoF_UPA_total_Scalar_Closed(mm,dd) = EDoF_UPA_Closed;


        %---ULA
        
        % Positions of transmitting antennas
        delta_t = L_t/M;
        delta_r = L_r/N;


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
        EDoF_ULA_total_Scalar(mm,dd) = (trace(R_ULA_Scalar)/norm(R_ULA_Scalar,'fro'))^2;

        %-Closed
        [EDoF_ULA_Closed] = functionComputerEDoF_ULA_Closed(L_r,L_t,M,N,delta_t,delta_r,d,wavenumber);
        EDoF_ULA_total_Scalar_Closed(mm,dd) = EDoF_ULA_Closed;



        clear H_UPA_Dyadic R_UPA_Dyadic R_UPA_Scalar H_UPA_Scalar R_ULA_Dyadic H_ULA_Dyadic R_ULA_Scalar H_ULA_Scalar

    end
end












