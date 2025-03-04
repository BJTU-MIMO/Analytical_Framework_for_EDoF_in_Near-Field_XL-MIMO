%---This file is used to show the impact of the mutual coupling on the EDoF
%performance

%%=============================================================
%This file generates the data applied in Figure 12 of the paper:
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
d_total = 1*wavelength; %Distance between the left-bottom antenna elements of the transmitter and the receiver
wavenumber = 2*pi/wavelength;

%--Transmitter
L_tv = 2*wavelength; %Vertical sige-length (meter)
L_th = 2*wavelength; %Horizontal sige-length (meter)

M_vv = 2:2:18;
M_hh = M_vv;
MM = M_vv.*M_hh; %Number of antennas

%--Receiver
L_rv = 2*wavelength; %Vertical sige-length (meter)
L_rh = 2*wavelength; %Horizontal sige-length (meter)

N_vv = 2:2:18;
N_hh = N_vv; %Horizontal number of antennas
NN = N_vv.*N_hh; %Number of antennas

%---UPA
EDoF_DIS_total_Dyadic = zeros(length(MM),1,length(d_total));
EDoF_DIS_total_Scalar = zeros(length(MM),1,length(d_total));
EDoF_DIS_total_Scalar_Coupling = zeros(length(MM),1,length(d_total));
EDoF_DIS_total_Dyadic_Coupling = zeros(length(MM),1,length(d_total));


% % % % % % % ----启动并行计算
core_number = 4;            %想要调用的处理器个数
parpool('local',core_number);


for dd = 1:length(d_total)

    d = d_total(dd);

    parfor mm = 1:length(M_vv)

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

        % Positions of transmitting antennas
        xt = zeros(M,1);
        yt = zeros(M,1);
        zt = zeros(M,1);


        for m = 1:M

            xt(m) = -L_th/2 + mod(m-1,M_h)*delta_th;
            yt(m) = -L_tv/2 + floor((m-1)/M_h)*delta_tv;

        end

        rt = [xt,yt,zt]';


        % Positions of receiving antennas
        
        xr = zeros(N,1);
        yr = zeros(N,1);
        zr = d*ones(N,1);

        for n = 1:N

            xr(n) = -L_rh/2 + mod(n-1,N_h)*delta_rh;
            yr(n) = -L_rv/2 + floor((n-1)/N_h)*delta_rv;

        end
        rr = [xr,yr,zr]';

        %---Mutual coupling

        [Z_Coupling] = functionMutualCoupling(wavelength,M_v,M_h,delta_tv,delta_th);

        Z_Coupling_Transmit = Z_Coupling;
        Z_Coupling_Receive = Z_Coupling;

        Z_Coupling_Diag = zeros(3*M,3*M);

        for ii = 1:3

            Z_Coupling_Diag((ii-1)*M+1:ii*M,(ii-1)*M+1:ii*M) = Z_Coupling;

        end


        %---Channel generation
        H_Dyadic = functionChannelGeneration(rr,rt,wavenumber,M,N,3);
        H_Dyadic_Coupling = Z_Coupling_Diag*H_Dyadic*Z_Coupling_Diag;


        H_Scalar = functionChannelGeneration_Scalar(rr,rt,wavenumber,M,N);
        H_Scalar_Coupling = Z_Coupling_Receive*H_Scalar*Z_Coupling_Transmit;


        %--Compute EDoF

        %-Dyadic
        R_Dyadic = H_Dyadic*H_Dyadic';
        EDoF_DIS_total_Dyadic(mm,dd) = (trace(R_Dyadic)/norm(R_Dyadic,'fro'))^2;

        R_Dyadic_Coupling = H_Dyadic_Coupling*H_Dyadic_Coupling';
        EDoF_DIS_total_Dyadic_Coupling(mm,dd) = (trace(R_Dyadic_Coupling)/norm(R_Dyadic_Coupling,'fro'))^2;

        
        %-Scalar
        R_Scalar = H_Scalar*H_Scalar';
        EDoF_DIS_total_Scalar(mm,dd) = (trace(R_Scalar)/norm(R_Scalar,'fro'))^2;

        R_Scalar_Coupling = H_Scalar_Coupling*H_Scalar_Coupling';
        EDoF_DIS_total_Scalar_Coupling(mm,dd) = (trace(R_Scalar_Coupling)/norm(R_Scalar_Coupling,'fro'))^2;


        disp([num2str(mm) '-th number of antennas ' num2str(length(M_vv))]);


    end
end


toc


