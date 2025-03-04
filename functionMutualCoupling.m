function [Z] = functionMutualCoupling(wavelength,N_x,N_y,AntennaSpacing_x,AntennaSpacing_y)

%%=============================================================
%The file generates the mutual coupling matrices of the paper:
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

% Parameters definitions
eta = 120*pi; % intrinsic impedance
gama_0 = 0.57721566490153; % Euler constant
wavenumber = 2*pi/wavelength;  % wavenumber
d_l = 0.1*wavelength; % dipole length of each antenna
a = 10^(-5)*wavelength; % the radius of the wire
Z_L = 50; % load impedance
N = N_x*N_y;

%% Antenna impedance Z_A

eta1 = eta/(2*pi*sin(wavenumber*d_l/2)^2);
eta2 = eta/(4*pi*sin(wavenumber*d_l/2)^2);

% Resistance
R_ZA = eta1*(gama_0 + log(wavenumber*d_l) - cosint(wavenumber*d_l) + sin(wavenumber*d_l)/2*(sinint(2*wavenumber*d_l) - 2*sinint(wavenumber*d_l))...
    + cos(wavenumber*d_l)/2*(gama_0 + log(wavenumber*d_l/2) + cosint(2*wavenumber*d_l) - 2*cosint(wavenumber*d_l)));

% Reactance
X_ZA = eta2*(2*sinint(wavenumber*d_l) + cos(wavenumber*d_l)*(2*sinint(wavenumber*d_l) - sinint(2*wavenumber*d_l)) ...
    - sin(wavenumber*d_l)*(2*cosint(wavenumber*d_l) - cosint(2*wavenumber*d_l) - cosint(2*wavenumber*a^2/d_l)));

Z_A = R_ZA + 1i*X_ZA;

%% Mutual impedance matrix

Z_C = zeros(N,N);


for p = 1:N_y
    for q = 1:N_y

        for n1 = 1:N_x
            for n11 = 1:N_x

                d_pqy = abs((p-q)*AntennaSpacing_y);
                d_n1n11x = abs((n1-n11)*AntennaSpacing_x);


                u0 = wavenumber*d_pqy;
                u1 = wavenumber*(sqrt(d_n1n11x^2 + d_pqy^2) + d_pqy);
                u11 = wavenumber*(sqrt(d_n1n11x^2 + d_pqy^2) - d_pqy);
                u2 = wavenumber*(sqrt(d_n1n11x^2 + (d_pqy-d_l)^2) + (d_pqy-d_l));
                u22 = wavenumber*(sqrt(d_n1n11x^2 + (d_pqy-d_l)^2) - (d_pqy-d_l));
                u3 = wavenumber*(sqrt(d_n1n11x^2 + (d_pqy+d_l)^2) + (d_pqy-d_l));
                u33 = wavenumber*(sqrt(d_n1n11x^2 + (d_pqy+d_l)^2) - (d_pqy-d_l));
                u4 = 2*wavenumber*(d_pqy+d_l);
                u5 = 2*wavenumber*(d_pqy-d_l);
                u6 = (d_pqy^2-d_l^2)/d_pqy^2;

                eta3 = -eta/(8*pi)*cos(u0);
                eta4 = eta/(8*pi)*sin(u0);

                if p == q && n1 == n11

                    Z_C((p-1)*N_x + n1,(q-1)*N_x + n11) = Z_A;

                else
                    if n1 == n11

                        R = eta3*(-2*cosint(2*u0) + cosint(u5) + cosint(u4) - log(u6))...
                            + eta4*(2*sinint(2*u0) - sinint(u5) - sinint(u4));

                        X = eta3*(2*sinint(2*u0) - sinint(u5) - sinint(u4))...
                            + eta4*(2*cosint(2*u0) - cosint(u5) - cosint(u4) - log(u6));

                        zz = R + 1i*X;


                        Z_C((p-1)*N_x + n1,(q-1)*N_x + n11) = zz;

                    else



                        R = eta3*(-2*cosint(u1) - 2*cosint(u11) + cosint(u2)  + cosint(u22) + cosint(u3) + cosint(u33))...
                            + eta4*(2*sinint(u1) - 2*sinint(u11) - sinint(u2) + sinint(u22) - sinint(u3) + sinint(u33));

                        X = eta3*(2*sinint(u1) + 2*sinint(u11) - sinint(u2) - sinint(u22) - sinint(u3) - sinint(u33))...
                            + eta4*(2*cosint(u1) - 2*cosint(u11) - cosint(u2) + cosint(u22) - cosint(u3) + cosint(u33));

                        zz = R + 1i*X;



                        % All sub-arrays have the same configurations so
                        % that mutual coupling matrices for all sub-arrays
                        % are same

                        Z_C((p-1)*N_x + n1,(q-1)*N_x + n11) = zz;

                    end
                end

            end
        end
    end
end

% Mutual coupling matrix
Z = (Z_A + Z_L)\(Z_C + Z_L*eye(N));
