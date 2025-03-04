function H = functionChannelGeneration(rr,rt,wavenumber,M,N,I)

%%=============================================================
%The file is used to generate the dyadic Green's function-based channel
%matrix of the paper:
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


for n = 1:N
    for m = 1:M

        G = functionDyadicGreen(rr(:,n),rt(:,m),wavenumber); 
        
        for i = 1:3
            for ii = 1:3
                H((i-1)*N+n,(ii-1)*M+m) = G(i,ii);
            end
        end
    end
end








