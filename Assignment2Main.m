%% Assignment 2 Sonya Stuhec-Leonard
clear
close all
%Part 1: solve by Finite difference a box of L X W
%number of nodes desired (mesh size)
L = 12; 
W = 8;

%Boundary condition values
V0 = 3;

%Create G matrix
G = sparse(L*W);

%B matrix defines boundary conditions
B = zeros(1, L*W);

for j = 1:L
    for i = 1:W
        %numbering scheme for G
        n = i+(j-1)*W;
        
        if j==1 
            %at x=0 Bc = V0, and diagonal of BC = 1
            %G(n, :) = V0;
            G(n, n) = V0;
            B(n) = V0;
            
        elseif j==L
            %at x=L BC=0, and diagonal of BC = 1
           % G(n, :) = 0;
            G(n, n) = 1;
            B(n) = 0;
            
        elseif i==0 || i==W
            %G(n, :) = 0;
            G(n, n) = 1;
            B(n) = 0;
            
        else
            %solve FD equations or put 1 an d4 in the row coresponding to
            %the i, j position
            nxMinus1 = i+(j-2)*W;
            nxPlus1 = i+ (j)*W;
            nyMinus1 = i-1+(j-1)*W;
            nyPlus1 = i+1+(j-1)*W;
            
            G(n, n)= -4;
            G(n, nxMinus1) = 1;
            G(n, nxPlus1) = 1;
            G(n, nyMinus1) = 1;
            G(n, nyPlus1) = 1;
            
        end
    end
end

%usful for quick verification of G
figure(1)
spy(G)
title('Spying on G matrix')

V = G\B'; %inverse(G)*B


% Define Eigen vector of nx by ny plane

for j=1:L
    for i=1:W
        n=i+(j-1)*W;
        %chaning first col of E (eigen vetor) to the nx X ny matrix for
        %plotting. Chaningn which col of E you convert will change withc
        %eigen vector you look at
        V_vec(j, i) = V(n); 
    end
end

[Ex, Ey] = gradient(V_vec);

figure
surf(V_vec)

% figure(3) 
% surf(Ex, Ey, '.')
% title('Plot of e-feild in x-y plane')

% 1b) anlytical serise solution for comparison
% x= L;
% y= W;
% a= 1;
% b= 1;
% 
% Sum_1 = cosh(n*pi.*x/a);
% Sum_2 = cosh(n*pi.*b/a);
% Sum_3 = cosh(n*pi.*y/a);

% V_ana = @(n) 4*V0/pi *sum((1/m*cosh(m*pi.*x/a))

