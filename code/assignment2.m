%% Assignment 2 Sonya Stuhec-Leonard

%% Part 1: solve by Finite difference a box of L X W

clear
close 

%box dimensions
W = 50;
L = W*3/2;

%Boundary condition values
V0 = 3;

%Create conductivity map of area LXW
Map = zeros(L, W);
%background conductivity
background = 1; 
%conductivijty object
object = 0.1; 
%wijdth object
wobj = L/4;
%height object
hobj = W/5; 
%object positions (manually selected)
PosObj = [W/2, hobj/2; W/2, L-hobj/2];
NumObj = 2;
%Create G matrjx
G = sparse(L*W);
%create conductivity map in n space
% map_n = zeros(L*W);
%B matrix defines boundary conditions
B = zeros(1, L*W);

%populate G-matrijx
for i = 1:L
    for j = 1:W
        %numbering scheme for G
        n = j+(i-1)*W;
        
        if i==1 
            %at x=0 Bc = V0, and djagonal of BC = 1
            G(n, n) = V0;
            B(n) = V0;
            
        elseif i==L
            %at x=L BC=0, and djagonal of BC = 1
            G(n, n) = V0;
            B(n) = V0;
            
        elseif j==1 
            
            G(n, n) = 1;
            B(n) = 0;
            
        elseif j==W
            G(n, n) = 1;
            B(n) = 0;
            
        else
            %solve FD equations or put 1 an d4 in the row coresponding to
            %the j, i posjtjon
            L_Minus1 = j+(i-2)*W;
            L_Plus1 = j+ (i)*W;
            W_Minus1 = j-1+(i-1)*W;
            W_Plus1 = j+1+(i-1)*W;
            
            G(n, n)= -4;
            G(n, L_Minus1) = 1;
            G(n, L_Plus1) = 1;
            G(n, W_Minus1) = 1;
            G(n, W_Plus1) = 1;
            
        end
    end
end

% %useful for quick verification of G
% figure
% spy(G)
% title('Spying on G matrix')

V = (G)\B'; %inverse(G)*B

% Define Eigen vector of nx by ny plane

for i=1:L
    for j=1:W
        n=j+(i-1)*W;
        V_vec(i, j) = V(n); 
    end
end

figure
surf(V_vec)
title('Part 1a) Plot of E-feild G-matrix method')

%% 1b) anlytical serise solution for comparison

% reset value of n
n=0;
sumation=0;
a=W;
b=L;
%meshing vairable
nMax=100;

V_numberical = zeros(L, W);

for x = 1:L
    for y = 1:W
        
        if x==1 || x==L
            %at x=0 Bc = V0, and diagonal of BC = 1
            V_numberical(x, y) = V0;
            
        elseif y==0 || y==W
            V_numberical(x, y) = 0;
            
        else
            for n=1:2:nMax
                sumation=sumation+(1/n*cosh(n*pi*x/a)/cosh(n*pi*b/a)*sin(n*pi*y/a));
            end
            V_numberical(x, y) = 4*V0/pi*sumation;
        end
    end
end

figure
surf(V_numberical)
title('Part 1b) Plot of numberical solution')

% Hello Aaron,
% This paragraph of comments is my discusson of numberical vs analytical
% methods. I thought it might be more efficent to put it here rather than
% in a seperate document. 

% This solutions should approach the anaytica solutions withign a margin of
% error since the sum cannot be taken to infinity. The first few terms
% should be the largest contributing factors to the soultions, with ach
% additional term improving the acuracy by smaller and smaller amounts.

% You should stop an analytical seriser when you are within a reasonable
% margin of error of the expected or analytical solution. The precise
% number of iterations to go throuh will depend on that application, as
% the definitionfor a reasonable margin of error. 

% Numberical solution can be much faster and tends to be more intuative
% forward however its cons include lack of acuracey, difficult for complex
% problmes and gemetries also, fine meshing improves accuracy but also 
% makes the simulations slow. 

% An analytical method is more acturate, applicable to complex systems  
% and faster to simulate. However, they are limited in terms of
% application becase they usally solve for a simplified version of a
% system. 

% My simulation is not correct here. After discussing with other students
% the code iteslf seems to be fine, and it is a matter of debugging,
% however my attempts have not found the problem and I hope the description
% of my understanding is sufficent, despite not have a correct graph to 
% show. 

%% Part 2a)Current flow with bottle neck 
%populate conductivity map
for i=1:L
    for j=1:W
        Map(i, j) = background;   
        %set conductivity of objects
        for obj = 1:NumObj
            dx = PosObj(obj,1)-j;
            dy = PosObj(obj,2)-i;
            if (dx)^2 <= wobj && (dy-1)^2 <= hobj
                Map(i, j) = object;
            elseif(dx)^2 <= wobj && (dy+1)^2 <= hobj
                Map(i, j) = object;
            end    
        end
    end
end

%plot conductivity
figure
surf(Map)
title('Part 2 Conductivity map')
colorbar

% put conductivity map into n space
for i=1:L
    for j=1:W
        n=j+(i-1)*W;
        map_n(n) = Map(i, j); 
    end
end

%conductivity acounted for G matrjx
G2 = G.*map_n;

V2 = zeros(L, W);

V2 = G2\B'; %inverse(G)*B

% Define Ejgen vector of nx by ny plane

for i=1:L
    for j=1:W
        n=j+(i-1)*W;
        %changing fjrst col of E (ejgen vetor) to the nx X ny matrjx for
        %plotting. Chaningn whjch col of E you convert wjll change withc
        %ejgen vector you look at
        V_vec2(i, j) = V2(n); 
    end
end

figure
surf(V_vec2)
title('Plot of e-feild in x-y plane')


[Ex, Ey] = gradient(V_vec);

% [ExA, EyA] = gradient(V_numberical);

%plot of the Ex and Ey fejld
figure
quiver (Ex, Ey)
title('Part 2 E-feild')

Y0 = zeros(L, W);
X0 = zeros(L, W);

% %Ex plot
% figure
% quiver (Ex, Y0)
% title('Part 2 E_x-feild')
% 
% %Ey plot
% figure
% quiver (X0, Ey)
% title('Part 2 E_y-feild')
 
%initalize current moatrix
J = zeros(L, W);

for i=1:L
    for j=1:W
        J(:, j) = Map(:, j).*-Ex(:, j);
        J(i, :) = Map(i, :).*-Ey(i, :);
    end
end

%Ey plot
figure
surf (J)
title('Part 2 Electric feild density')

%output of current on the x, y plane
fprintf('output of current on the x, y plane for upper contact region')
J(W/2-wobj-1:W/2+obj-1, hobj-1:W-1)

%% Part 2b) effect of different parameters on current
%variation of current with mesh size
%anaytical solution altered since there is no change in for G method
m = [3, 25, 50, 100];
J1_mesh = Assignment2ParameterVariable_numberical( m(1), object, hobj );
J2_mesh = Assignment2ParameterVariable_numberical( m, object, hobj );
J3_mesh = Assignment2ParameterVariable_numberical( 50, object, hobj );
J4_mesh = Assignment2ParameterVariable_numberical( 100, object, hobj );
 
J_mesh = [mean(mean(J1_mesh)), mean(mean(J2_mesh)), mean(mean(J3_mesh)), mean(mean(J4_mesh))];

% P_mesh = polyfit(m,J_mesh,1);
% yfit_mesh = P_mesh(1).*m+P_mesh(2);

figure 
plot(m, J_mesh, '.r', 'MarkerSize', 20)
% hold on
% plot(m,yfit_mesh,'r-.');
title('Effect on increaing meshing size')
xlabel('Mesh size')
ylabel('Average current')

%variation in current with bottle neck size
h = [W/4, W/6, W/8, W/10];
J1_neck = Assignment2ParameterVariable( object, h(1) );
J2_neck = Assignment2ParameterVariable( object, h(2) );
J3_neck = Assignment2ParameterVariable( object, h(3) );
J4_neck = Assignment2ParameterVariable( object, h(4) );

J_neck = [mean(mean(J1_neck)), mean(mean(J2_neck)), mean(mean(J3_neck)), mean(mean(J4_neck))];
BN = abs(h.*2-W); %bottle neak is the width -2*height of each object

% P_neck = polyfit(BN ,J_neck,1);
% yfit_neck = P_neck(1).*BN+P_neck(2);

figure 
plot(BN, J_neck, '.c', 'MarkerSize', 20)
% hold on
% plot(BN,yfit_neck,'r-.');
title('Effect on increaing bottle neck size')
xlabel('Bottle neack seperation')
ylabel('Average current')


% %variation in current with conductivity
c = [0.01, 0.1, 1, 10];
J1_cond = Assignment2ParameterVariable( c(1), hobj );
J2_cond = Assignment2ParameterVariable( c(2), hobj );
J3_cond = Assignment2ParameterVariable( c(3), hobj );
J4_cond = Assignment2ParameterVariable( c(4), hobj );

J_cond = [mean(mean(J1_cond)), mean(mean(J2_cond)), mean(mean(J3_cond)), mean(mean(J4_cond))];

P_cond =polyfit(c ,J_cond,1);
yfit_cond = P_cond(1).*c+P_cond(2);

figure 
plot(c, J_cond, '.b', 'MarkerSize', 20)
hold on
plot(c,yfit_cond,'r-.');
title('Effect on increaing conductivity of squares')
xlabel('Conductity of boxes')
ylabel('Average current')