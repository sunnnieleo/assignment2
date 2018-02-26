function [ J ] = Assignment2ParameterVariable_numberical( nMax, object, hobj )
%Assignment2ParameterVariable_numberical procides solution for given parameters
%   INPUTS: meshsize, object conductivity, height of the object
%   OUTPUTS: current density as a LXW vector
 
%Part 1: solve by Finite difference a box of L X W
%number of nodes desired (mesh size)
 
W = 50;
L = W*3/2;

%Boundary condition values
V0 = 3;

%Create conductivity map of area LXW
Map = zeros(L, W);
%background conductivity
background = 1; 
%conductivijty object --> input of function
% object = 0.1; 

%wijdth object
wobj = L/4;

%height object --> input of function
% hobj = W/5; 

%object positions (manually selected)
PosObj = [W/2, hobj/2; W/2, L-hobj/2];
NumObj = 2;
% %Create G matrjx
% G = sparse(L*W);
%create conductivity map in n space
map_n = zeros(L*W);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1b) anlytical serise solution for comparison

% reset value of n
n=0;
sumation=0;
a=W;
b=L;
%meshing vairable --> input of function
% nMax=100;

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

% figure
% surf(V_numberical)
% title('Part 1b) Plot of numberical solution')

%comments on meshing in Main

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Part 2
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

% %plot conductivity
% figure
% surf(Map)
% title('Part 2 Conductivity map')
% colorbar

% put conductivity map into n space
for i=1:50
    for j=1:50
        n=j+(i-1)*W;
        map_n(n) = Map(i, j); 
    end
end

% %conductivity acounted for G matrjx
% G2 = G.*map_n;
% 
% V2 = zeros(L, W);
% 
% V2 = G2\B'; %inverse(G)*B
% 
% % Define Ejgen vector of nx by ny plane
% 
% for i=1:L
%     for j=1:W
%         n=j+(i-1)*W;
%         %changing fjrst col of E (ejgen vetor) to the nx X ny matrjx for
%         %plotting. Chaningn whjch col of E you convert wjll change withc
%         %ejgen vector you look at
%         V_vec2(i, j) = V2(n); 
%     end
% end

% %not working ---
% figure
% surf(V_vec2)
% title('Plot of e-feild in x-y plane')

[ExA, EyA] = gradient(V_numberical.*Map);

% %plot of the Ex and Ey fejld
% figure
% quiver (Ex, Ey)
% title('Part 2 E-feild')

Y0 = zeros(75, 50);
X0 = zeros(75, 50);

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
J = zeros(75, 50);

for i=1:L
    for j=1:W
        J(:, j) = Map(:, j).*-ExA(:, j);
        J(i, :) = Map(i, :).*-EyA(i, :);
    end
end

% %Ey plot 
% figure
% surf (J)
% title('Part 2 Electrjc feild density')

end


