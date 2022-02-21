clear all; close all; clc

%% Material data
E_f = 350*10^9;
E_m = 3.5*10^9;
alpha_m = 50*10^-6;
alpha_f = -1*10^-6;
V_f = 0.55;
V_m = 1-V_f;
v_f = 0.2;
v_m = 0.35;
G_f = E_f/(2*(1+v_f));
G_m = E_m/(2*(1+v_m));

%% Convert to LT system
%halpin tsai
xi = 2;
eta = ((E_f/E_m) - 1)/(E_f/E_m + xi);
E_T = (1+xi*eta*V_f)*E_m/(1-eta*V_f);

%ROM
E_L = E_f*V_f + E_m*V_m;

eta_G = ((G_f/G_m) - 1)/((G_f/G_m) + 1);
G_LT = (1+eta_G*V_f)*G_m/(1-eta_G*V_f);

v_LT = v_f*V_f+V_m*v_m;
v_TL = v_LT*E_T/E_L;

alpha_L = 1/E_L*(alpha_f*E_f*V_f+alpha_m*E_m*V_m);
alpha_T = (1+v_f)*alpha_f*V_f + (1+v_m)*alpha_m*V_m-alpha_L*v_LT;

%% Other question specific data
theta = [0, 45*pi/180, 90*pi/180, -45*pi/180 , -45*pi/180, 90*pi/180 , 45*pi/180, 0]; % [rad]
del_T = -90;
eps_0 = 0.001;

M = [60 , 60 , 0]'; % [Nm/m]

% matrix such that each coloumn represents the start and end point for
% ply(k)
height(1,1:8) = linspace(0.6*10^-3,-0.6*3/4*10^-3,8);
height(2,1:8) = linspace(0.6*3/4*10^-3,-0.6*10^-3,8);

%% Calculate ply specific data
for i = 1:length(theta)
    [T_1{i},T_2{i},Q_bar{i}] = matrix_func(theta(i),E_L, E_T, v_LT,v_TL,G_LT);
end

%% Caluculate laminate specific data
[A,B,D] = lamina_func(height,Q_bar);

%% Calculate force from temperature change in each ply
cooling_force = 0;
for i = 1:length(theta)
    alpha_xy{i} = inv(T_2{i})*[alpha_L , alpha_T , 0]';
    cooling_force = cooling_force + Q_bar{i}*alpha_xy{i}*del_T*(height(1,i) - height(2,i));
end
% temperature related strain
eps_0_T = inv(A)*cooling_force;


%% Part 2, applied load
% solve linear system to obtain y and xy strains from total load
syms N eps_y_N eps_xy_N
sol = solve(A*[eps_0, eps_y_N, eps_xy_N]' == ([N, 0 , 0]' + cooling_force));
eps_y_N = double(sol.eps_y_N);
eps_xy_N = double(sol.eps_xy_N);

k = inv(D)*M;

z = linspace(height(1,1),height(2,end),10000);
height_limit = height(2,1);
j = 1;
for i = 1:length(z)
    if z(i) < height_limit
        % change ply index when z < than the end height of current ply
        j = j + 1 ;
        height_limit = height(2,j);
    end
    % stress from temperature load only
    stress_T(i,:) = Q_bar{j}*(eps_0_T-alpha_xy{j}*del_T); 
   
    % stress from both temperature and applied load
    stress_xy{i} = Q_bar{j}*[eps_0 + k(1)*z(i), eps_y_N + k(2)*z(i) , eps_xy_N + k(3)*z(i)]' - Q_bar{j}*alpha_xy{j}*del_T;
    stress_LT(i,:) = T_1{j}*stress_xy{i};
end

%% display max and min stress values for temperature and applied load
disp("maximum longitudinal stress: " + max(stress_LT(:,1))*10^-6 + " [MPa]")
disp("maximum transverse stress: " + max(stress_LT(:,2)*10^-6) + " [MPa]")
disp("maximum shear stress : " + max(stress_LT(:,3))*10^-6 + " [MPa]")
disp("minimum longitudinal stress: " + min(stress_LT(:,1))*10^-6 + " [MPa]")
disp("minimum transverse stress: " + min(stress_LT(:,2)*10^-6) + " [MPa]")
disp("minimum shear stress : " + min(stress_LT(:,3))*10^-6 + " [MPa]")


%% Plots
% note the postive and negative direction of x and y axis. 
figure
plot(stress_T, z)
hold on
title("Only temperature load")
legend("\sigma_x","\sigma_y","\tau_x_y")
set(gca,'xdir','reverse','ydir','reverse')
ylabel("z [m]")
xlabel("stress [Pa]")

figure
plot(stress_LT, z)
hold on
title("Temperature and applied load")
legend("\sigma_L","\sigma_T","\tau_L_T")
set(gca,'xdir','reverse','ydir','reverse')
ylabel("z [m]")
xlabel("stress [Pa]")


function[T_1,T_2,Q_bar] = matrix_func(theta, E_L, E_T, v_LT,v_TL,G_LT)
%% calculates ply matrixes from material data in LT system
% input: angle [rad], E_L, E_T, G_LT [Pa]
% output: T_1,T_2,Q_bar
T_1 = [cos(theta)^2, sin(theta)^2, 2*sin(theta)*cos(theta);
       sin(theta)^2, cos(theta)^2, -2*sin(theta)*cos(theta);
      -sin(theta)*cos(theta), sin(theta)*cos(theta), (cos(theta)^2 - sin(theta)^2)];
T_2 = [cos(theta)^2, sin(theta)^2, sin(theta)*cos(theta) ;
       sin(theta)^2, cos(theta)^2, -sin(theta)*cos(theta);
       -2*sin(theta)*cos(theta), 2*sin(theta)*cos(theta), (cos(theta)^2-sin(theta)^2)];

Q = [E_L/(1-v_LT*v_TL), v_TL*E_L/(1-v_LT*v_TL), 0; 
     v_LT*E_T/(1-v_LT*v_TL), E_T/(1-v_LT*v_TL), 0; 
     0,  0, G_LT];
Q_bar = inv(T_1)*Q*T_2;
end

function[A,B,D] = lamina_func(height,Q_bar)
%% calculate laminate specific matrixes A,B,D
A = zeros(3,3);
B = zeros(3,3);
D = zeros(3,3);
for i=1:length(height)
    A = A + Q_bar{i}*(height(1,i) - height(2,i));
    B = B + (1/2)*Q_bar{i}*(height(1,i)^2 - height(2,i)^2);
    D = D + (1/3)*Q_bar{i}*(height(1,i)^3 - height(2,i)^3); 
end
end
