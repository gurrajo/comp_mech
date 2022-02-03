clear all; close all; clc

theta = [0, 45*pi/180, 90*pi/180, -45*pi/180 , -45*pi/180, 90*pi/180 , 45*pi/180, 0];
E_f = 350*10^9;
alpha_m = 50*10^-6;
alpha_f = -1*10^-6;
E_m = 3.5*10^9;
V_f = 0.55;
V_m = 1-V_f;
v_f = 0.2;
v_m = 0.35;
rho_f = 2.0; %[g/cm^3]
G_f = E_f/(2*(1+v_f));
G_m = E_m/(2*(1+v_m));
eta_G = ((G_f/G_m) - 1)/((G_f/G_m) + 1);
G_LT = (1+eta_G*V_f)*G_m/(1-eta_G*V_f);



xi = 2;
eta = ((E_f/E_m) - 1)/(E_f/E_m + xi);
E_T = (1+xi*eta*V_f)*E_m/(1-eta*V_f);
E_L = E_f*V_f + E_m*V_m;

v_LT = v_f*V_f+V_m*v_m;
v_TL = v_LT*E_T/E_L;

alpha_L = 1/E_L*(alpha_f*E_f*V_f+alpha_m*E_m*V_m);
alpha_T = (1+v_f)*alpha_f*V_f + (1+v_m)*alpha_m*V_m-alpha_L*v_LT;

del_T = -90;
eps_0 = 0.001;

M = [60 , 60 , 0]'; % [Nm/m]

height(1,1:8) = linspace(4*1.2*10^-3,-3*1.2*10^-3,8);
height(2,1:8) = linspace(3*1.2*10^-3,-4*1.2*10^-3,8);

for i = 1:length(theta)
    [T_1{i},T_2{i},Q_bar{i}] = matrix_func(theta(i),E_L,v_LT,v_TL,E_T,G_LT);
end

[A,B,D] = lamina_func(height,Q_bar);

cooling_force = 0;
for i = 1:length(theta)
    alpha_xy{i} = inv(T_2{i})*[alpha_L , alpha_T , 0]';
    cooling_force = cooling_force + Q_bar{i}*alpha_xy{i}*del_T*(height(1,i) - height(2,i));
end
eps_0_T = inv(A)*cooling_force;

for i = 1:length(theta)
    strain_mech_T{i} = eps_0_T-alpha_xy{i}*del_T;
    stress_T{i} = Q_bar{i}*strain_mech_T{i};
end


k = inv(D)*M;

for i = 1:length(theta)
end

z = linspace(height(1,1),height(2,end),100);
height_limit = height(2,1);
j = 1;
for i = 1:length(z)
    if z(i) < height_limit
        j = j + 1 ;
        height_limit = height(2,j);
    end
    stress_xy{i} = Q_bar{j}*[eps_0 + k(1)*z(i), eps_0_T(2) + k(2)*z(i) , eps_0_T(3) + k(3)*z(i)]' - Q_bar{j}*alpha_xy{j}*del_T;
    stress_LT(i,:) = T_1{j}*stress_xy{i};
end

disp("maximum longitudinal stress: " + max(stress_LT(:,1))*10^-6 + " [MPa]")
disp("maximum transverse stress: " + max(stress_LT(:,2)*10^-6) + " [MPa]")
disp("maximum shear stress : " + max(stress_LT(:,3))*10^-6 + " [MPa]")
disp("minimum longitudinal stress: " + min(stress_LT(:,1))*10^-6 + " [MPa]")
disp("minimum transverse stress: " + min(stress_LT(:,2)*10^-6) + " [MPa]")
disp("minimum shear stress : " + min(stress_LT(:,3))*10^-6 + " [MPa]")

function[T_1,T_2,Q_bar] = matrix_func(theta, E_L, v_LT,v_TL, E_T,G_LT)
T_1 = [cos(theta).^2, sin(theta)^2, 2*sin(theta)*cos(theta); sin(theta)^2, cos(theta)^2, -2*sin(theta)*cos(theta); -sin(theta)*cos(theta), sin(theta)*cos(theta), (cos(theta)^2 - sin(theta)^2)];
T_2 = [cos(theta)^2 sin(theta)^2 sin(theta)*cos(theta)
       sin(theta)^2 cos(theta)^2 -sin(theta)*cos(theta)
       -2*sin(theta)*cos(theta) 2*sin(theta)*cos(theta) cos(theta)^2-sin(theta)^2];

Q = [E_L/(1-v_LT*v_TL) v_TL*E_L/(1-v_LT*v_TL) 0 
     v_TL*E_L/(1-v_LT*v_TL) E_T/(1-v_LT*v_TL) 0 
     0  0   G_LT];
Q_bar = T_2^-1*Q*T_1;
end
function[A,B,D] = lamina_func(height,Q_bar)
A = zeros(3,3);
B = zeros(3,3);
D = zeros(3,3);
for i=1:length(height)
    A = A + Q_bar{i}*(height(1,i) - height(2,i));
    B = B + 1/2*Q_bar{i}*(height(1,i)^2 - height(2,i)^2);
    D = D + 1/3*Q_bar{i}*(height(1,i)^3 - height(2,i)^3); 
end
end
