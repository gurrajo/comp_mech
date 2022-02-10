clear all; close all; clc
E_L = 181*10^9;
E_T = 10.3*10^9;
v_LT = 0.28;
v_TL = v_LT*E_T/E_L;
G_LT = 7.17*10^9;
G_TT = 3.5*10^9;
alpha_L = 1*10^-6;
alpha_T = 50*10^-6;

theta = [0,pi/2, 0, pi/2, pi/2, 0, pi/2, 0];
p0 = -2*10^3; % [Pa]
a = 0.6;
b = 0.5;
h = 1*10^-3;
t = h/8;
honey_t = t*0; % honeycomb thickness
h_tot = h + honey_t;
height(1,1) = h_tot/2;
height(2,1) = height(1,1) - t;
for i = 2:length(theta)+1
    if i == length(theta)/2 + 1
        t = honey_t;
    else
        t = h/8;
    end
    height(1,i) = height(2,i-1);
    height(2,i) = height(1,i) - t;
end

%% Calculate ply specific data
for i = 1:length(theta)
    [T_1{i},T_2{i},Q_bar{i}] = matrix_func(theta(i),E_L, E_T, v_LT,v_TL,G_LT);
end
for i = 1:length(theta)+1
    if i == 5
        % honeycomb bending stiffnes not included
        Q_bar_new{i} = zeros(3,3); 
    elseif i > 5
        Q_bar_new{i} = Q_bar{i-1};
    else
        Q_bar_new{i} = Q_bar{i};
    end
end

%% Caluculate laminate specific data
[A,B,D] = lamina_func(height,Q_bar_new);

%% Fourier analysis
m_max = 10;
n_max = 10;


w_xy = 0;
x = a/2;
y = b/2;
for m = 1:m_max
    for n = 1:n_max
        p_mn = 16*p0/(m*n*pi^2);
        w_mn = p_mn/(pi^4*(D(1,1)*(m/a)^4 + 2*(D(1,2) + 2*D(3,3))*(m/a)^2*(n/b)^2 + D(2,2)*(n/b)^4));
        w_xy = w_xy + w_mn*sin(m*pi*x/a)*sin(n*pi*y/b);
    end
end

disp("Displacemnt of middle point: " + w_xy*10^3 + " [mm]")
disp("Thickness of honeycomb section: " + honey_t*10^3 + " [mm]")
disp("Thickness of single ply: " + h/8*10^3 + " [mm]")

%% draw cross section
figure
hold on
for i = 1:length(height)
    if mod(i,2) == 0
        colour = 'r';
    elseif i == length(theta)/2 + 1
        colour = 'y';
    else
        colour = 'g';
    end
    fill([0,0,1,1], [height(:,i); flip(height(:,i))],colour)
end




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