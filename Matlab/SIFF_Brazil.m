% Martin Raming
% Calcultae stress intensity factors for brizil disk for Mode I and ModeII
% 4/25/2018

a = 4; % [mm] Half of the intial crack length
b = 3; % [mm] thickness of specimen
r = 12.5; % [mm] radius of disk
% Ai theta and Bi thetas constants from Atkinson
T_3 = [1.135551 0.533477 0.39164 .0325033]; % a/r = .3
T_4 = [1.24314 0.559734 0.404603 0.408597 0.334831];%a/r = .4
S_3= [1.089702 0.522272 0.386086 0.387518 0.320834];% a/r = .3
S_4= [1.160796 0.539824 0.394822 0.397403 0.327411];;%a/r = .4
% Interploate for actual a/r

ratio = a/r;
T = (T_4(1)-T_3(1))/(.4-.3)*(ratio-.3) + T_3(1)