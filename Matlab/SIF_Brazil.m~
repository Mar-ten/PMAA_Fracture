% Martin Raming
% Calcultae stress intensity factors for brizil disk for Mode I and ModeII
% 4/25/2018
% P is force at crack propagation
% theta is inclination angle in degrees
function [K_I,K_II,Mix_I,Mix_II] = SIF_Brazil(P,theta)

a = 4; % [mm] Half of the intial crack length
b = 3; % [mm] thickness of specimen
r = 25/2; % [mm] radius of disk
ratio = a/r;

% Ai theta and Bi thetas constants from table found in Atkinson
T_3 = [1.135551 0.533477 0.39164 0.393835 0.325033]; % a/r = .3
T_4 = [1.24314 0.559734 0.404603 0.408597 0.334831];%a/r = .4
S_3= [1.089702 0.522272 0.386086 0.387518 0.320834];% a/r = .3
S_4= [1.160796 0.539824 0.394822 0.397403 0.327411];%a/r = .4

% Interploate for actual a/r
T= zeros(5);
S=zeros(5);
for i = 1:length(T_3)
    T(i) = (T_4(i)-T_3(i))/(.4-.3)*(ratio-.3) + T_3(i);
    S(i) = (S_4(i)-S_3(i))/(.4-.3)*(ratio-.3) + S_3(i);
end

s = sind(theta); % not sure if this actually supposed to be theata????
c = cosd(theta);
%Calculate A's
A(1) = 1-4.*s.^2;
A(2) = 8.*s.^2.*(1-4*c^2);
A(3) = -4.*s.^2.*(3-36.*c.^2+48.*c.^4);
A(4) = -16*s.^2.*(-1 +24*c.^2 - 80*c.^4+64*c.^6);
A(5) = -20*s.^2.*(1-40*c.^2+240*c.^4-448*c.^6+256*c.^8);
%Calculate B's
B(1) = 1;
B(2) = -5 +8*c.^2;
B(3) = -3+8*(1-2*c.^2).*(2-3*c.^2);
B(4) = 3+16*(1-2*c.^2)-12*(1-2*c.^2).^2-32*(1-2*c.^2).^3;
B(5) = 5-16*(1-2*c.^2)-60*(1-2*c.^2).^2+32*(1-2*c.^2).^3+80*(1-2*c.^2).^4;

% Calculate N's
for k = 1:length(T)
    
    N_i(k) = T(k).*ratio.^(2*k - 2)*A(k);
    N_ii(k) = S(k)*ratio.^(2*k-2).*B(k);
end
N_I = sum(N_i); % sum all N's for SIF mode I
N_II = 2*sind(2.*theta)*sum(N_ii);% sum all N's for SIF mode II

K_I = P.*sqrt(a).*N_I/(sqrt(pi)*r*b); % Mode I SIF
K_II = P.*sqrt(a).*N_II/(sqrt(pi)*r*b);% Mode II SIF
Mix_I = K_I/(K_II+K_I); % Percent Mode I
Mix_II = K_II/(K_II+K_I); % Percent Mode II

end





