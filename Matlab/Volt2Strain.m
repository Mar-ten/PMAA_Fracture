function [strain]=Volt2Strain(V_ch)
% Martin Raming 4/2018
%this function calculates strain in aluminum for a full-bridge configuration

    V_ex=12; %excitation Voltage
    GF =2.125; % gauge Factor
    nu =0.345;  % poissons ratio
    u = 20.83333; %[mV/V}dewtron Setting
    Gain = 5/((V_ex/1000)*u);
    V_f = V_ch./(V_ex*Gain);

    strain = -2*V_f./(GF*(nu+1)-V_f.*(nu-1));
end 