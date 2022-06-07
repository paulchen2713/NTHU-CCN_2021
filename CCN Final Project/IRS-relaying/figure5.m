%
% Figure 5 - The energy efficiency as a function of the rate R¯
% ref. E. Björnson, Ö. Özdogan and E. G. Larsson, "Intelligent Reflecting Surface 
%      Versus Decode-and-Forward: How Large Surfaces are Needed to Beat Relaying?," 
%      in IEEE Wireless Communications Letters, vol.9, no.2, pp.244-248, Feb.2020
%
% close all;
clear;
clc;
%
% set simulation parameters
%
% carrier frequency f_c = 3GHz
%
fc = 3;
%
% bandwidth B = 10 MHz
%
B = 10e6;
%
% noise figure in dB
%
noiseFiguredB = 10;
%
% compute the noise power in dBm
%
sigma2dBm = -174 + 10*log10(B) + noiseFiguredB;
sigma2 = db2pow(sigma2dBm);
%
% define the channel gain functions based on the 3GPP Urban Micro (UMi)
%
% pathloss_3GPP_LOS  = @(x) db2pow(-28.0 - 20*log10(fc) - 22.0*log10(x));
% pathloss_3GPP_NLOS = @(x) db2pow(-22.7 - 26*log10(fc) - 36.7*log10(x));
%
% define the antenna gains at the source, relay/IRS, and destination. 
% (the numbers are in linear scale)
%
Gs = db2pow(5); % antenna gains at the source
Gr = db2pow(5); % antenna gains at the relay/IRS
Gd = db2pow(0); % antenna gains at the destination
%
% set the amplitude reflection coefficient
%
alpha = 1;
%
% set the range of rate values
%
Rbar = [0.01 0.1:0.1:10];
%
% set parameters related to circuit power consumption
%
Ps = 100; % Power dissipation in the transceiver hardware of the source
Pd = 100; % Power dissipation in the transceiver hardware of the destination
Pe = 5;   % Power dissipation per element in the IRS (mW)
Pr = 100; % Power dissipation in the transceiver hardware of the relay
nu = 0.5; % Efficiency of the power amplifier at the source
%
% define distances in simulation setup
%
d_SR = 80; % distance between the source and IRS/relay
dv = 10;   % minimum distance between destination and the IRS/relay
%
% define the range of d1 values in the simulation setup
%
d1 = 70;
%
% copmute distance between the source and destination
%
d_SD = sqrt(d1^2+dv^2);
%
% compute distance between the IRS/relay and destination
%
d_RD = sqrt((d1-d_SR)^2+dv^2);
%
% compute the channel gains using the 3GPP models and antenna gains
%
betaSR = pathloss_3GPP_LOS(d_SR, fc) * Gs * Gr;  % β_sr
betaRD = pathloss_3GPP_LOS(d_RD, fc) * Gr * Gd;  % β_rd
betaSD = pathloss_3GPP_NLOS(d_SD, fc) * Gs * Gd; % β_sd
betaIRS = betaSR * betaRD;  % β_IRS = β_sr * β_rd
%
% prepare to save simulation results
%
EE_SISO = zeros(length(Rbar), 1);
EE_IRS  = zeros(length(Rbar), 1);
EE_DF   = zeros(length(Rbar), 1);
Nopt    = zeros(length(Rbar), 1);
%
% go through all rate values
%
for i = 1:length(Rbar)
    %
    % compute required SINR values
    %
    SINR = 2^(Rbar(i)) - 1;        % SISO and IRS
    SINR_DF = 2^(2 * Rbar(i)) - 1; % DF relaying
    %
    % compute the transmit power in the SISO case, using Eq.(17)
    %
    P_SISO = SINR * sigma2 / betaSD;
    %
    % compute the energy efficiency in the SISO case, using B * R¯ / P_total
    % 
    % EE = B * R¯ / P_total (the factor 1000 is used to convert mW to W)
    EE_SISO(i) = 1000 * B * Rbar(i) / (P_SISO/nu + Ps + Pd);
    %
    % compute the transmit power in the DF relaying case, using Eq.(19)
    %
    P_DF = SINR_DF * sigma2 * (betaSR + betaRD - betaSD) / (2 * betaIRS);
    %
    % compute the energy efficiency in the DF relaying case
    %
    % EE = B * R¯ / P_total (the factor 1000 is used to convert mW to W)
    EE_DF(i) = 1000 * B * Rbar(i) / (P_DF/nu + Ps/2 + Pd + Pr);
    %
    % compute the power-minimizing number of reflecting elements, using Eq.(22)
    %
    Nopt(i) = (2 * SINR * sigma2 / (alpha^2 * betaIRS * Pe))^(1/3) ...
                - sqrt(betaSD / betaIRS) / alpha;
    if Nopt(i) < 0
        Nopt(i) = 0;
    end
    %
    % compute the transmit power in the IRS case, using Eq.(18)
    %
    P_IRS = SINR * sigma2 ./ (sqrt(betaSD) + Nopt(i) * alpha * sqrt(betaIRS)).^2;
    %
    % compute the energy efficiency in the IRS case, using B * R¯ / P_total
    %
    % EE = B * R¯ / P_total (the factor 1000 is used to convert mW to W)
    EE_IRS(i) = 1000 * B * Rbar(i) / (P_IRS/nu + Ps + Pd + Nopt(i)*Pe);
end
%
% plot simulation results
%
figure;
hold on; box on; grid on;
%
plot(Rbar, EE_DF/1e6, 'b-.', 'LineWidth', 2);
plot(Rbar, EE_IRS/1e6, 'r-', 'LineWidth', 2);
plot(Rbar, EE_SISO/1e6, 'k--', 'LineWidth', 2);
%
title('Figure 5');
xlabel('Achievable rate [bit/s/Hz]');
ylabel('Energy efficiency [Mbit/Joule]');
%
legend('DF relay', 'IRS', 'SISO', 'Location', 'NorthWest');
set(gca,'fontsize', 12);
%
% define the channel gain functions based on the 3GPP Urban Micro (UMi)
%
function out = pathloss_3GPP_LOS(x, fc)
    % x is measured in m, antenna gains are included separately in the code
    out = db2pow(-28.0 - 20*log10(fc) - 22.0*log10(x));
end
function out = pathloss_3GPP_NLOS(x, fc)
    % x is measured in m, antenna gains are included separately in the code
    out = db2pow(-22.7 - 26*log10(fc) - 36.7*log10(x));
end
% 
