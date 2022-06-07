%
% Figure 2 - Typical channel gains as a function of the distance, when 
%            including the antenna gains Gt = Gr = 5 dBi.
% ref. E. Björnson, Ö. Özdogan and E. G. Larsson, "Intelligent Reflecting Surface 
%      Versus Decode-and-Forward: How Large Surfaces are Needed to Beat Relaying?," 
%      in IEEE Wireless Communications Letters, vol.9, no.2, pp.244-248, Feb.2020
%
% close all;
clear;
clc;
%
% set parameter values
%
% carrier frequency f_c = 3GHz
%
fc = 3; 
%
% distances in meter
%
d = 10:0.1:100;
%
% define the antenna gains at the transmitter and receiver
%
Gt = 5; % Tx antenna Gain in dBi
Gr = 5; % Rx antenna Gain in dBi
%
% compute channel gains based on the 3GPP Urban Micro
% Note that the antenna gains are included.
%
beta_3GPP_LOS  = db2pow(Gt + Gr - 28.0 - 20 * log10(fc) - 22.0 * log10(d));
beta_3GPP_NLOS = db2pow(Gt + Gr - 22.7 - 26 * log10(fc) - 36.7 * log10(d));
%
% plot simulation results
%
figure;
hold on; box on; grid on;
%
plot(d, 10 * log10(beta_3GPP_LOS), 'k-.', 'LineWidth', 2);
plot(d, 10 * log10(beta_3GPP_NLOS), 'r--', 'LineWidth', 2);
%
title('Figure 2');
xlabel('Distance d [m]');
ylabel('Channel gain \beta(d) [dB]');
%
legend({'UMi-LOS', 'UMi-NLOS'});
set(gca, 'fontsize', 13);
xlim([10 100]);
ylim([-110 -50]);
%
%
