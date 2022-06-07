%
% Figure 4 - The transmit power needed to achieve the rate R¯ in the scenario
%            shown in Fig. 3, as a function of the distance d1
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
% compute the noise power in dBm, σ^2
%
sigma2dBm = -174 + 10*log10(B) + noiseFiguredB; 
sigma2 = db2pow(sigma2dBm);
%
% define the channel gain functions based on the 3GPP Urban Micro (UMi)
%
% pathloss_3GPP_LOS  = @(x) db2pow(-28.0 - 20*log10(fc) - 22.0*log10(x));
% pathloss_3GPP_NLOS = @(x) db2pow(-22.7 - 26*log10(fc) - 36.7*log10(x));
% 
% define the antenna gains at the source, Relay/IRS, and destination. 
% the numbers are in linear scale
%
Gs = db2pow(5); % antenna gains at the source
Gr = db2pow(5); % antenna gains at the relay/IRS
Gd = db2pow(0); % antenna gains at the destination
%
% set the amplitude reflection coefficient
%
alpha = 1;
%
% define the range of different number of reflection elements in the IRS
%
N = [25 50 100 150];
%
% set the range of rate values
%
Rbar = [4 6];
%
% define distances in simulation setup
%
d_SR = 80; % distance between the source and IRS/relay
dv = 10;   % minimum distance between destination and the IRS/relay
%
% define the range of d1 values in the simulation setup
%
d1range = 40:100;
%
% prepare to save simulation results
%
powerFractionRelay = zeros(length(d1range), length(Rbar));
Nmin = zeros(length(Rbar), 1);
%
% go through all rate values, which is 4 and 6 bit/s/Hz
%
for i = 1:length(Rbar)
    %
    % prepare to save transmit powers
    %
    P_IRS  = zeros(length(d1range), length(N));
    P_DF   = zeros(length(d1range), 1);
    P_SISO = zeros(length(d1range), 1);
    %
    % compute required SINR values
    %
    SINR = 2^(Rbar(i)) - 1;      % SISO and IRS
    SINR_DF = 2^(2*Rbar(i)) - 1; % DF relaying
    %
    % go through all values of d1
    %
    for k = 1:length(d1range)
        %
        % extract value of d1
        %
        d1 = d1range(k);
        %
        % compute distance between the source and destination
        %
        d_SD = sqrt(d1^2 + dv^2);
        % 
        % compute distance between the IRS/relay and destination
        %
        d_RD = sqrt((d1 - d_SR)^2 + dv^2);
        %
        % compute the channel gains using the 3GPP models and antenna gains
        %
        betaSR = pathloss_3GPP_LOS(d_SR, fc) * Gs * Gr;  % β_sr
        betaRD = pathloss_3GPP_LOS(d_RD, fc) * Gr * Gd;  % β_rd
        betaSD = pathloss_3GPP_NLOS(d_SD, fc) * Gs * Gd; % β_sd
        betaIRS = betaSR * betaRD;  % β_IRS = β_sr * β_rd
        %
        % compute the transmit power in mW in the SISO case, using Eq.(11)
        %
        % p_SISO = (2^R¯ - 1) * σ^2 / β_sd;
        P_SISO(k) = SINR * sigma2 / betaSD;
        %
        % compute the transmit power in mW in the IRS case, using Eq.(12)
        %
        % p_IRS(N) = (2^R¯ - 1) * σ^2 / (√β_sd + N*α*√β_IRS)^2;
        P_IRS(k, :) = SINR * sigma2 ./ (sqrt(betaSD) + N * alpha * sqrt(betaIRS)).^2;
        %
        % compute the transmit power in mW in the DF relaying case, using Eq.(14)
        %
        if betaSR >= betaSD
            % SINR_DF = 2^(2*R¯) - 1;
            bSum = betaSR + betaRD - betaSD;
            P_DF(k) = SINR_DF * sigma2 * bSum / (2 * betaRD * betaSR);
            powerFractionRelay(k, i) = 2 * (betaSR - betaSD) / bSum;
        else
            P_DF(k) = SINR_DF * sigma2 / betaSD;
            powerFractionRelay(k, i) = 0;
        end
        %
        % compute the number of reflecting elements needed to get a lower
        % transmit power with the IRS than with DF relaying, using Eq.(15)
        %
        if d1 == d_SR
            % ρ = p / σ^2, β_IRS = β_sr * β_rd
            p = P_DF(k); bIRS = betaIRS; bSum = betaSR + betaRD - betaSD;
            Nmin(i) = (sqrt((sqrt(1 + (2*p*bIRS / bSum / sigma2)) - 1)*(sigma2/p))...
                - sqrt(betaSD)) / (alpha * sqrt(bIRS));
        end
    end
    %
    % plot simulation results
    %
    figure;
    hold on; box on; grid on;
    %
    plot(d1range, 10*log10(P_SISO), 'k--', 'LineWidth', 2);
    plot(d1range, 10*log10(P_IRS(:,1)), 'r-', 'LineWidth', 2);
    plot(d1range, 10*log10(P_DF), 'b-.', 'LineWidth', 2);
    for n = 2:length(N)
        plot(d1range, 10*log10(P_IRS(:,n)), 'r-','LineWidth', 2);
    end
    % 
    title('Figure 4');
    xlabel('Distance d_1 [m]');
    ylabel('Transmit power [dBm]');
    %
    legend('SISO', 'IRS', 'DF relay', 'Location', 'NorthWest');
    set(gca, 'fontsize', 12);
    xlim([40 100]);
end
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
