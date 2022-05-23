%
% CCN HW2
%
close all
clc
clear
tic
T = 1;
pout_s = zeros(3,16);
pout_t = zeros(3,16);
X = 0:15;
Maxbit = 1e5;

SNR_db = 0:1:30;           % SNR range = [0, 30]
SNR_points = size(SNR_db); % number of points to be ploted

for x = 0:15
    i = 0;
    err1 = 0;
    err2 = 0;
    err3 = 0;
    while(i < Maxbit)
        h = 1 / sqrt(2) * ( randn(1,3) + 1j * randn(1,3) ); % Rayleigh
        a1 = ( abs(h(1)) ) ^ 2;
        a2 = ( abs(h(2)) ) ^ 2;
        a3 = ( abs(h(3)) ) ^ 2;
        SNR = 10 ^ ( x / 10 );
        
        
% N=1
        r1 = SNR * a1;
        if(r1 < T)
            err1 = err1 + 1;
        end
% N=2
        a12 = max(a1,a2);
        r2 = SNR * a12;
        if(r2 < T)
            err2 = err2 + 1;
        end
% N=3
        A = [a1,a2,a3];
        a13 = max(A);
        r3 = SNR * a13;
        if(r3 < T)
            err3 = err3 + 1;
        end
        i = i + 1;
    end
    pout_s(1,x+1) = err1 / Maxbit;
    pout_t(1,x+1) = ( T / SNR ) ^ 1;
    pout_s(2,x+1) = err2 / Maxbit;
    pout_t(2,x+1) = ( T / SNR ) ^ 2;
    pout_s(3,x+1) = err3 / Maxbit;
    pout_t(3,x+1) = ( T / SNR ) ^ 3;
    
    for i = 1:SNR_points(2)
        SNR_temp(i) = (10^(SNR_db(i) / 10));
        % 
        % No Direct Link 
        %
        G(i) = (1 / (a1^2 + 1 / SNR_temp(i)));
    end
end
figure();
line1 = semilogy(X,pout_s(1,:),'-x');
line1.Color = [1 0 0];
line1.LineWidth = 1.5;
hold on;
line2 = semilogy(X,pout_s(2,:),'-x');
line2.Color = [0 0 1];
line2.LineWidth = 1.5;
line3 = semilogy(X,pout_s(3,:),'-x');
line3.Color = [1 0 1];
line3.LineWidth = 1.5;
line4 = semilogy(X,pout_t(1,:),'-o');
line4.Color = [0 1 0];
line4.LineWidth = 1.5;
line5 = semilogy(X,pout_t(2,:),'-o');
line5.Color = [0 1 1];
line5.LineWidth = 1.5;
line6 = semilogy(X,pout_t(3,:),'-o');
line6.Color = [0 0 0];
line6.LineWidth = 1.5;

xlabel('SNR (dB)')
ylabel('P_{outage}')
legend('simulation N=1','simulation N=2','simulation N=3', ...
       'theoretical N=1','theoretical N=2','theoretical N=3')
grid on;
hold off;

figure();
line7 = semilogy(SNR_db, G(:),'-x');
line7.Color = [1 0 0];
line7.LineWidth = 1;
xlabel('SNR (dB)')
ylabel('Average squared gain G^2')
grid on;
toc




