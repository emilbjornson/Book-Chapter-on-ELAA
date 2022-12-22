% This code corresponds to Figure 8 of the book chapter entitled "Near-Field Beamforming and Multiplexing Using Extremely Large Aperture Arrays"
%written by Parisa Ramezani and Emil Björnson.

close all;
clear;

%Wavelength
lambda = 0.1;

%Propagation distance
d = 20;

%Number of receive antennas
N = logspace(0,10,100);

%Area of each antenna
A = (lambda/4)^2;

%Compute the diagonal of the array for different number of antennas
beta_d = A/(4*pi*d^2);


%Compute exact total channel gain
PrxPtx_exact = (N*beta_d)./(3*(N*beta_d*pi+1).*sqrt(2*N*beta_d*pi+1)) + 2/(3*pi)*atan(N*beta_d*pi./sqrt(2*N*beta_d*pi+1));


%Compute SNR scaling to get 0 dB SNR with N = 1
SNR0scaling = 1/min(PrxPtx_exact);


%Define the scaling exponents
rho = [0 1/2 1];


%Compute the SNRs with the different power scaling laws
SNR = zeros(length(PrxPtx_exact),length(rho));

for j = 1:length(rho)
    
    SNR(:,j) = (SNR0scaling./N.^rho(j)).*PrxPtx_exact;
    
end


%% Plot the simulation results
set(groot,'defaultAxesTickLabelInterpreter','latex');  

figure;
hold on; box on; grid on;
plot(N,SNR(:,1),'r','LineWidth',2);
plot(N,SNR(:,2),'b--','LineWidth',2);
plot(N,SNR(:,3),'k-.','LineWidth',2);
set(gca,'XScale','log');
set(gca,'YScale','log');
xlabel('Number of antennas ($N$)','Interpreter','Latex');
ylabel('SNR','Interpreter','Latex');
legend({'$\rho=0$','$\rho=1/2$','$\rho=1$'},'Interpreter','Latex','Location','SouthWest');
set(gca,'fontsize',18);
ylim([1e-4 1e7])
xticks(10.^(0:2:10));
yticks(10.^(-4:2:6));

