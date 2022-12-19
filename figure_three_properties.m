% This code corresponds to Figure 6 of the book chapter entitled "Near-Field Beamforming and Multiplexing Using Extremely Large Aperture Arrays"
%written by Parisa Ramezani and Emil Björnson. 

close all;
clear;

%Wavelength
lambda = 0.1;

%Propagation distance
d = 20;

%Number of antennas
N = logspace(3,9,100);

%Area of each antenna
A = (lambda/4)^2;

%Compute the diagonal of the array for different number of antennas
diagonals = sqrt(2*N*A);

%Compute free-space far-field channel gain
beta_d = A/(4*pi*d^2);

%Compute exact total channel gain
PrxPtx_exact = (N*beta_d)./(3*(N*beta_d*pi+1).*sqrt(2*N*beta_d*pi+1)) + 2/(3*pi)*atan(N*beta_d*pi./sqrt(2*N*beta_d*pi+1));


%Compute total channel gains using integrals by considering all or some
%of the near-field properties

fun1 = @(x,y) 1./(4*pi*(x.^2 + y.^2 + d^2));
fun2 = @(x,y) d./sqrt( x.^2 + y.^2 + d^2) .* 1./(4*pi*(x.^2 + y.^2 + d^2));
fun3 = @(x,y) d./sqrt( x.^2 + y.^2 + d^2) .* (x.^2 + d^2)./(x.^2 + y.^2 + d^2) .* 1./(4*pi*(x.^2 + y.^2 + d^2));

PrxPtx_property1 = zeros(1,length(N));
PrxPtx_properties12 = zeros(1,length(N));
PrxPtx_properties123 = zeros(1,length(N));

for n = 1:length(N)
    
    PrxPtx_property1(n) = integral2(fun1,-sqrt(A*N(n))/2,sqrt(A*N(n))/2,-sqrt(A*N(n))/2,sqrt(A*N(n))/2);
    PrxPtx_properties12(n) = integral2(fun2,-sqrt(A*N(n))/2,sqrt(A*N(n))/2,-sqrt(A*N(n))/2,sqrt(A*N(n))/2);
    PrxPtx_properties123(n) = integral2(fun3,-sqrt(A*N(n))/2,sqrt(A*N(n))/2,-sqrt(A*N(n))/2,sqrt(A*N(n))/2);
    
end




%% Plot the simulation results
set(groot,'defaultAxesTickLabelInterpreter','latex');  

figure;
hold on; box on;
plot(diagonals,N*beta_d,'k:','LineWidth',2);
plot(diagonals,PrxPtx_property1,'b--','LineWidth',2);
plot(diagonals,PrxPtx_properties12,'k-.','LineWidth',2);
plot(diagonals,PrxPtx_exact,'r','LineWidth',2);
xlabel('Diagonal $\sqrt{2NA}$ of the array [m]','Interpreter','Latex');
ylabel('Channel gain','Interpreter','Latex');
legend({'Far-field approx','Property 1','Properties 1,2','Exact (Properties 1,2,3)'},'Interpreter','Latex','Location','NorthWest');
set(gca,'fontsize',18);
xlim([0 50]);
ylim([0 0.25]);
