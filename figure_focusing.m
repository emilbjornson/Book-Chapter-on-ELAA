% This code corresponds to Figure 11 of the book chapter entitled "Near-Field Beamforming and Multiplexing Using Extremely Large Aperture Arrays"
%written by Parisa Ramezani and Emil Björnson. 

close all
clear;

%Carrier frequency
f_c = 3e9;

%Wavelength
lambda = 3e8/f_c;

%Side length of each receive antenna in fraction of wavelengths
sidelength = 1/4; 

%Diagonal of each receive antenna
D_antenna = sqrt(2)*sidelength*lambda;

%Number of receive antennas per dimension
Ndim = 100;

%Total number of receive antennas
N = Ndim^2;

%Compute the Fraunhofer distance
fraunhoferDistanceAntenna = 2*D_antenna^2/lambda;

%Björnson distance in number of Fraunhofer distance of a single antenna
distance_B = (2*D_antenna*Ndim)/(fraunhoferDistanceAntenna);

%Define the range of points along the horizontal axis
relativeRange = sort([logspace(1,5,300) distance_B N/10 ]);

%Determine the range of distances to be considered
zRange = relativeRange*fraunhoferDistanceAntenna;


%Prepare to save results
channelVectors = zeros(N,length(relativeRange));


%Diagonal of entire array
D_array = Ndim*D_antenna;


%Center points of the antennas in x and y dimensions
gridPoints = (-(Ndim-1)/2:1:(Ndim-1)/2)*D_antenna/sqrt(2);

%Go through all distances
for m = 1:length(zRange)
    
    disp(['Distance ' num2str(m) ' out of ' num2str(length(zRange))]);
    
    %Extract distance to transmitter
    z = zRange(m);
    
    %Define E(x,y) and its absolute value squared
    %(we have removed E_0/sqrt(4*pi) since it will cancel out)
    E_fun = @(x,y) sqrt((z.*(x.^2+z.^2))./((x.^2+y.^2+z.^2).^(5/2))).*exp(-1j*(2*pi).*sqrt((x.^2+y.^2+z.^2))/lambda);
    E2_fun = @(x,y) (z.*(x.^2+z.^2))./((x.^2+y.^2+z.^2).^(5/2));
    
    
    %%Array Gain
    
    %compute the integral in the demominator
    denominator_exact =  integral2(E2_fun, -D_antenna/sqrt(8),D_antenna/sqrt(8),-D_antenna/sqrt(8),D_antenna/sqrt(8));
    
    %compute all the terms in the numerator
    numerator_exact = zeros(Ndim,Ndim);
    
    
    for xdim = 1:Ndim
        
        for ydim = 1:Ndim
            
            numerator_exact(xdim,ydim) = integral2(E_fun, gridPoints(xdim)-D_antenna/sqrt(8),gridPoints(xdim)+D_antenna/sqrt(8),gridPoints(ydim)-D_antenna/sqrt(8),gridPoints(ydim)+D_antenna/sqrt(8));
            
        end
        
    end
    
    channelVectors(:,m) = sqrt((2/D_antenna^2)/(N*denominator_exact))*numerator_exact(:);
    
end


%Determine which points the array focuses at
focusIndex = [find(relativeRange==distance_B) find(relativeRange==N/10) ]; 


%Determine the normalized antenna array gains that is achieved at different
%distances
arrayGain = zeros(length(relativeRange),length(focusIndex));

for n = 1:length(focusIndex)
    
    combining = channelVectors(:,focusIndex(n))/norm(channelVectors(:,focusIndex(n)));
    
    arrayGain(:,n) = abs(channelVectors'*combining).^2;
    
end


matchedfilteringGains = sum(abs(channelVectors).^2,1)';
arrayGainFarfield = abs(sum(channelVectors,1)').^2/N;

set(groot,'defaultAxesTickLabelInterpreter','latex');


%% Plot the simulation results
figure;
hold on; box on; grid on;
plot(relativeRange, matchedfilteringGains,'k:', 'Linewidth', 2);
plot(relativeRange, arrayGain(:,1),'k-', 'Linewidth', 2);
plot(relativeRange, arrayGain(:,2),'b--', 'Linewidth', 2);
plot(relativeRange, arrayGainFarfield,'r-.', 'Linewidth', 2);
set(gca,'XScale','log');
xticks([10 100 1000 10000 100000])
xticklabels({'$10 d_F$','$10^2 d_F$','$10^3 d_F$','$10^4 d_F$','$10^5 d_F$'})
ylim([0 1]);
yticks([0 0.25 0.5 0.75 1]);
legend({'Maximum gain', 'Focus at $z=d_B$','Focus at $z=d_{FA}/10$','Focus at $z=\infty$'},'Interpreter','Latex','Location','southeast');
xlabel('Propagation distance','Interpreter','Latex');
ylabel('Normalized antenna array gain','Interpreter','Latex');
set(gca,'fontsize',18);
axis tight
