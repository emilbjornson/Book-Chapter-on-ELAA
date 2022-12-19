% For Figure 10 of the Book Chapter
close all
clear;

%Carrier frequency
f_c = 3e9;

%Wavelength
lambda = 3e8/f_c;

%Changing this factor will not affect the result since everything becomes
%scaled accordingly
scalefactor = sqrt(2)/4;

%Diagonal of each receive antenna
D_antenna = scalefactor*lambda;

%Number of receive antennas per dimension
Ndim = 100;

%Total number of receive antennas
N = Ndim^2;

%Diagonal of entire array
D_array = Ndim*D_antenna;

%Center points of the antennas in x and y dimensions
gridPoints = (-(Ndim-1)/2:1:(Ndim-1)/2)*D_antenna/sqrt(2);

relativeRange = logspace(0,4.3,300);


%Compute the Fraunhofer distance
fraunhoferDistanceAntenna = 2*D_antenna^2/lambda;

%Determine the range of distances to be considered
zRange = relativeRange*fraunhoferDistanceAntenna;


%Prepare to save results
G_exact = zeros(length(relativeRange),1);
G_bound = zeros(length(relativeRange),1);

%Go through all distances
for m = 1:length(zRange)
        
    disp(['Distance ' num2str(m) ' out of ' num2str(length(zRange))]);
        
    %Extract distance to transmitter
    z = zRange(m);
        
    %Define the integrand E(x,y) in (5) and its absolute value square
    %(we have removed E_0/sqrt(4*pi) since it will cancel out)
    E_fun = @(x,y) sqrt((z.*(x.^2+z.^2))./((x.^2+y.^2+z.^2).^(5/2))).*exp(-1j*(2*pi).*sqrt((x.^2+y.^2+z.^2))/lambda);
    E2_fun = @(x,y) (z.*(x.^2+z.^2))./((x.^2+y.^2+z.^2).^(5/2));
        
    %Define the integrand described in the text on Page 4 and its absolute
    %value square (we have removed E_0/sqrt(4*pi) since it will cancel out)
    E_fun2 = @(x,y) sqrt(1./((x.^2+y.^2+z.^2))).*exp(-1j*(2*pi).*sqrt((x.^2+y.^2+z.^2))/lambda);
    E2_fun2 = @(x,y) (1./((x.^2+y.^2+z.^2)));
        
        
    numerator_exact = zeros(Ndim,Ndim);
    denominator_exact =  integral2(E2_fun, -D_antenna/sqrt(8),D_antenna/sqrt(8),-D_antenna/sqrt(8),D_antenna/sqrt(8));
    
    for xdim = 1:Ndim
            
        for ydim = 1:Ndim
                
            numerator_exact(xdim,ydim) = abs(integral2(E_fun, gridPoints(xdim)-D_antenna/sqrt(8),gridPoints(xdim)+D_antenna/sqrt(8),gridPoints(ydim)-D_antenna/sqrt(8),gridPoints(ydim)+D_antenna/sqrt(8))).^2;
                 
        end
            
    end
        
    G_exact(m) = (2/(N*D_antenna^2))*sum(numerator_exact(:))/denominator_exact;
        
    alpha_z = D_antenna^2/(8*z^2);
    
    numerator_bound = N*alpha_z/(2*(N*alpha_z+1)*sqrt(2*N*alpha_z+1)) + atan(N*alpha_z/sqrt(2*N*alpha_z+1));
    denominator_bound = alpha_z/(2*(alpha_z+1)*sqrt(2*alpha_z+1)) + atan(alpha_z/sqrt(2*alpha_z+1));
        
    G_bound(m) = numerator_bound/(N*denominator_bound);
    
end

set(groot,'defaultAxesTickLabelInterpreter','latex');

%% Plot the simulation results
figure;
semilogx(relativeRange, G_bound,'r-.', 'Linewidth', 2);
hold on; grid on;
semilogx(relativeRange, G_exact,'b-', 'Linewidth', 2);
xlim([1 10^4.3]);
xticks([1 10 100 1000 10000])
xticklabels({'$d_F$','$10 d_F$','$100 d_F$','$1000 d_F$','$10000 d_F$'})
ylim([0 1]);
yticks([0 .2 .4 .6 .8 1])
legend({'Upper bound','Exact'},'Interpreter','Latex','Location','SouthEast');
xlabel('Propagation distance','Interpreter','Latex');
ylabel('Normalized antenna array gain','Interpreter','Latex');
set(gca,'fontsize',15);


% For the small figure which zooms in a particular range of z
figure
indexOfInterest = (relativeRange < 1.2 * 10^4) & (relativeRange > 200);
semilogx(relativeRange(indexOfInterest) , G_bound(indexOfInterest),'r-.', 'Linewidth', 2)
hold on; grid on;
semilogx(relativeRange(indexOfInterest) , G_exact(indexOfInterest),'b-', 'Linewidth', 2)
set(gca,'fontsize',25);
xticks ([200*sqrt(2) 1000 4000 10000])
xticklabels({'$d_B$','$1000 d_F$','$4000 d_F$','$d_{FA}$'})
yticks([.92 .94 0.96 .98 1])
set(gca, 'XLimSpec', 'Tight');

