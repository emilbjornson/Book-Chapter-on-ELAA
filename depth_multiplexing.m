close all
clear;

%Carrier frequency
f_c = 3e9;

%Wavelength
lambda = 3e8/f_c;

%Antenna spacing in fraction of wavelengths
scalefactor = sqrt(2)/4;

%Diagonal of each receive antenna
D_antenna = scalefactor*lambda;

%Number of receive antennas per dimension
Ndim = 100;

%Total number of receive antennas
N = Ndim^2;

%Compute the Fraunhofer distance
fraunhoferDistanceAntenna = 2*D_antenna^2/lambda;

%Compute the Fraunhofer array distance
fraunhoferArrayDistance = N*fraunhoferDistanceAntenna;

%Define the focal points of the receivers
focalPoints = [fraunhoferArrayDistance fraunhoferArrayDistance/20 fraunhoferArrayDistance/40 fraunhoferArrayDistance/60 fraunhoferArrayDistance/80];

%Number of users
K = length(focalPoints);

%Prepare to save results
H = zeros(N,K);

%Diagonal of entire surface
D_RIS = Ndim*D_antenna;


%Center points of the antennas in the x and y dimensions
gridPoints = (-(Ndim-1)/2:1:(Ndim-1)/2)*D_antenna/sqrt(2);

%Go through all distances
for k = 1:K

    disp(['User ' num2str(k) ' out of ' num2str(K)]);

    %Extract the distance to receiver
    z = focalPoints(k);

    E_fun = @(x,y) sqrt((z.*(x.^2+z.^2))./((x.^2+y.^2+z.^2).^(5/2))).*exp(-1j*(2*pi).*(sqrt(x.^2+y.^2+z.^2))/lambda);
    %E2_fun = @(x,y) (z.*(x.^2+z.^2))./((x.^2+y.^2+z.^2).^(5/2));

    %Compute all the terms in channel expression
    numerator_exact = zeros(Ndim,Ndim);
    %denominator_exact =  integral2(E2_fun, -D_antenna/sqrt(8),D_antenna/sqrt(8),-D_antenna/sqrt(8),D_antenna/sqrt(8));


    for xdim = 1:Ndim

        for ydim = 1:Ndim

            numerator_exact(xdim,ydim) = integral2(E_fun, gridPoints(xdim)-D_antenna/sqrt(8),gridPoints(xdim)+D_antenna/sqrt(8),gridPoints(ydim)-D_antenna/sqrt(8),gridPoints(ydim)+D_antenna/sqrt(8));

        end

    end

    %H(:,k) = sqrt((2/D_antenna^2)/(N^2*denominator_exact))*numerator_exact(:);
    H(:,k) = sqrt(2/D_antenna^2)*numerator_exact(:);

end


%Normalize the channel matrix to have SNR 0 dB at d_FA
gainAtdFA = norm(H(:,1)).^2;
H = H / sqrt(gainAtdFA);

%Compute the Grammian of the channel matrix
HH = H'*H;


%Compute zero forcing with normalized columns
W_ZF = H/HH;
for k = 1:K
    W_ZF(:,k) = W_ZF(:,k)/norm(W_ZF(:,k));
end

%Extract the channel gains after ZF
HW = abs(diag(H'*W_ZF)).^2;




%Select the range of SNR values at d_FA
SNRdB = -20:20;
SNR = db2pow(SNRdB);


%Compute the sum rate with ZF
sumrate_ZF = zeros(length(SNR),1);

for s = 1:length(SNR)
    
    %Compute the waterfilling power allocation
    poweralloc = waterfill(SNR(s),1./HW');
    
    %Compute the sum rate
    sumrate_ZF(s) = sum(log2(1+poweralloc.*HW'));

end


%Benchmark where the users take turns in round robin fashion
sumrate_scheduling = sum(log2(1+SNR'*sum(abs(H).^2,1)),2)/K;



%Plot simulation results
set(groot,'defaultAxesTickLabelInterpreter','latex');

figure;
hold on; box on;
plot(SNRdB,sumrate_ZF,'b-','LineWidth',2);
plot(SNRdB,sumrate_scheduling,'k--','LineWidth',2);
xlabel('SNR at $d_{\mathrm{FA}}$ [dB]','Interpreter','latex');
ylabel('Sum SE [bit/s/Hz]','Interpreter','latex');
set(gca,'fontsize',18);
legend({'ZF','Scheduling'},'Interpreter','latex','Location','NorthWest');
grid on;
ylim([0 70]);



