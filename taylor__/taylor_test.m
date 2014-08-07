%the following code depicts the taylor line source method for antenna synthesis

% Program written by: Sathvik N. Prasad
% Date              : 18/06/2014 


len=input('Specify length of line source (in wavelengths)\n');
nbar=input('Specify number of side lobes of the same constant level (nbar)\n');
side_lobe=input('Specify desired constant side lobe level (in -dB)\n');

Ntheta=360;

theta=linspace(0,pi,Ntheta);
u=pi*len*cos(theta);

thres=-35;

R0=10^(-side_lobe/20);
A=1/pi*acosh(R0);
sig=nbar/sqrt(A^2+(nbar-0.5)^2);

n=1:nbar-1;

un = pi * sig * sqrt(A^2 +(n-.5).^2);

SF_rec=sinc(u/pi);   

for p=1:nbar-1,
SF_rec=SF_rec.*(1-(u/un(p)).^2)./(1-(u/(p*pi)).^2);
end;

SF_recdb=20*log10(abs(SF_rec)/max(abs(SF_rec)));
SF_recdb(SF_recdb<=thres)=thres;
%SF_recdb = SF_recdb - thres;

figure(1);
plot(theta*180/pi,SF_recdb,'linewidth',2); grid on;
xlim([0 180]);
xlabel('\theta (in degrees)');
ylabel('Normalized Space Factor (dB)');
title('Taylor method based on Tschebyscheff error');