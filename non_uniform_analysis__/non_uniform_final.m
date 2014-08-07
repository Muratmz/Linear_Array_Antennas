% This matlab code plots the various radiation pattern ( relative directivity)
% of non uniform array antennas
% 
% The methods considered are:
%     1. Binomial linear array
%     2. Dolph Chebyshev array
%By SATHVIK N. PRASAD
%DATE: 10/06/2014
%SERC, IISc


clc;
clear;
close all;
disc=181;
MM=disc;
NN=disc;
theta_low=0;
theta_up=180;
phi_low=0;
phi_up=360;

M=1800;
k=2*pi;
theta=linspace(0,pi,M+1);
theta3=linspace(theta_low*pi/180,theta_up*pi/180,MM+1);
phi3=linspace(phi_low*pi/180,phi_up*pi/180,NN+1);
[THETA,PHI]=meshgrid(theta3,phi3);
dtheta=pi/M;

disp('NON-UNIFORM LINEAR ARRAY');
disp('');


disp('OPTION 1 - Binomial linear array');
disp('OPTION 2 - Dolph Chebyshev array');
option_e=input('Enter the option:=');

if option_e==1, % BINOMIAL
Nelem=0;   
while (Nelem<1),   
Nelem=floor(input('NUMBER OF ELEMENTS ='));
end 
d=input('SPACING d BETWEEN THE ELEMENTS (IN WAVELENGTHS) =');
beta=0;
for i=1:M+1
[AF,Ncoef,Coef]=bin(theta,Nelem,d);
end
for i=1:MM+1
[AF3,Ncoef3,Coef3]=bin(THETA,Nelem,d);
end

   
elseif option_e==2, % DOLPH-TSCHEBYSCHEFF
Nelem=0;   
while (Nelem<1),   
Nelem=floor(input('NUMBER OF ELEMENTS ='));
end 
d=input('SPACING d BETWEEN THE ELEMENTS (IN WAVELENGTHS) =');
beta=0;
RdB=input('SIDE LOBE LEVEL (IN dB) =');
for i=1:M+1;
[AF,Ncoef,Coef]=tscheby(theta,Nelem,d,RdB);
end
for i=1:MM+1;
[AF3,Ncoef3,Coef3]=tscheby(THETA,Nelem,d,RdB);
end

end

Coef=Coef(1:Ncoef);
Ncoef=Coef(1:Ncoef)/Coef(Ncoef);

U=(abs(AF)./max(abs(AF))).^2;
Prad=2*pi*sum(U.*sin(theta).*dtheta);
D=4*pi*U/Prad;
DdB=10.*log10(D+eps);
Do=max(D);
DodB=max(DdB);
   
figure;
polar_dB(theta*180/pi,DdB-max(DdB),max(-60,6*floor(min(DdB)/6)),0,12,'-')
hold on
polar_dB(-theta*180/pi,DdB-max(DdB),max(-60,6*floor(min(DdB)/6)),0,12,'-')
title('Polar plot of Relative Directivity (0< \phi <360 degrees)','Fontsize',15)




%Spherical Plot3D
D3=4*pi*abs(AF3).^2/Prad;
D3dB=10.*log10(D3);
D3dB=D3dB-min(min(D3dB));
DD=D3dB;
disc=size(DD,1);
spherical_plot(DD,THETA,PHI,disc)
ss=title('3D Spherical plot of Directivity (dB)','Fontsize',15);

fff = 20*log10(AF);
figure;
plot(1:1800,fff(1:1800));