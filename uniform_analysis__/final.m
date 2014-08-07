% This matlab code plots the various radiation pattern ( relative directivity)
% of uniform array antennas
% 
% The methods considered are:
%     1. Broadside
%     2. Ordinaary end-fire
%     3. Hansel Woodyard end-fire
%     4. Phase(scanning)
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

disp('UNIFORM LINEAR ARRAY');
disp('');
option_c = 1;
%disp('Enter the type of array:=');
if (option_c==1)
    disp('OPTION 1-  BROADSIDE');
    disp('OPTION 2 - ORDINARY ENDFIRE');
    disp('OPTION 3 - HANSEL WOODYARD');
    disp('OPTION 4 - SCANNING');
    option_d = input('Enter the type of uniform array:=');
    
    if (option_d==1)%broadside
        Nelem = input('Enter the number of elements for broadside:=');
        d = input('Spacing between elements in terms of lambda for broadside:=');
        beta=0;
        psi=k.*d.*cos(theta)+beta;
        psi3=k.*d.*cos(THETA)+beta;
        AF=sinc((Nelem.*psi./2)/pi)./sinc((psi./2)/pi);   
        AF3=sinc((Nelem.*psi3./2)/pi)./sinc((psi3./2)/pi);
    end
        
    if (option_d == 2)%ordinary end-fire
        Nelem = input('Enter the number of elements for ordinary end-fire:=');
        d = input('Spacing between elements in terms of lambda for ordinary end-fire:=');
        disp('');
        thmax=input('Value of Thete where maximum should occur (THETA = 0 OR 180 DEG.:)=');
        disp('');
        if abs(thmax)<eps,beta=-k*d;
        elseif abs(thmax-180)<eps,beta=k*d;
        end
        psi=k*d*cos(theta)+beta;
        psi3=k*d*cos(THETA)+beta;
        AF=sinc((Nelem.*psi./2)/pi)./sinc((psi./2)/pi);
        AF3=sinc((Nelem.*psi3./2)/pi)./sinc((psi3./2)/pi);
    end
    if option_d==3, % HANSEN-WOODYARD END-FIRE 
        Nelem = input('Enter the number of elements for Hansel-Woodyard end-fire:=');
        disp('');
        d = input('Spacing between elements in terms of lambda for Hansel-Woodyard end-fire:=');
        disp('');
        thmax=input('Value of Thete where maximum should occur (THETA = 0 OR 180 DEG.:)=');
        disp('');
        if abs(thmax)<eps,beta=-(k*d+pi/(Nelem));
        elseif abs(thmax-180)<eps,beta=k*d+pi/Nelem; 
        end
        psi=k*d*cos(theta)+beta;
        psi3=k*d*cos(THETA)+beta;
        AF=sinc((Nelem.*psi./2)/pi)./sinc((psi./2)/pi);
        AF3=sinc((Nelem.*psi3./2)/pi)./sinc((psi3./2)/pi);
    end
    
    if option_d==4      %scanning
        Nelem = input('Enter the number of elements for Scanning array:=');
        disp('');
        d = input('Spacing between elements in terms of lambda for Scanning array:=');
        disp('');
        thmax=input('Value of Thete where maximum should occur ( from THETA = 0 to 180 DEG.:)=');
        disp('');
        beta=-k*d*cos(thmax*pi/180);
        psi=k*d*cos(theta)+beta;
        psi3=k*d*cos(THETA)+beta;
        AF=sinc((Nelem.*psi./2)/pi)./sinc((psi./2)/pi);
        AF3=sinc((Nelem.*psi3./2)/pi)./sinc((psi3./2)/pi);   
    end;

        
      
        
end

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