% this program takes the number of elements as one of the
% input and also the shaping curve i.e. array factor as an input that is hard coded in the program
% it uses the fourier transform synthesis method to calculate the excitation coeff and 
% hence the array factor
% the plots are of both current values and the array factor
% 
% Program written by: Sathvik N. Prasad
% Date              : 18/06/2014 

clc;
close all; 
clear all;
Nel=input('Enter number of elements in the array (odd number):=');
    warning off;
%%
 t  = [-pi: pi/180:pi-pi/180 ] ;
 SFAF = sin(t);
%SFAF = [ zeros(1,45) ones(1,90) zeros(1,45)];


%%
d=0.5;
Ntheta=180;
theta=linspace(0,pi,Ntheta);
psi=pi*cos(theta);
dth=theta(2)-theta(1);
M=(Nel-1)/2;
m=-M:M;
 for ind=1:length(m),
    a(ind)=1/2*sum(SFAF.*exp(-j*m(ind)*pi*cos(theta)).*sin(theta)*dth);
 end;
 SFAF_rec=zeros(size(theta));
 for ind=1:length(m),
    SFAF_rec=SFAF_rec+a(ind)*exp(j*m(ind)*psi);
 end;               
%  p = length(theta)
  q = length(a)
figure(1);
    plot(theta*180/pi,abs(SFAF),'color','b'); hold on;
    
    
    
    %%
    plot(theta*180/pi,abs(SFAF_rec),'color','r');
    
    legend('Desired','Linear array'); grid on;
    
    xlabel('\theta (in degrees)');
    ylabel('Array Factor')
    title(' Synthesis using Fourier method');
    
    figure(2);
    s = 1: length(a);
    stem(s,a/abs(max(a))); 
    ylim([-1.1 1.1])
    hold on;
    k = zeros(1,Nel);
    plot(s,k,'color','r');
     xlabel('\theta (in degrees)');
    ylabel('Array Factor')
    title(' Synthesis using Fourier method');
   ;
    
