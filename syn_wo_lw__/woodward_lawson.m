%The following code illustrates the use of woodward lawson method for 
%antenna synthesis
%% Program written by: Sathvik N. Prasad
% Date              : 25/06/2014 
clc;
clear;
close all;
P=input('enter number of elements:');
la = 1; 
ss = input('enter the spacing innterms of lambda:');
d=la*ss; %spacing btw elements
ki=2*pi/la; % 2*pi/lambda
M=(P*d)/la; %total length
m=-M:M;  % array to represent the location of each element
wn=m*(la/(P*d)); %representation of arccos to 
wn=acos(wn);     %find the sampling location
%balanis 7.23

deg=wn*180/pi; %convert theta_m to degrees
%disp(length(tt));
%af =[ zeros(1,45) ones(1,90) zeros(1,45)  ]; %required pattern/array factor
 tt = 0:pi/180:pi-pi/180;
 af = sin(tt);




p = ceil(deg)+1; 
for ii = 1:length(deg)
    if p(ii)>=180
        p(ii) = 180;
    end
end

for ii = 1:length(deg) %get the samples
an(ii)=af(p(ii)); 
end
k=length(an); %length of sampled data
th=0:0.02:pi; 
th1 = th/pi*180;
afb=0;

figure(1);
subplot(1,2,1);
%plot(th,era,'ro'); 
plot(1:length(af),af, 'color','r');
grid on
xlabel('\theta in degrees');
ylabel('amplitude');
title('Required patterm');
ylim([0 1.2]);

subplot(1,2,2);
%plot(th,era,'ro'); 
plot(1:length(af),20*log10(abs(af/max(abs(af)))), 'color','r');
grid on
xlabel('\theta in degrees');
ylabel('amplitude in dB (Normalized)');
title('Required patterm');




figure(2);
hold on 
%array factor eq. balanis 7.21
%add array factor values at the sampling location to show the final graph
yy = zeros(1,158);
for lk=1:k 
    
afb=afb+an(lk)*(sin((P/2)*ki*d*(cos(th)-cos(wn(lk))))./(P*sin((1/2)*ki*d*(cos(th)-cos(wn(lk)))))); 
 yy =  afb + yy;

plot(th,double(afb)) 
 afb =0;

end 
xlabel('\theta in radians');
ylabel('amplitude- normalised');
title('Plot of various composing functions');

y=abs(yy);
y=double(y);
for ii = 1:length(th1)
    if th1(ii)>=180
        th1(ii)= 180;
    end
    if th1(ii)==0
        th1(ii) = 1;
    end
end



figure(3);
%obtained plot in black
plot(th,y,'color','black','LineWidth',2);

hold on
%figure(1);
%plot(th,era,'ro');
%desired plot in red
kk = 1:length(af);
kk = kk .* pi/180;
plot(kk,af, 'color','red','LineWidth',2);
grid on
xlabel('\theta in degrees');
ylabel('amplitude');
%title('Required patterm');


hold off
grid on 
xlabel('\theta in Radians','fontsize',18); 

ylabel('Magnitude','fontsize',18); 
title('Synthesised array pattern using Woodward-Lawson Method (N=15 , d = 0.4 \lambda)','fontsize',18); 
xlim([0 pi])



%plot in dB



figure(4);
%obtained plot in black
plot(th,20*log10(abs(y/max(abs(y)))),'color','black');

hold on

%desired plot in red
kk = 1:length(af);
kk = kk .* pi/180;
plot(kk,20*log10(abs(af/max(abs(af)))), 'color','red');
grid on
