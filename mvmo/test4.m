%this program takes the values of the current in each element of the linear
%array and the successive phase of the elements, for EVEN number of
%elements and then, calculates the array factor for the linear array
%antenna configuration



%case 1: linear array : 32 elements
%      uniformly separated by lambda/2 spacing
%fgure 7 in MVMO paper

clc;
clear all;
close all;

M = 32;
m = [ -1*M+1 : 2 : M-1];
m = m/2;
k = 2 *pi;
lambda = 1;
d = lambda/2;
kd= k*d;
samplerate = 1800;
theta = -pi:pi/samplerate:pi-pi/1800;
% u = kd*m*cos(theta)+del

magf = load('fig7_im.txt');
magf = magf';
magf = [fliplr(magf) magf];
phaf = load('fig7_degree.txt');
phaf = phaf';
phaf = [-1.*fliplr(phaf) phaf];

u = zeros(length(m),length(theta));
for ii = 1 :length(m)
    del(ii,:) = phaf(ii).*ones(1,length(theta)).*pi./180;
    %u(ii,:) = kd.*m(ii).*sin(theta)+del(ii,:);
 
end

for ii = 1 :length(m)
   
u(ii,:) = kd.*m(ii).*sin(theta)+del(ii,:);

end

af = zeros(1,length(theta));
for ii = 1 :length(m)
   
af = af + magf(ii).*exp(j*u(ii,:));

end
afdb = 20*log10(abs(af)/max(abs(af)));
plot(sin(theta),afdb);
hold on;

%case 2: linear array : 32 elements
%      uniformly separated by lambda/2 spacing
%fgure 8 in MVMO paper



M = 32;
m = [ -1*M+1 : 2 : M-1];
m = m/2;
k = 2 *pi;
lambda = 1;
d = lambda/2;
kd= k*d;
theta = -pi:pi/1800:pi-pi/1800;
% u = kd*m*cos(theta)+del

magf = load('fig8_im.txt');
magf = magf';
magf = [fliplr(magf) magf];
phaf = load('fig8_degree.txt');
phaf = phaf';
phaf = [-1.*fliplr(phaf) phaf];

u = zeros(length(m),length(theta));
for ii = 1 :length(m)
    del(ii,:) = phaf(ii).*ones(1,length(theta)).*pi./180;
    %u(ii,:) = kd.*m(ii).*sin(theta)+del(ii,:);
 
end

for ii = 1 :length(m)
   
u(ii,:) = kd.*m(ii).*sin(theta)+del(ii,:);

end

af = zeros(1,length(theta));
for ii = 1 :length(m)
   
af = af + magf(ii).*exp(j*u(ii,:));

end
afdb = 20*log10(abs(af)/max(abs(af)));
plot(sin(theta),afdb,'color','red');

a = -1:2/(2*samplerate):1-2/(2*samplerate);
plot(a, -25.*ones(1,2*samplerate));
plot(a, -35.*ones(1,2*samplerate))
plot(a, -40.*ones(1,2*samplerate))
hold off;

