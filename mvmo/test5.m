%figure 6 of MVMO paper ---- plots the 20 element array pattern clc;
clear all;
close all;

M = 10;
m = [ -1*M+1 : 2 : M-1];
m = m/2;
k = 2 *pi;
lambda = 1;
d = lambda/2;
kd= k*d;
samplerate = 1800;
theta = -pi:pi/samplerate:pi-pi/samplerate;
% u = kd*m*cos(theta)+del

magf = load('fig6_im.txt');
magf = magf';
magf = [fliplr(magf) magf];
phaf = load('fig6_degree.txt');
phaf = phaf';
phaf = [-1.*fliplr(phaf) phaf];
d = load('fig6_d.txt');
d = d' .* lambda/2;
d = [-1*fliplr(d) d];

u = zeros(length(m),length(theta));
for ii = 1 :length(m)
    del(ii,:) = phaf(ii).*ones(1,length(theta)).*pi./180;
 
end

for ii = 1 :length(m)
   
u(ii,:) = kd.*d(ii).*sin(theta)+del(ii,:);

end

af = zeros(1,length(theta));
for ii = 1 :length(m)
   
af = af + magf(ii).*exp(j*u(ii,:));

end
afdb = 20*log10(abs(af)/max(abs(af)));
%plot([-90:180/length(afdb):90-90/length(afdb)],afdb);
plot(1:length(afdb),afdb)
