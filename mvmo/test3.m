clc;
clear;
close all;

%%
% theta = 0:pi/180:pi-pi/180;
% afd = sin(theta);
% %length(afd)
% afd = [ .1*ones(1,90) afd .1*ones(1,90)];
% plot(1:length(afd),afd);
% psi=k*d*cos(theta)+beta;
% AF=sinc((Nelem.*psi./2)/pi)./sinc((psi./2)/pi);
%%

lambda = 1;
d = lambda/2;
M = 32;
m = M/2;
theta = 0:pi/1800:pi-pi/1800;
%%
magf = load('fig7_im.txt');
phaf = load('fig7_degree.txt');
phaf = phaf.*pi/180;
magf = ones(1,16);
% phaf = zeros(1,16);

%am = exp(j*phaf) .* magf;
am = magf;
am = am';


%%

m = [1:2:M-1];
u = pi*d/lambda*sin(theta);
costerm = zeros(length(m),length(u));

for ii = 1 :length(m)
   pha(ii,:) = phaf(ii).*ones(1,length(theta));
    
end

for ii = 1:length(m)
costerm(ii,:) = cos(m(ii).*u + pha(ii,:));
end

af = zeros(1,length(u));
for ii = 1:length(m)

    af = af + costerm(ii,:);
    
end

 afdb = 20*log10(abs(af)/max(abs(af)));
% %plot(theta,af);
 plot(theta,afdb);

%%
% NN = 181;
% theta=linspace(0,pi,M+1);
% phi3=linspace(0*pi/180,360*pi/180,NN+1);
% k=2*pi;
% 
% Nelem = 32;
% psi=k.*d.*cos(theta);%+beta;
% 
% AF=sinc((Nelem.*psi./2)/pi)./sinc((psi./2)/pi)