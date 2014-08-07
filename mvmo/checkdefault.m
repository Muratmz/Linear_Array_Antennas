%%

%plots the graph of the MVMO output in the paper by taking the values from
%the text files

function afdb = checkdefault()
%%
M = 32;
m = [ -1*M+1 : 2 : M-1];
m = m/2;
k = 2 *pi;
lambda = 1;
d = lambda/2;
kd= k*d;
samplerate = 180;
theta = -pi/2:pi/samplerate:pi/2-pi/samplerate;
% u = kd*m*cos(theta)+del8
magf = load('fig7_im.txt');
magf = magf';
magf = [fliplr(magf) magf];
%phaf = load('fig7_degree.txt');
phaf = load('fig7_rad.txt');
phaf = phaf';
phaf = [-1.*fliplr(phaf) phaf];

u = zeros(length(m),length(theta));
for ii = 1 :length(m)
    %uses phase in degrees from the file
    del(ii,:) = phaf(ii).*ones(1,length(theta)).*pi./180;
    %uses phase in radians from the file
    del(ii,:) = phaf(ii).*ones(1,length(theta));
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
%plot(sin(theta),afdb);