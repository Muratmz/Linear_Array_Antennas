
function val =  fin_fig_test(pop)
%%
% a = load('fig7_im.txt');
% a = a';
% b = load('fig7_degree.txt');
% b = b';
% pop = [a b];






%%
M = 32;
m = [ -1*M+1 : 2 : M-1];
m = m/2;
k = 2 *pi;
lambda = 1;
d = lambda/2;
kd= k*d;
samplerate = 1800;
theta = -pi:pi/samplerate:pi-pi/1800;


% magf = load('fig7_im.txt');
% magf = magf';
magf = pop(1:M/2);
magf = [fliplr(magf) magf];


%  phaf = load('fig7_degree.txt');
%  phaf = phaf';
 phaf = pop(M/2+1:M);
phaf = [-1.*fliplr(phaf) phaf];


% u = kd*m*cos(theta)+del
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

%%


theta = -pi:pi/1800:pi-pi/1800;

u = sin(theta);

jj = 1;
for ii = 1 :length(u)
    
    if(u(ii)>=-1 && u(ii)<=0)
        af(ii) = -25;
    elseif (u(ii)>0 && u(ii)<0.1)
        af(ii) = 0;
    elseif (u(ii)>=0.1 && u(ii)<=0.5)
        af(ii) = 1/sin(u(ii));
        aftn(jj) = af(ii);
        jj = jj +1;
        
    elseif (u(ii)>0.5 &&u(ii)<0.6)
        af(ii) = -13.54;
    else af(ii) = -25;
    end
   
            
%     if(sdvs)
%         sdvVD
%     else if(u(ii)>=.1 && u(ii)<=.5)
%         af(ii) = 1/sin(u(ii));
%     else if (u(ii)>=-1 && u(ii) <= -.35)
%         af(ii) = -25
%     
%     end
    
end

for ii = 1 :length(u)
if(u(ii)>=0.1 && u(ii)<=0.5)
        af(ii) = 20*log10(af(ii)/max(abs(aftn)));
        
end 
end





