%this function uses the pop value and plots the graph
%used to check the value of the solution got after the optimization
function  mvmo_fitness(pop)







M = 32;
m = [ -1*M+1 : 2 : M-1];
m = m/2;
k = 2 *pi;
lambda = 1;
d = lambda/2;
kd= k*d;
samplerate = 1800;
theta = -pi/2:pi/samplerate:pi/2-pi/samplerate;



magf = pop(1:M/2);
magf = [fliplr(magf) magf];



 phaf = pop(M/2+1:M);
phaf = [-1.*fliplr(phaf) phaf];


% u = kd*m*cos(theta)+del
u = zeros(length(m),length(theta));
for ii = 1 :length(m)
    
    del(ii,:) = phaf(ii).*ones(1,length(theta));
    
  
 
end
sintheta = sin(theta);
for ii = 1 :length(m)
   
u(ii,:) = kd.*m(ii).*sintheta+del(ii,:);

end

af = zeros(1,length(theta));
for ii = 1 :length(m)
   
af = af + magf(ii).*exp(j*u(ii,:));

end
afdb = 20*log10(abs(af)/max(abs(af)));
plot(sin(theta),afdb);


             
