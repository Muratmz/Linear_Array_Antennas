del = load('fig7_degree.txt');
im = load('fig7_im.txt');
theta = -pi : pi/180: pi-pi/180;

del = [del del];
im = [im im];
m1 = [-16:-1];
m2 = [1:16];
m = [m1 m2]
lambda =1;

for ii = 1: length(theta)
    for jj = m
        af(ii) = af(ii) + im(jj+16)*exp(i*(2*pi/lambda*lambda*sin(theta(ii) + del(jj+16)*pi/180 )));
        
    end
    
    
end



