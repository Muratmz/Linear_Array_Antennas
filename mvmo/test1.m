

clc;
clear;
close all;

theta = -pi:pi/1800:pi-pi/1800;

u = sin(theta);
size(u)
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

plot(u,af);
ylim([-50 5])

% afn = af./(max(af));
% for ii = 1: length(af)
%     if(afn(ii)==0)
%         afn(ii) = afn(ii) + 0.0000000001;
%     end
% end
% figure(1);
% plot(u,afn)
% figure(2);
% afndb = 10*log10(afn);
% afndbsl = afndb;
% for ii = 1:length(u)
%     if(u(ii)>-1 && u(ii)<0)
%         afndbsl(ii) = -25;
%     elseif (u(ii)>0.6 && u(ii)< 1)
%         afndbsl(ii) = -25;
%     end
%        
% end
% 
% plot(u,afndbsl)

