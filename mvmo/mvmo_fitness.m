function val =  mvmo_fitness(pop)



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

%maximum side-lobe level
msll = -25;



magf = pop(1:M/2);
magf = [fliplr(magf) magf];



 phaf = pop(M/2+1:M);
phaf = [-1.*fliplr(phaf) phaf];


% u = kd*m*cos(theta)+del
u = zeros(length(m),length(theta));
for ii = 1 :length(m)
    
    %this uses the phase values in defgrees from pop.
    %del(ii,:) = phaf(ii).*ones(1,length(theta)).*pi./180;
    
    %this uses phase values in radians
    del(ii,:) = phaf(ii).*ones(1,length(theta));
    
    
    %u(ii,:) = kd.*m(ii).*sin(theta)+del(ii,:);
 
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
%plot(sin(theta),afdb);

%%


% theta = -pi:pi/1800:pi-pi/1800;
% 
% u = sin(theta);
% 
% jj = 1;
% for ii = 1 :length(u)
%     
%     if(u(ii)>=-1 && u(ii)<=0)
%         af(ii) = -25;
%     elseif (u(ii)>0 && u(ii)<0.1)
%         af(ii) = 0;
%     elseif (u(ii)>=0.1 && u(ii)<=0.5)
%         af(ii) = 1/sin(u(ii));
%         aftn(jj) = af(ii);
%         jj = jj +1;
%         
%     elseif (u(ii)>0.5 &&u(ii)<0.6)
%         af(ii) = -13.54;
%     else af(ii) = -25;
%     end
%    
%             
%     
% end
% 
% for ii = 1 :length(u)
% if(u(ii)>=0.1 && u(ii)<=0.5)
%         af(ii) = 20*log10(af(ii)/max(abs(aftn)));
%         
% end 
% end
afdbref = checkdefault();
val = 0;
for ii = 1:length(theta)
    
    
    
   if(sintheta(ii)>=0.1 && sintheta(ii)<=0.5)
       val = val + (abs(afdb(ii))-abs(afdbref(ii))).^2;
   %elseif (sintheta(ii)>=0.5892 && sintheta(ii)<=-0.0226 && afdb(ii)>=msll)
   elseif (sintheta(ii)>0.6 && sintheta(ii)<0 && afdb(ii)<msll)
       val = val + (abs(afdb) - abs(msll)).^2;
   end




%    if(sintheta(ii)>=0.1 && sintheta(ii)<=0.5)
%        val = val + (abs(afdb(ii))-abs(afdbref(ii))).^2;
%    %elseif (sintheta(ii)>=0.5892 && sintheta(ii)<=-0.0226 && afdb(ii)>=msll)
%    elseif (afdb(ii)<msll)
%        if(sintheta(ii)<0.1)
%        val = val + (abs(afdb) - abs(msll)).^2;
%        end
%        
%        if(sintheta(ii)>0.5)
%        val = val + (abs(afdb) - abs(msll)).^2;
%        end
%    end
%    
    
   
end



%val = sum ( (abs(afdb) - abs(afdbref)).^2);
             




