%%
%This code demonstrates the implementation of schelkenoff's, which is based 
%on the information about the null possition on the final radiation pattern.

% Program written by: Sathvik N. Prasad
% Date              : 23/06/2014 


close all;
clc;
clear all;


%%
fprintf('Determination of linear array antenna directivity by given null points location \n');
d=input('Specify spacing between elements (in wavelengths)\n');
beta=0;
beta=beta*pi/180;
Nnul=input('Specify number of desired nulls\n');
Nel=Nnul+1;
nul=zeros(size(Nnul));

	for n=1:Nnul,    
     	    nul(n)=input(['Specify position of null #',num2str(n), ...
                   ' (in degrees)\n']);
    end;

   nul=nul*pi/180;   % convert to rads


   Ntheta=360; 
   thres=-120;
   
   %used for plotting
   
   
theta=linspace(0,2*pi,2*Ntheta);
dth=theta(2)-theta(1);   % theta in radians

psi=2*pi*d*cos(theta)+beta;

	     psi_nul=2*pi*d*cos(nul)+beta;
	     a=poly(exp(j*psi_nul))  % coefficients of z polynomial
        
    
        
	     z=exp(j*psi);

	     for n=1:Nel,
	        temp(n,:)=z.^(n-1);
	     end;

	     AF=fliplr(a)*temp;  % final array factor
        	
        for mu=1:ceil(d*(1-min(cos(nul)))),    % extra code to cover the case
					       % d>0.5 lambda  
           theta_nul1(mu,:)=acos(cos(nul)+mu/d)*180/pi;
           theta_nul2(mu,:)=acos(cos(nul)-mu/d)*180/pi;
        end;

	     theta_nul=[theta_nul1(:) theta_nul2(:)];
        theta_nul(abs(imag(theta_nul))>=1e-4)=[];
	     nul_tot=[nul*180/pi theta_nul];

	     for w=1:length(nul_tot),   %  remove multiple null
                                    %  occurences
	         nul_final(w)=nul_tot(1);
            nul_tot(abs(nul_tot-nul_final(w))<=1e-2)=[];
	         if (isempty(nul_tot)),
                break;
	         end;
	     end;
        
% 
% 	     disp('Array Factor nulls found at:');
% 		  disp(strcat(num2str(sort(nul_final)','%6.2f'), ...
%              repmat(' deg.',length(nul_final),1)))                  
%         fprintf(2,'\nExcitation coefficients:\n');  
%         fprintf(2,'   Real part \t    Imag part\n');
%         fprintf(2,'%10.4f \t %10.4f\n',[real(fliplr(a)); imag(fliplr(a))]);
  
          
       
    

	AF_db=20*log10(abs(AF)/max(abs(AF)));
	AF_db(AF_db<=thres)=thres;
   
%% directivity computation

U=(abs(AF(theta<pi))).^2;
Prad=2*pi*sum(U.*sin(theta(theta<pi))*dth);
D=4*pi*max(U)/Prad;
D_db=10*log10(D);
   



%% Figure 
figure(1);
set(gca,'fontname','timesnewroman','fontsize',14); 
hl=plot(theta*180/pi,AF_db,'LineWidth',3,'color','black'); grid; zoom; thres=get(gca,'ylim'); grid on;


%%Description regarding the directivity plot
set(gca,'fontname','timesnewroman','fontsize',14,'xlim',[0 180]);
set(hl,'linewidth',2); 
xlabel('\theta (in degrees)','fontname','timesnewroman','fontsize',18, ...
       'verticalalignment','top');
ylabel('Normalized Array Factor (in dB)','fontname','timesnewroman', ...
       'fontsize',18,'verticalalignment','bottom');
title(['Synthesized Array Factor using Schelkunoff polynomial(N = ',num2str(Nel),', d = ',num2str(d),'\lambda)'], ...
      'fontname','timesnewroman','fontsize',18);


