% TSCEBY(THETA,NELEM,D,RDB)
function[AF,Ncoef,Coef]=tscheby(theta,Nelem,d,RdB);
Ro=10^(RdB/20);
P=Nelem-1;
Zo=0.5*((Ro+sqrt(Ro^2-1))^(1/P)+(Ro-sqrt(Ro^2-1))^(1/P));
if 2*floor(Nelem/2)==Nelem,
   Ncoef=Nelem/2;
   M=Ncoef;
for i=1:M;
   Coef(i)=0;
   for j=i:M;
      Coef(i)=Coef(i)+(-1)^(M-j)*Zo^(2*j-1)*fact(j+M-2)*(2*M-1)/(fact(j-i)*fact(j+i-1)*fact(M-j));
   end
end
elseif 2*floor((Nelem+1)/2)==Nelem+1,
   Ncoef=(Nelem+1)/2;
   M=Ncoef-1;
for i=1:M+1;
   Coef(i)=0;
   for j=i:M+1;
      if i==1,EN=2;
      else EN=1;
      end
      Coef(i)=Coef(i)+(-1)^(M-j+1)*Zo^(2*(j-1))*fact(j+M-2)*2*M/(EN*fact(j-i)*fact(j+i-2)*fact(M-j+1));
   end
end
end
u=pi*d*cos(theta);
if 2*floor(Nelem/2)==Nelem,
   AF=0;
   for i=1:Ncoef;
      AF=AF+Coef(i)*cos((2*i-1)*u);
   end
elseif 2*floor((Nelem+1)/2)==Nelem+1, 
   AF=0;
   for i=1:Ncoef;
      AF=AF+Coef(i)*cos(2*(i-1)*u);
   end
end