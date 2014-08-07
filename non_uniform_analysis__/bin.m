% BINOMIAL
function[AF,Ncoef,Coef]=bin(theta,Nelem,d)
%Check if nelement is even or odd and make ncoef = .5 * nelemen
if 2*floor(Nelem/2)==Nelem,Ncoef=Nelem/2;
else Ncoef=(Nelem+1)/2;
end

%create the pascal's triangle
for i=1:Ncoef;
   Coef(i)=1;
   for j=1:Ncoef-i;
      Coef(i)=Coef(i).*(Nelem-j)./j;
   end
end

%pascal triangle entries stored in Coef
if 2*floor(Nelem/2)~=Nelem,Coef(1)=Coef(1)/2;
end
u=pi*d*cos(theta);
%implement array factor equation 6-66 pg 294 balanis
if 2*floor(Nelem/2)==Nelem,
   AF=0;
   for i=1:Ncoef;
      AF=AF+Coef(i).*cos((2.*i-1).*u);
   end
else AF=0;
   for i=1:Ncoef;
      AF=AF+Coef(i).*cos(2.*(i-1).*u);
   end
end
