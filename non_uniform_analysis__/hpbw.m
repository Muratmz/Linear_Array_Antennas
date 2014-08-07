% HPBWCALC
function[hp,thmax]=hpbw(U,M)
tol=0.001;
imax=0;
j=0;
for i=1:M+1;
   if abs(U(i)-1)<tol & floor((j+2)/2)==imax+1,
      imax=imax+1;
      thmax(imax)=(i-1)/10;
   end
   if i>1 & abs(U(i)-1)<tol & U(i)>U(i-1) & j~=0,
      thmax(imax)=(i-1)/10;
   end
   if i>1,
      y(1)=U(i)-0.5;
      y(2)=U(i-1)-0.5;
      x(1)=(i-1)/10;
      x(2)=(i-2)/10;
      sign=y(1)*y(2);
      if sign<0,
         j=j+1;
         root(j)=x(2)-y(2)*(x(2)-x(1))/(y(2)-y(1));
         if j>=2 & y(2)>y(1),
            hp(imax)=root(j)-root(j-1);
         elseif j==1 & y(2)>y(1),
            hp(imax)=2.*root(j);
         end
      end
   end
end
if thmax(imax)>root(j),
   hp(imax)=2.*(180-root(j));
end
