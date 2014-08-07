% spherical plot
function spherical_plot(r,THETA,PHI,disc)
%theta = linspace(theta_low,theta_up,disc);
%phi   = linspace(phi_low,phi_up,disc);

%[THETA,PHI] = meshgrid(theta,phi);

% spherical to rectangular conversion
x = abs(r).*sin(THETA).*cos(PHI);
y = abs(r).*sin(THETA).*sin(PHI);
z = abs(r).*cos(THETA);

% do the plot
figure; surf(x,y,z); view(135,20);
C = [.8 .8 .8]; colormap(C); axis off equal;

% Draw x, y, and z axes
set(line([1e-8;max(max(x))+3],[1e-8;1e-8],[1e-8;1e-8]),'Color','r');
set(line([1e-8;1e-8],[1e-8;max(max(y))+3],[1e-8;1e-8]),'Color','r');
set(line([1e-8;1e-8],[1e-8;1e-8],[1e-8;max(max(z))+3]),'Color','r');

% Label x, y, and z axes
text(max(max(x))+4,0,0,'x','FontSize',14,'FontName','Times','FontAngle','italic','Color','r');
text(0,max(max(y))+4,0,'y','FontSize',14,'FontName','Times','FontAngle','italic','Color','r');
text(0,0,max(max(z))+4,'z','FontSize',14,'FontName','Times','FontAngle','italic','Color','r');

% Fill surface using patches
patch_1 = zeros(3,disc+1);  patch_2 = zeros(3,disc+1);
patch_1(1,1:disc) = x(1,:); patch_2(1,1:disc) = x(disc,:);
patch_1(2,1:disc) = y(1,:); patch_2(2,1:disc) = y(disc,:);
patch_1(3,1:disc) = z(1,:); patch_2(3,1:disc) = z(disc,:);
patch(patch_1(1,:),patch_1(2,:),patch_1(3,:),C);
patch(patch_2(1,:),patch_2(2,:),patch_2(3,:),C);
