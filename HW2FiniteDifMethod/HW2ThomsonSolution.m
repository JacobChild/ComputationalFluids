% ME 541
% Scott Thomson
clear all
% Properties
hcoeff = 1000;
tinf = 300;
% Generate grid using Practice B
NCV = 979; % Number of control volumes
xL = 0; xR = 0.05; % Domain size
dx = (xR - xL)/NCV;
xCV = xL:dx:xR; % Control volume placement
% Place nodes at boundaries & interior
x(1) = 0;
for i = 2:NCV + 1
 x(i) = 0.5*(xCV(i) + xCV(i - 1));
end
x(NCV+2) = xR;
% Compute geometric parameters
dx(1) = 0;
dxw(1) = 0;
dxe(1) = x(2) - x(1);
fw(1) = 0;
fe(1) = 1;
for i = 2:length(x) - 1
 dx(i) = xCV(i) - xCV(i-1);
 dxw(i) = x(i) - x(i-1);
 dxe(i) = x(i+1) - x(i);
 fw(i) = (xCV(i-1) - x(i-1))/(x(i) - x(i-1));
 fe(i) = (x(i+1) - xCV(i))/(x(i+1)-x(i));
end
dx(length(x)) = 0;
dxw(length(x)) = x(length(x)) - x(length(x)-1);
dxe(length(x)) = 0;
fw(length(x)) = 1;
fe(length(x)) = 0;
% Generate k & S values at each node
for i = 1:length(x)
 if x(i)<0.03
 k(i) = 15;
 s(i) = 4e6;
 else
% k(i) = 60;
 k(i) = 137*exp(25*x(i) - 2);
 s(i) = 0;
 end
end
% Initialize coefficient matrix
A = zeros(length(x));
% Fill coefficient matrix for interior nodes
for i = 2:length(x) - 1
 ke = 1/((1 - fe(i))/k(i) + fe(i)/k(i+1))
 kw = 1/((1 - fw(i))/k(i) + fw(i)/k(i-1));

 aw = kw/dxw(i);
 ae = ke/dxe(i);
 ap = ae + aw;
 b(i,1) = s(i)*dx(i);

 A(i,i-1:i+1) = [-aw ap -ae];
end
% Left side boundary condition (Neumann)
disp('first ke ')
ke = 1/((1 - fe(1))/k(1) + fe(1)/k(2))
ae = ke/dxe(1);
ap = ae;
b(1,1) = 0;
A(1,1:2) = [ap -ae];
% Right side boundary condition (mixed)
%disp('last kw ')
kw = 1/((1 - fw(length(x)))/k(length(x)) + fw(length(x))/k(length(x)-1));
aw = kw/dxw(length(x));
ap = aw + hcoeff;
b(length(x),1) = hcoeff*tinf;
A(length(x),length(x)-1:length(x)) = [-aw ap];
%%% Direct solver
T = (A^-1)*b;
% Print & plot results
T
h = plot(x,T,'k');
% Format plot for inserting in document
set(gca,'FontSize',10,'FontName','Times')
grid on
xlabel('{\itx} (m)'); ylabel('{\itTemperature} (K)')
xlim([0 0.05]); ylim([400 600]);
set(gca,'XMinorTick','on','YMinorTick','on','GridLineStyle','--')
set(gcf,'PaperPosition',[1 1 4 3]);
hold all;
print('-dtiff','-r400','figure.tif');
hold off;
