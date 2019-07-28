%%%equilibrium charge for f1 =f2=1;
clear all;
close all;

n = 4;
h = 1/(n+1);
x = h:h:1-h;
Ef = 4;
t0 = 1;
zplus=i*1e-9;

K1D = spdiags(ones(n,1)*[-1 2 -1],-1:1,n,n);  % 1d Poisson matrix




Emin = -6;
Emax = Ef;
eecentricity = 0;
R = (Emax-Emin)/2;
R2= (Emax+Emin)/2;
g1 = @(theta) (R*cos(theta)+1i*(1-eecentricity)*R*sin(theta))+R2;
g1prime = @(theta) -R*sin(theta)+1i*(1-eecentricity)*R*cos(theta); 



sig1=zeros(n);sig2=zeros(n);



myfun=@(theta) myintd(theta, sig1, sig2,zplus,n,K1D,t0); %integrand

myQ1= integral(@(t) myfun(g1(t))*g1prime(t),0,pi,'AbsTol',1e-9,'ArrayValued', true);
% myQ2= quadv(@(t) myfun(g1(t))*g1prime(t),0,pi,'AbsTol',1e-9);
%quadv will have some problem

gn1 = (-1/pi)*1i*diag(myQ1); %real part is electron density
% gn2 = (-1/pi)*1i*diag(myQ2); %real part is electron density


