function out= myintd(theta, sig1, sig2,zplus,n,K1D,t0)





ck=1-((theta+zplus)/(2*t0));
ka=acos(ck);
sig1(1,1)=-t0*exp(i*ka);
gam1=i*(sig1-sig1');

ck=1-((theta+zplus)/(2*t0));
ka=acos(ck);
sig2(n,n)=-t0*exp(i*ka);
gam2=i*(sig2-sig2');

G=inv(((theta+zplus)*eye(n))-K1D-sig1-sig2);


out =G;  




end