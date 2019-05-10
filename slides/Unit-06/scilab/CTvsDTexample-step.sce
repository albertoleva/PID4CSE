clear; clc;

T   = 1;
mu  = 1;
Ts  = [0.1,0.2,0.5,1];

t   = 0:T/10:8*T;
yCT = mu*(1-exp(-t/T));

hf  =scf(0);
hf.figure_size = [850,500]
clf;
prows = ceil(sqrt(length(Ts)));
pcols = prows;

for i=1:length(Ts)
    rho = mu*Ts(i)/(T+Ts(i));
    p   = T/(T+Ts(i));
    P   = syslin('d',%z*rho/(%z-p));
    k   = 0:Ts(i):8*T;
    y   = dsimul(tf2ss(P),ones(k));
    subplot(pcols,prows,i)
    plot(t,yCT,'k:');
    plot2d2(k,y,2);
    //plot(k,y,'b.')
    title(sprintf("Ts = %.1f",Ts(i)));
end



