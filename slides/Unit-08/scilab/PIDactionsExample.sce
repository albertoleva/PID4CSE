clear; clc; z=%z;

mu     = 1;
T      = 10;
lambda = 4;
Beta   = 0.6;
Ts     = 1;
tfin   = 50;

K      = T/lambda/mu;
Ti     = 4*T;
Td     = 0.2*Ti;

P      = syslin(Ts,z*Ts/(T+Ts)/(z-T/(T+Ts)));
Cp     = K;
Ci     = syslin(Ts,K*Ts/Ti*z/(z-1));
Cd     = syslin(Ts,K*Td*(1-Beta)/Ts*(z-1)/(z-Beta));
C      = Cp+Ci+Cd;
w2y    = tf2ss(C*P/(1+C*P));
d2y    = tf2ss(P/(1+C*P));
w2u    = tf2ss(C/(1+C*P));
d2u    = tf2ss(-C*P/(1+C*P));
w2up   = tf2ss(Cp/(1+C*P));
w2ui   = tf2ss(Ci/(1+C*P));
w2ud   = tf2ss(Cd/(1+C*P));
d2up   = tf2ss(-P*Cp/(1+C*P));
d2ui   = tf2ss(-P*Ci/(1+C*P));
d2ud   = tf2ss(-P*Cd/(1+C*P));


t      = 0:Ts:tfin;
wd     = ones(t);
wd(1)  = 0;
yw     = dsimul(w2y,wd)';
yd     = dsimul(d2y,wd)';
uw     = dsimul(w2u,wd)';
ud     = dsimul(d2u,wd)';
upidw  = dsimul([w2up;w2ui;w2ud],wd)';
upidd  = dsimul([d2up;d2ui;d2ud],wd)';

hf=scf(0);clf;hf.figure_size = [1200,700];
subplot(321);
 title('Responses to a w unit step'); 
 ylabel('y'); 
 plot(0,0,'k');
 plot(t',ones(t'),'k:');
 plot(t',yw,'b','linewidth',4);
subplot(322);
 title('Responses to a d unit step'); 
 plot(t',zeros(t'),'k:');
 plot(t',yd,'b','linewidth',4);
subplot(323);
 xlabel('time');
 ylabel('u,up,ui,ud (g,b,r,m)');
 plot(0,0,'k');
 plot(t',uw,'g','linewidth',4);
 plot(t',upidw(:,1),'b','linewidth',2);
 plot(t',upidw(:,2),'r','linewidth',2);
 plot(t',upidw(:,3),'m','linewidth',2);
subplot(324);
 xlabel('time');
 plot(0,0,'k');
 plot(t',ud,'g','linewidth',4);
 plot(t',upidd(:,1),'b','linewidth',2);
 plot(t',upidd(:,2),'r','linewidth',2);
 plot(t',upidd(:,3),'m','linewidth',2);
