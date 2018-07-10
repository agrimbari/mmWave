figure 
hold on 

subplot(2,2,1) 
d = linspace(1,100);
SF = normrnd(0,4,1,100);
PL1 = 32.4+21*log10(d) +20*log10(28);
plot(d,PL1+SF(1,:))
title('Pathloss and Shadowing for LOS ')
xlabel('Distance(m)')
ylabel('Power(dBm)')
subplot(2,2,2)
plot(d,PL1)
title('Path loss for LOS')
xlabel('Distance(m)')
ylabel('Power(dBm)')
subplot(2,2,3)
SF_2 = 35*((rand(1,100)));
PL_pNLOS = 35.3*log10(d)+22.4+21.3*log10(28)-0.3*(1.4-1.5);
plot(d,PL_pNLOS+SF_2(1,:))
title('Pathloss and Shadowing for NLOS ')
xlabel('Distance(m)')
ylabel('Power(dBm)')
subplot(2,2,4)
plot(d,PL_pNLOS)
title('Pathloss for NLOS ')
xlabel('Distance(m)')
ylabel('Power(dBm)')