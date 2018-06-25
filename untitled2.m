
figure 
hold on 
subplot(3,1,1);
d = linspace(1,100);
SF = normrnd(0,4,1,100);

% for i= 1:100 
% PL1 = 32.4+21*log10(i) +20*log10(28)+SF(i,1);
% end

PL1 = 32.4+21*log10(d) +20*log10(28)
plot(d,PL1+SF(1,:));

subplot(3,1,2)
d = linspace(1,100);
SF_2 = normrnd(0,7.82,1,100);
PL_pNLOS = 35.3*log10(d)+22.4+21.3*log10(28)-0.3*(1.4-1.5);

plot(d,PL_pNLOS+SF_2(1,:));

