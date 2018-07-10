function loss = channelModel3(xT,yT,considerNLOS)
%Aug 02: Removed the part of code that consider correlation in SF for speed
%July 17: Don't consider correlation in SHADOW FADING
%July 6: Updated 
%July 4: 3GPP UMi channel model
%Input: locBS, locUE: location of BS and UE. considerNLOS
%Output: PL, SF: single PL and SF value for location provided


% if(nargin<1)
%     %     clear all
%     %     close all
%     xSize =64;%100
%     ySize = 24;%100
%     locBS = [xSize/2, -10, 10];
%     locUE = [xSize/4,ySize/4,1.5];%; xSize/4,ySize/4+2,1.5]; %vector
% 
% end % is it relevant to my case ?

xBS = xT; yBS = yT;
xUE = 0; yUE = 0; 
hBS = 5; hUE = 1.4;
fc = 28; % frequency of operation, we assume no doppler effect  
sigma_LOS = 4;
sigma_NLOS = 7.82;

d2D = sqrt((xUE(end)-xBS)^2+(yUE(end)-yBS)^2);
d3D = sqrt(d2D^2+(hBS-hUE)^2);
% implementing the model for path loss based on the 2d and 3d distances ...
% do we want to keep that or bring a change in that since different papers
% site different models for the same 
% Calculate dBPprime, PL1,PL2 as mentioned in 3GPP UMi model
dBPprime = 4*(hBS-1)*(hUE-1)*fc*10/3; %1680m. So, not useful in our scenarios considered 
PL1 = 32.4+21*log10(d3D) +20*log10(fc);

PL2 = 32.4+40*log10(d3D) +20*log10(fc) - 9.5*log10(dBPprime^2+(hBS-hUE)^2);
if(d2D<10), PL_LOS  = 80; %%Err: PL in 3GPP, starts with 10m, what abt d<10m
elseif(d2D<dBPprime), PL_LOS  = PL1;
else, PL_LOS  = PL2;
end

%NLOS PL calculation
PL_pNLOS = 35.3*log10(d3D)+22.4+21.3*log10(fc)-0.3*(hUE-1.5);
if(d2D<10), PL_NLOS  = 85;
else, PL_NLOS  = max(PL_LOS , PL_pNLOS);
end

%Shadowing term
% SF_LOS  =  sigma_LOS*randn(1,1);
% SF_NLOS  = sigma_NLOS*randn(1,1);
 SF_LOS = normrnd(0,sigma_LOS,1,1);
 SF_NLOS  =normrnd(0,sigma_NLOS,1,1); 
%%Calculate probability of LOS channel
%  if(d2D<=18), pLOS  = 1;
%  else, pLOS  = 18/d2D + exp(-d2D/36)*(1-18/d2D);
%  end
%  if (considerNLOS == 0), pLOS = 1; end

%%Use LOS model with prob pLOS, as in 3GPP
% if(rand<pLOS )
%     PL  = PL_LOS ;
%     SF  = SF_LOS ;
% else
%     PL  = PL_NLOS ;
%     SF  = SF_NLOS ;
% end
% CHANGES MADE HERE
 
if(considerNLOS~=0)
    PL  = PL_NLOS ;
    SF  = 35*((rand(1,1))*(considerNLOS));
%     SF = SF_NLOS;
else
    PL = PL_LOS;
    SF = SF_LOS;
end
loss= PL+SF;

% CHANGES MADE HERE 
% % % % These plots does make sense in channelModel.m ;)
% % % figure;plot(d2Drange, pLOS, 'lineWidth',2); xlabel('d2D');ylabel('LOS Probability ')
% % % figure;plot(d2Drange, PL_LOS, 'lineWidth',2); xlabel('d2D');ylabel('LOS PL ')
% % % figure;plot(d2Drange, PL_NLOS, 'lineWidth',2); xlabel('d2D');ylabel('NLOS PL ')
% % % figure;plot(d2Drange, PL+SF, 'lineWidth',2); xlabel('d2D');ylabel('PL + SF ')
% % % figure;plot(d2Drange, corrLOS, 'lineWidth',2); xlabel('d2D');ylabel('Auto-correlation sequence for d2D=100 ')
% % % figure;plot(d2Drange, PL+SFnew, 'lineWidth',2); xlabel('d2D');ylabel('PL + SF (ARIMA)')
 end