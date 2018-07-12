% Simulation_Feb20
%Apr3: Generate BS number and location once here and run for various
%blocker density and various omega. Earlier it was di=one in BlockageSimFn

% Mar14: Added code for self blockage (omega) and modified bl/AP densities
% Update Feb20: Generate AP location using PPP

% Update Feb17: Transferred to BlockageSimFn_Feb17.m
% Random Way-point mobility model for blockers
% Simulating Blockage of nT number of APs.
% Generate time sequence of blocked/unblocked periods

close all;
clear;

%----Play-with-values---------------------------------------
aID = getenv('SLURM_ARRAY_TASK_ID');
rng('shuffle');
wannaplot=0;% keep it as zero if you don't want to generate the plots of time sequences
V = 1; %velocity m/s
hb = 1.8;
hr = 1.4;
ht = 5;
frac = (hb-hr)/(ht-hr);
simTime =100; %sec, keep the simulation time at least 100 seconds to get results, otherwise, the code doesn't run 
% simTime = 3*60*60; %sec Total Simulation time
tstep = 0.05; %(sec) time step, to get a better granuality decrease the step size
mu = 2; %Expected bloc dur =1/mu
R = 100; %m Radius
densityBL = 0.005; %(in bl/m^2)
densityAP = 200*10^(-6);% what is the unit of this base station density that is being specified here 
%unit is number_of_BS/m^2
omegaVal = 0;
% densityBL = [0.005,0.01];
% densityAP = [50,100,200,300,400,500]*10^(-6);%(1:1:10)/10^4;
% omegaVal = [0, pi/3];
% for the blocker 
% the following code models the random way point model and returns time
% instances of blockages 
s_input = cell(1,2);
s_mobility = cell(1,2);
for indB=1:length(densityBL)
s_input{indB} = struct('V_POSITION_X_INTERVAL',[-R R],...%(m)
    'V_POSITION_Y_INTERVAL',[-R R],...%(m)
    'V_SPEED_INTERVAL',[10 10],...%(m/s)
    'V_PAUSE_INTERVAL',[0 0],...%pause time (s) % don't want the blocker to intermittently stop 
    'V_WALK_INTERVAL',[1.00 6],...%walk time (s)
    'V_DIRECTION_INTERVAL',[-180 180],...%(degrees)
    'SIMULATION_TIME',simTime,...%(s)
    'NB_NODES',ceil(pi*R^2*densityBL(indB))); % defines the number of blockers in a rectangle of 2r*2r 
s_mobility{indB} = Generate_Mobility(s_input{indB});

for alpha = 1:s_mobility{indB}.NB_NODES

ANGLE{alpha} = atan2(s_mobility{indB}.VS_NODE(alpha).V_SPEED_Y,s_mobility{indB}.VS_NODE(alpha).V_SPEED_X);
 s_mobility{indB}.VS_NODE(alpha).ANGLE = ANGLE{alpha};
end 
end

% disp('finished loop')
% finaldata = zeros(5,length(densityAP),length(densityBL),length(omegaVal)); % making an array of size 5*density of AP, NOT RELEVANT TO THE CAPACITY ANALYSIS  
for indT = 1:length(densityAP) % so as to run it for various values of density of AP and density of BL 
        rhoT = densityAP(indT);% as per the notation of the paper this is lambda_t
        nTorig = poissrnd(rhoT*pi*R^2); %original AP number (without self-block)
        nTorig = 9;
        rT = R*sqrt(rand(nTorig,1)); %location of APs, an array of size nTorig*1 is returned 
        alphaT = 2*pi*rand(nTorig,1);%location of APs, an array of size nTorig*1 is returned 
        alphai=(2/pi)*(0.1)*frac*rT;

    for indB = 1:length(densityBL)
        
        
        for indO = 1:length(omegaVal)
            omega = omegaVal(indO);
            rhoB = densityBL(indB);%0.65;%Rajeev calculated central park
            nB = ceil(pi*R^2*rhoB)%=4000; %number of blockers, since we are concerned about the density of 
            % blockers inside the circle 
            
            AP_input = struct('WANNAPLOT',wannaplot,...
                'RADIUS_AROUND_UE',R,...
                'DENSITY_AP',rhoT,...
                'SIMULATION_TIME',simTime,...
                'TIME_STEP',tstep,...
                'MU',mu,...
                'FRACTION',frac,...
                'SELF_BL_ANGLE_OMEGA',omega,...
                'Original_NUM_AP',nTorig,...
                'LOC_AP_DISTANCE', rT,...
                'LOC_AP_ANGLE',alphaT,...
                'NUM_BL',nB);
            % CHANGES MADE HERE
%             disp('Entering blockage simulation')
              BlockageSimFn_Feb17(s_mobility{indB},AP_input);
%             finaldata(:,indT,indB,indO) = output;
            %         output is [avgFreq,avgDur,probAllBl,th_freqBl,th_durBl,th_probAllBl];
        end
    end
end
% CHANGES MADE HERE 
% csvwrite(strcat('output',num2str(aID),'.csv'),finaldata)

