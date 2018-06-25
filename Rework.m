mu=2;
R = 100; %m Radius
densityBL = 0.05; %(in bl/m^2)
densityAP = 200*10^(-6);% what is the unit of this base station density that is being specified here 
%unit is number_of_BS/m^2
frac = (1/9);    rhoT = 200*10^(-6);% as per the notation of the paper this is lambda_t
        nTorig = poissrnd(rhoT*pi*R^2); %original AP number (without self-block)
        rT = R*sqrt(rand(nTorig,1)); %location of APs, an array of size nTorigorig*1 is returned 
        alphaT = 2*pi*rand(nTorig,1);
V = 1;
simTime =100; %sec% keep this as this so as to run the simulations for the generate mobility to work 
nB = 4*R^2*densityBL;
tstep = 0.1;
% simTime = 3*60*60; %sec Total Simulation time
% nT = poissrnd(densityAP*pi*R^2);
s_input = struct('V_POSITION_X_INTERVAL',[-R R],...%(m)
    'V_POSITION_Y_INTERVAL',[-R R],...%(m)
    'V_SPEED_INTERVAL',[V V],...%(m/s)
    'V_PAUSE_INTERVAL',[0 0],...%pause time (s) % don't want the blocker to intermittently stop 
    'V_WALK_INTERVAL',[1.00 60.00],...%walk time (s)
    'V_DIRECTION_INTERVAL',[-180 180],...%(degrees)
    'SIMULATION_TIME',simTime,...%(s)
    'NB_NODES',4*R^2*densityBL ); % defines the number of blockers in a rectangle of 2r*2r 
s_mobility = Generate_Mobility(s_input);

dataAP = cell(nTorig,1); %contain array of timestamps for all APs no matter which blocker
xT = rT.*cos(alphaT);%location of APs
yT = rT.*sin(alphaT);%location of APs
xTfrac = frac*xT; %blockage zone around UE for each APs
yTfrac = frac*yT;
locT = [xTfrac';yTfrac']; %2 rows for x and y, nT columns

for indB = 1:nB %for every blocker
    
    for iter =1:(length(s_mobility.VS_NODE(indB).V_POSITION_X)-1)
        
        % for every time blocker changes direction
        loc0 = [s_mobility.VS_NODE(indB).V_POSITION_X(iter);...
            s_mobility.VS_NODE(indB).V_POSITION_Y(iter)];
        loc1 = [s_mobility.VS_NODE(indB).V_POSITION_X(iter+1);...
            s_mobility.VS_NODE(indB).V_POSITION_Y(iter+1)];
        start_time = s_mobility.VS_NODE(indB).V_TIME(iter);
        velocity = sqrt((s_mobility.VS_NODE(indB).V_SPEED_X(iter))^2+ ...
            (s_mobility.VS_NODE(indB).V_SPEED_Y(iter))^2);
        for indT = 1:nTorig %for every AP
            distance_travelled = find_blockage_distance([loc0,loc1],locT(:,indT),alphaT(indT));
            timeToBl = distance_travelled/velocity; %time to blocking event
            timestampBl = start_time+timeToBl; %timestamp of blockage event
            if(distance_travelled>=0 && timestampBl<=simTime)
                %                 data{indB,indT} = [data{indB,indT},start_time+blockage_time];
                dataAP{indT} = [dataAP{indT}, timestampBl];
                
            end
        end
    end
    
end


totaltime = (simTime)/tstep;
binary_seq = zeros(nTorig,totaltime);
allBl = ones(1,totaltime); %binary seq of all blocked
losses_seq = zeros(nTorig,totaltime); % to store the value of all the losses in a given time sequence
power_received_byvirtueofanyparticular_basestation = zeros(nTorig,totaltime); 

for indT = 1:nTorig
    len =length(dataAP{indT});
    dataAP{indT}(2,:) =  exprnd(1/mu,1,len);
end

for indT = 1:nTorig
    %     blDur  = exprnd(1/mu);
    for timestamp = 1:size(dataAP{indT},2)
        blDur  = ceil(dataAP{indT}(2,timestamp)/tstep);
        blTime = ceil(dataAP{indT}(1,timestamp)/tstep);
        if(blTime+blDur<=simTime/tstep)%avoid excess duration
            binary_seq(indT, blTime+1:1:(blTime+blDur))=binary_seq(indT, blTime+1:1:(blTime+blDur))+1;
        end
    end
    for particular_inst = 1:1:totaltime
            xT = rT(indT).*cos(alphaT(indT));%location of APs
            yT = rT(indT).*sin(alphaT(indT));%location of APs
            if(binary_seq(indT,particular_inst)~=0)
                loss = channelModel3(xT,yT,1);
            else
                loss = channelModel3(xT,yT,0);
            end
            losses_seq(indT,particular_inst) = loss;
            power_received_byvirtueofanyparticular_basestation(indT,particular_inst) = 200-losses_seq(indT,particular_inst); 
    end
end

figure;
hold on
for i = 1:nTorig
     subplot(nTorig,1,i);plot(tstep:tstep:simTime,losses_seq(i,:));
end