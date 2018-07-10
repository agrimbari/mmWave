function BlockageSimFn_Feb17(s_mobility,AP_input)
%Update Apr13: Move s_mobility to Simulation.m
% Update Feb 28: Save whole DataAP for all AP and BL density for all aID
%               : Then write separate code for analysis and theoretical
%               plots
% BlockageSim_Feb17
% Random Way-point mobility model for blockers
% Simulating Blockage of nT number of APs.
% Generate time sequence of blocked/unblocked periods
% AND operator on all those nT time sequences.
% FInally, report the blockage freq and duration of all nT APs
%%frequency = (per sec)count transitions from 0 to 1 devided by simTime
%%duration: (sec) count total # of 1's multiplied by tstep/blockageCount


%----Play-with-values-here--------------------------------------
wannaplot = AP_input.WANNAPLOT; %1;

nB = AP_input.NUM_BL;%4*R^2*rho_b;%=4000; %number of blokers

nTorig = AP_input.Original_NUM_AP
rT =AP_input.LOC_AP_DISTANCE; %location of APs
alphaTorig = AP_input.LOC_AP_ANGLE;%location of APs
% alphai = AP_input.DENSITY_BLCK;
frac = AP_input.FRACTION;
omega = AP_input.SELF_BL_ANGLE_OMEGA;  

tempInd =  find(alphaTorig>=omega); % for evaluating the number of blockers that are outside the self blockage range 
% xT = rT(tempInd).*cos(alphaTorig(tempInd));%location of APs
% yT = rT(tempInd).*sin(alphaTorig(tempInd));%location of APs
nT = length(tempInd)
% nT=0
if(nT==0) 
    b= 100  
    return;
end % Dealing zero APs
alphaT = zeros(length(tempInd),1);
B = zeros(length(tempInd),1);
for i = 1:length(tempInd)
alphaT(i) = alphaTorig(tempInd(i));
B(i,1)= rT(tempInd(i));
end
[arranged_rT,I] = sort(B)
xT = B.*cos(alphaT);%location of APs
yT = B.*sin(alphaT);%location of APs
xTfrac = frac*xT; %blockage zone around UE for each APs
yTfrac = frac*yT;
locT = [xTfrac';yTfrac']; %2 rows for x and y, nT columns
simTime = AP_input.SIMULATION_TIME; %sec Total Simulation time
tstep = AP_input.TIME_STEP; %(sec) time step
mu = AP_input.MU; %Expected bloc dur =1/mu

%---------------I am moving this to Sim...m--------
% locT = AP_input.T_EFF_LOCATION; %AP location
% alphaT = AP_input.T_ANGLE;

% s_mobility = Generate_Mobility(s_input);  

%------------------------------------------------------

dataAP = cell(nT,1); %contain array of timestamps for all APs no matter which blocker
info_vector= [];
check =1; 
for indB = 1:nB %for every blocker
    disp('Running for blocker:')
    disp(indB)
    for iter =1:(length(s_mobility.VS_NODE(indB).V_POSITION_X)-1)
        
        % for every time blocker changes direction
        loc0 = [s_mobility.VS_NODE(indB).V_POSITION_X(iter);...
            s_mobility.VS_NODE(indB).V_POSITION_Y(iter)];
        loc1 = [s_mobility.VS_NODE(indB).V_POSITION_X(iter+1);...
            s_mobility.VS_NODE(indB).V_POSITION_Y(iter+1)];
        start_time = s_mobility.VS_NODE(indB).V_TIME(iter);
        velocity = sqrt((s_mobility.VS_NODE(indB).V_SPEED_X(iter))^2+ ...
            (s_mobility.VS_NODE(indB).V_SPEED_Y(iter))^2);
        for indT = 1:nT %for every AP
            distance_travelled = find_blockage_distance([loc0,loc1],locT(:,indT),alphaT(indT));
            timeToBl = distance_travelled/velocity; %time to blocking event
            timestampBl = start_time+timeToBl; %timestamp of blockage event
            if(distance_travelled>=0 && timestampBl<=simTime)
                %                 data{indB,indT} = [data{indB,indT},start_time+blockage_time];
               dataAP{indT} = [dataAP{indT}, timestampBl];% to append the new timestamps of blockage arrival in every iteration for every base station
               
            end
        end
        
    end
  info_vector= [info_vector,indB];
  for i=1:length(dataAP)
    info_vector = [info_vector, dataAP{i}];
  end
end
x = info_vector(find(1):find(2))

time_blocker= cell(nB,nT,1) ;
for i=1:nB
for j= 1:nT
    len = length(dataAP{j});
    temp_var =find(i);
%     time_blocker{i,j}
end 
end 
disp('finished with blockers')

totaltime = (simTime)/tstep;
binary_seq = zeros(nT,totaltime);% to store whether the particular base station is blocked or not 
losses_seq = zeros(nT,totaltime); % to store the value of all the losses in a given time sequence
power_received_byvirtueofanyparticular_basestation = zeros(nT,totaltime); % to store the power received by virtue of any base station 
% allBl = ones(1,totaltime); %binary seq of all blocked
% Tval = tstep:tstep:totaltime*tstep; %run simulation till tdur with step of tstep
% if(wannaplot),figure; hold on; end
for indT = 1:nT 
    alpha = dataAP{indT};
end 
for indT = 1:nT
    len =length(dataAP{indT});
    r = 0 + (2*pi).*rand(1,len);
    dataAP{indT}(2,:) =  0.5*abs(sec((r));
end
 % MADE CHANGES HERE 
% indT = plot_input.indT;
% indB = plot_input.indB;
% aID = plot_input.aID;
% save(strcat('dataAP_',num2str(aID),...
%     '_',num2str(indB),...
%     '_',num2str(indT),'.mat'),'dataAP')
% csvwrite(strcat('output',num2str(aID),'.csv'),finaldata)
% arrangement done for applying the association rule based on distance 
for indT = 1:nT
    disp('Running for AP:')
    disp(indT)
    %     blDur  = exprnd(1/mu);
    for timestamp = 1:size(dataAP{indT},2)
        blDur  = ceil(dataAP{indT}(2,timestamp)/tstep);
        blTime = ceil(dataAP{indT}(1,timestamp)/tstep);
        if(blTime+blDur<=simTime/tstep)%avoid excess duration
            binary_seq(indT, blTime+1:1:(blTime+blDur))=binary_seq(indT, blTime+1:1:(blTime+blDur))+1;
        end
        % MADE CHANGES HERE 
        for particular_inst = 1:1:totaltime
            xT = B(indT).*cos(alphaT(indT,1));%location of APs
            yT = B(indT).*sin(alphaT(indT,1));%location of APs
            if(binary_seq(indT,particular_inst)~=0)
                loss = channelModel3(xT,yT,1);
            else
                loss = channelModel3(xT,yT,0);
            end
            losses_seq(indT,particular_inst) = loss;
            power_received_byvirtueofanyparticular_basestation(indT,particular_inst) = 35-losses_seq(indT,particular_inst); 
        end
        % MADE CHANGES HERE 
    % call the function that tells you about the path loss and shadow
    % fading for a given base station, and now based on the association
    % rule you will select a different base station that is not inside the
    % self blockage and for that you will compute the path loss and self
    % blockage and after step size you will again check through the use of
    % association rule which is the nearest base station that is available
    % for you to use and in the case that all the base stations are blocked
    % then for the nearest base station compute the nlos path loss and keep
    % that as the path and shadow fading value provided we dont assume that
    % the base station cant multiplex and thus increase the capacity
    % thereby 
    % general purpose make an array which can for every base station based
    % on whether the base station is active or not stores the value path
    % loss plus shadow fading depending on whether you are making a
    % consideration for nlos or los and this is done for every base
    % station that is outside the self blockage range 
%     allBl = allBl & binary_seq(indT,:);
%     if(wannaplot)
%         subplot(nT+1,1,indT)
%         plot(binary_seq(indT,1:10/tstep), 'lineWidth',4)
%         set(gca,'xtick',[])
%         set(gca,'xticklabel',[])
%         set(gca,'ytick',[])
%         set(gca,'yticklabel',[])
%     end
    end
end
% to plot the sequences of losses observed for any of the base station 
figure;
hold on
for i = 1:nT
     subplot(nT,1,i)
     plot(tstep:tstep:simTime,power_received_byvirtueofanyparticular_basestation(I(i),:))
end   
% for i= 1:nT
%     subplot(nT,2,i+nT)
%     plot(tstep:tstep:simTime,losses_seq(i,:))
% end 
% losses_seq(i,:)
% power_received_byvirtueofanyparticular_basestation(indT,particular_inst)
final_power_received = zeros(totaltime,1);
% CHANGES MADE HERE % APPLYING THE ASSOCIATION RULE 
% keep it set to the value 0 if for any of the given time stamp the power doesn't go down below this min_value
a=1;% to work here when you get back !! 
% Association rule 1
% for indTime = 1:totaltime 
%      while binary_seq(I(a),indTime)
%          if(a==length(I))
%              break
%          end
%          a=a+1;   
%      end
%  final_power_received(indTime) = power_received_byvirtueofanyparticular_basestation(I(a),indTime);
%  
% %  csvwrite(strcat('final_power_received',num2str(indTime),'.csv'),final_power_received);
%   a=1;
% end
% association rule 2 
max_power=-120;
for indTime = 1:totaltime
 for indT = 1:nT
     if(power_received_byvirtueofanyparticular_basestation(indT,indTime)>=max_power)
          max_power = power_received_byvirtueofanyparticular_basestation(indT,indTime);
     end
 end 
 final_power_received(indTime) =max_power;
     max_power= -120;
end
aggregate= 0 ; 
for i = 1:indTime
    aggregate = final_power_received(i,1)+aggregate;
end
aggregate = aggregate/totaltime
figure 
hold on 
plot(tstep:tstep:simTime,final_power_received)
%CHANGES MADE HERE

% if(wannaplot)
%     subplot(nT+1,1,nT+1)
%     plot(Tval(1:10/tstep),allBl(1:10/tstep),'r-', 'lineWidth',4)
%     xlabel('Time (sec)')
% end
% plot(binary_seq(1,:), 'lineWidth',4)

%%Evaluate frequency and average duration of blockage
% CHANGES MADE HERE
% avgFreq = sum(diff(allBl)>0)/simTime;
% avgDur = sum(allBl)*tstep/sum(diff(allBl)>0);
% probAllBl = sum(allBl)*tstep/simTime;
% CHANGES MADE HERE 
% ? HOW TO STORE AN ARRAY or pass it to the main code so that it can be
% stored in a csv file 
% output=[avgFreq,avgDur,probAllBl,nTorig,nT];
nTorig
durations = [];
for i=1:length(dataAP)
    durations = [durations, dataAP{i}(2,:)];
end
% to plot ECDF of durations, uncomment next line
% figure; ecdf(durations);
end