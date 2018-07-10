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
[arranged_rT,I] = sort(B);
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
length_of_blockers = zeros(nB,1);
for i=1:nB
%     length_of_blockers(i,1)= exprnd(20,1,1);
    length_of_blockers(i,1)= 13.6;
end     
dataAP = cell(nT,2,1);
dataAP_test = cell(nT,1);%contain array of timestamps for all APs no matter which blocker
info_vector= [];
check =1; 
duration = 0;
loop=0;
for indB = 1:nB %for every blocker
%     disp('Running for blocker:')
%     disp(indB)
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
            phi(indB,iter) =  s_mobility.VS_NODE(indB).ANGLE(iter);
          
            theta = alphaT(indT,1);
            angle =(theta-phi(indB,iter));
            
             if ((angle>(-pi/2))&&(angle<(pi/2)))
                 
                    duration = length_of_blockers(indB,1)*(abs(sec(angle)))*0.1;
                    if (duration>=20)
                        duration =20; 
                    end 
             else 
                     duration = 0;
             end 
             if(distance_travelled>=0 && timestampBl<=simTime)
                %                 data{indB,indT} = [data{indB,indT},start_time+blockage_time];
                
                dataAP{indT,1} = [dataAP{indT,1}, timestampBl] ;% to append the new timestamps of blockage arrival in every iteration for every base station
            
                dataAP{indT,2} = [dataAP{indT,2},duration];
               
                dataAP_test{indT} = [dataAP_test{indT},timestampBl];
                check = check+1;
%                data,  dataAP{indT}(2,:) =  0.5*(abs(sec(r)));
             end
        end
      
    end
    
%     for indT = 1:nT
%     dataAP_test{indT} = [dataAP_test{indT},indB];
%     end
end

% disp('finished with blockers')

% time_of_blockers = cell(nB,nT,1);
% check = 0; 
% for indB = 1:nB
%     for indT= 1:nT
%         if((find(dataAP_test{indT}==indB+1)-find(dataAP_test{indT}==indB))==1)
%         time_of_blockers{indB,indT}= 0;
%         end 
%         time_of_blockers{indB,indT}= dataAP_test{indT}((find(dataAP_test{indT}==indB))+1:(find(dataAP_test{indT}==indB+1))-1);
%     end 
% end
% check = zeros(nB,1);
% for indB = 1:nB
%     disp(indB)
%     
%     for indT= 1:nT
%        beta = time_of_blockers{indB,indT}
%        if(beta~=0)
%            check{indB} = check(indB,1)+1;
%        end 
%        disp(indT)
%     end 
% end
% for indB = 1:nB
% check(indB,1)
% end 
% loop_count=0;
% for indT = 1:nT
%     len =length(dataAP{indT})
%     r = 0 + (2*pi).*rand(1,len);
%     dataAP{indT}(2,:) =  0.5*(abs(sec(r)));
% %     dataAP{indT}(2,:) =  0.5 ;
% end
% for indB= 1:nB
% for i = 1:nT
%     if(i==nT)
%         loop_count= loop_count+1;
%         break 
%     end 
%     for j= i+1:nT
%         alpha = time_of_blockers{indB,i};
%         beta = time_of_blockers{indB,j};
%         loop_count= loop_count+1; 
%         if(alpha==0) 
%             break 
%         end 
%         if(beta==0)
%             break 
%         end 
%         tweak_parameter= 5;
%         len1= length(alpha);
%         len2 = length(beta);
%         if(len1>=len2)
%            for m = 1:len1
%                for n = 1:len2
%                    difference = abs(alpha(m)-beta(n));
%                    loop_count= loop_count+1;
%                    if(difference<=tweak_parameter)
%                        check(indB,1) = check(indB,1)+1;
%                         
%                         dataAP{indT}(2,find(alpha(m)==dataAP{i})) = length_of_blockers(indB,1);
% %                         x= dataAP{indT}(2,find(alpha(m)==dataAP{i}))
%                         dataAP{indT}(2,find(beta(n)==dataAP{j})) = length_of_blockers(indB,1);
% %                         y= dataAP{indT}(2,find(beta(n)==dataAP{j}))
%                    end 
%                end 
%            end 
%         end 
%         if(len2>len1)
%            for n = 1:len2
%                for m = 1:len1
%                    difference = abs(alpha(m)-beta(n));
%                    loop_count= loop_count+1;
%                    if(difference<=tweak_parameter)
%                        check(indB,1) = check(indB,1)+1;
%                        dataAP{indT}(2,find(alpha(m)==dataAP{i})) = length_of_blockers(indB,1);
% %                         x= dataAP{indT}(2,find(alpha(m)==dataAP{i}))
%                         dataAP{indT}(2,find(beta(n)==dataAP{j})) = length_of_blockers(indB,1);
% %                         y= dataAP{indT}(2,find(beta(n)==dataAP{j}))
%                    end 
%                end 
%            end 
%         end  
%     end
% end
% end 
% alpha=0;
% for i = 1:nB
%     alpha= check(i,1)+alpha;
% end 
% loop_count;
% alpha;
totaltime = (simTime)/tstep;
binary_seq = zeros(nT,totaltime);% to store whether the particular base station is blocked or not 
losses_seq = zeros(nT,totaltime); % to store the value of all the losses in a given time sequence
power_received_byvirtueofanyparticular_basestation = zeros(nT,totaltime); % to store the power received by virtue of any base station 
% allBl = ones(1,totaltime); %binary seq of all blocked
% Tval = tstep:tstep:totaltime*tstep; %run simulation till tdur with step of tstep
% if(wannaplot),figure; hold on; end
% for indT = 1:nT 
%     alpha = dataAP_test{indT};
%      dataAP_test{indT}((find(dataAP_test{indT}==100))+1:(find(dataAP_test{indT}==101))-1);
% %      disp(indT)
% end 

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
  
        blDur_array  = dataAP{indT,2};
        blTime_array = dataAP{indT,1};
        for timeseq = 1:length(dataAP{indT,1})
            blDur = ceil(blDur_array(1,timeseq)/tstep);
            blTime = ceil(blTime_array(1,timeseq)/tstep);
        if(blTime+blDur<=simTime/tstep&&blDur~=0)%avoid excess duration
            binary_seq(indT, blTime+1:1:(blTime+blDur))=binary_seq(indT, blTime+1:1:(blTime+blDur))+1;
        end
        % MADE CHANGES HERE 
        for particular_inst = 1:1:totaltime
            xT = B(indT).*cos(alphaT(indT,1));%location of APs
            yT = B(indT).*sin(alphaT(indT,1));%location of APs
            if(binary_seq(indT,particular_inst)~=0)
                loss = channelModel3(xT,yT,binary_seq(indT,particular_inst));
            else
                loss = channelModel3(xT,yT,0);
            end
            losses_seq(indT,particular_inst) = loss;
            power_received_byvirtueofanyparticular_basestation(indT,particular_inst) = 84-losses_seq(indT,particular_inst);
%     x= power_received_byvirtueofanyparticular_basestation(indT,particular_inst)
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
end
% to plot the sequences of losses observed for any of the base station 
% figure 
% hold on 
% plot(1:nB,check)
% figure;
% hold on
%  title('Power Received by virtue of any of the Base Station')
figure 
for i = 1:nT
     subplot(nT,1,i)
     plot(tstep:tstep:simTime,power_received_byvirtueofanyparticular_basestation(i,:))
%      plot(tstep:tstep:simTime,power_received_byvirtueofanyparticular_basestation(I(i),:))
end 
figure 
for i= 1:nT
    subplot(nT,1,i)
    plot(tstep:tstep:simTime,binary_seq(i,:))
end 
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
last_power =0;
handover=0;
loop=0;
base_station = zeros(totaltime,1);
last_power_base_station = 1;
a_max_power_base_station = 1;
for indTime = 1:totaltime
 for indT = 1:nT
     
     if(power_received_byvirtueofanyparticular_basestation(indT,indTime)>=max_power)
          max_power = power_received_byvirtueofanyparticular_basestation(indT,indTime);
          
          a_max_power_base_station = indT;
     end
     
     loop= loop+1;
 end
 final_power_received(indTime) =max_power;
%  if(last_power~=max_power)
%          handover= handover+1;
%          final_power_received(indTime)= -80;
%  end 
lower_bound = last_power-5;
upper_bound = last_power+5;
if ((max_power>last_power+2) && (last_power<-15))
    handover= handover+1;
     final_power_received(indTime)= - 100;
      base_station(indTime,1)=  a_max_power_base_station;
end 
% if (last_power~=max_power)
%     handover= handover+1;
%      final_power_received(indTime)= - 80;
% end  
     last_power=max_power;
     max_power= -120;
end
handover
loop
aggregate= 0 ; 
for i = 1:indTime
    aggregate = final_power_received(i,1)+aggregate;
end
aggregate = aggregate/totaltime
capacity = 400*log(1+aggregate/7)
figure 
plot(tstep:tstep:simTime,final_power_received)
title(' Final Power Being Received at the UE Terminal ')
xlabel('Time(sec)')
ylabel('Power(dBm)')
figure 
plot(tstep:tstep:simTime,base_station)
title('Handovers Between the Different Base Stations as time proceeds')
xlabel('Time(sec)')
ylabel('Base Station Number')
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
% nTorig;
% durations = [];
% for i=1:length(dataAP)
%     durations = [durations, dataAP{i}(2,:)];
% end
% to plot ECDF of durations, uncomment next line
% figure; ecdf(durations);
end