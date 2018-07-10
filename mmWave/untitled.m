
mu = 2;
tstep=0.1;
totaltime =600 ; % must be an integer
R = 100; frac = (1/9);    rhoT = 200*10^(-6);% as per the notation of the paper this is lambda_t
        nTorig = poissrnd(rhoT*pi*R^2) %original AP number (without self-block)
        rT = R*sqrt(rand(nTorig,1)); %location of APs, an array of size nTorigorig*1 is returned 
        alphaT = 2*pi*rand(nTorig,1);%location of APs, an array of size nTorigorig*1 is returned 
        x=(2/pi)*(0.1)*frac*rT; 
            alphai =  1./(x);
 beta = cell(nTorig,1);
for indT = 1:nTorig
    beta{indT} =  exprnd(x(indT),1,totaltime);
end 
dataAP= cell(nTorig,1);
simTime =60;
start_time = unifrnd(0,simTime,nTorig,1);
alpha =  exprnd(1/mu,1,totaltime);
for indT = 1:nTorig
    
    beta{indT};
    
end
for indT = 1:nTorig
        for i = 1:totaltime
            if (i==1)
                timestampBl = beta{indT}(1);
            else
                timestampBl = timestampBl+alpha(i)+ beta{indT}(i);
            end
           if (timestampBl<=simTime)
               dataAP{indT} = [dataAP{indT}, timestampBl];
           end
        end
    
end
for i= 1:nTorig
    dataAP{i};
    length(dataAP{i});
end
for indT = 1:nTorig
    len =length(dataAP{indT});
    dataAP{indT}(2,:) =  exprnd(1/mu,1,len);
end
binary_seq = zeros(nTorig,totaltime);
losses_seq = zeros(nTorig,totaltime); % to store the value of all the losses in a given time sequence
power_received_byvirtueofanyparticular_basestation = zeros(nTorig,totaltime); 
for indT = 1:nTorig
    %     blDur  = exprnd(1/mu);
    for timestamp = 1:length(dataAP{indT})
        blDur  = ceil(dataAP{indT}(2,timestamp)/tstep);
        blTime = ceil(dataAP{indT}(1,timestamp)/tstep);
        if(blTime+blDur<=simTime/tstep)%avoid excess duration
            binary_seq(indT, blTime+1:1:(blTime+blDur))=binary_seq(indT, blTime+1:1:(blTime+blDur))+1;
        end
        
         for particular_inst = 1:1:(totaltime)
            xT = rT(indT).*cos(alphaT(indT));%location of APs
            yT = rT(indT).*sin(alphaT(indT));%location of APs
            if((binary_seq(indT,particular_inst)~=0)&&(blTime+blDur<=simTime/tstep))
                loss = channelModel3(xT,yT,1);
            elseif (blTime+blDur<=simTime/tstep)
                loss = channelModel3(xT,yT,0);
            else 
                continue 
            end
            
            losses_seq(indT,particular_inst) = loss;
            power_received_byvirtueofanyparticular_basestation(indT,particular_inst) = 200-losses_seq(indT,particular_inst); 
        end
    end
end
figure;
hold on
for i = 1:nTorig
     subplot(nTorig,1,i);plot(0.001:tstep:simTime,losses_seq(i,:));
end














