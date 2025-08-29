%type 1,2,3 chignolin unbiased
%type 11 trpcage unbiased

type = 3;

if(type == 1)
   filez = sprintf('/mnt/wdred/Folding_paper/chign_charmm22s/unbiased_Nicola/CVS_r2');
    col = 2; 
   fol_lim = 0.15;
    unfol_lim = 0.25;   
    print_stride = 0.01; %ns/frame
elseif(type == 2)
    filez = sprintf('/mnt/wdred/Folding_paper/chign_charmm22s/unbiased_Nicola/CVS_r3');
    col = 2; 
    fol_lim = 0.15;
    unfol_lim = 0.25;    
    print_stride = 0.01; %ns/frame
elseif(type == 3)    
    filez = sprintf('/mnt/wdred/Folding_paper/chign_charmm22s/unbiased_Nicola/CVS_r4');
    col = 2; 
    fol_lim = 0.15;
    unfol_lim = 0.25;
    print_stride = 0.01; %ns/frame
elseif(type == 11)    
    filez = sprintf('/mnt/wdred/Folding_paper/TRPcage/unbiased_K8A_folded/COLVAR_200us');
    col = 2; 
    fol_lim = 0.4;
    unfol_lim = 0.6;
    print_stride = 1.0; %ns/frame
end
matrix = dlmread(filez,'',1,0);

time_runavg = 10.0; %ns time window 
frames_runavg = time_runavg/print_stride;

num_data=size(matrix,1);

clear Q;
Q(1:(num_data),1:8) = 0;
%forward direction
for i = (frames_runavg+1):num_data
    Q(i,1) = i*print_stride; %real time
    Q(i,2) = matrix(i,col);
    Q(i,3) = mean(matrix((i-frames_runavg):i,col)); %running average of Q, forward direction
    if(Q(i,3) < fol_lim) 
        Q(i,4) = 1;  %FOLDED
    end
    if(Q(i,3) > unfol_lim) 
        Q(i,4) = 2;  %UNFOLDED
    end    
end
%backward direction
for i = (num_data-(frames_runavg+1)):(-1):1
    Q(i,1) = i*print_stride; %real time
    Q(i,2) = matrix(i,col);
    Q(i,5) = mean(matrix(i:(i+frames_runavg),col)); %running average of Q, backward direction
    if(Q(i,5) <= fol_lim) 
        Q(i,6) = 1;  %FOLDED
    end
    if(Q(i,5) >= unfol_lim) 
        Q(i,6) = 2;  %UNFOLDED
    end           
end
%initial condition 
%1 for folded
%2 for unfolded
if(Q(1,2) <= fol_lim) 
    Q(i,8) = 1;  %FOLDED
else
    Q(i,8) = 2;  %UNFOLDED    
end

%final check of the two directions
clear timeQ;
timeQ(1:50,1:3) = 0;
counter=1;
for i = 2:num_data
    if(Q(i,4) == Q(i,6)) 
        Q(i,7) = Q(i,4); %only if the forward and backward averages agree, we can assign the folded or unfolded label
    end
         
    state = Q(i,7); %I have 3 kinds of states 1 folded, 2 unfolded, 0 none  
    
    %if(state == Q(i-1,7)) %accumulate time, no change in the state at all, including state0 TS        
    if((state == 0) || (state == Q(i-1,8))) %accumulate time, no change in the state or just state 0 (TS)
        timeQ(counter,3) = timeQ(counter,3) + print_stride; %transition time
        Q(i,8) = Q(i-1,8); %macro state does not change, remains 1 or 2 as it was
    else     
        Q(i,8) = Q(i,7); %macro state changes
        timeQ(counter,2) = Q(i-1,8); %label of event just finished
        timeQ(counter,3) = timeQ(counter,3)/1000; %convert the time of last event in us
        counter = counter+1; %new event
        timeQ(counter,1) = i;%*print_stride; %timestamp of start of new event, in units of row
    end
end
%last point of the list
timeQ(counter,2) = Q(i-1,8); %label of event just finished
timeQ(counter,3) = timeQ(counter,3)/1000; %convert the time of last event in us

%split data
clear timeF;
clear timeU;
timeF(1:10,1:3) = 0;
timeU(1:10,1:3) = 0;
counterF=1;
counterU=1;
for j=1:size(timeQ,1)
    if(timeQ(j,2) == 1)    
        timeF(counterF,:) = timeQ(j,:);
        counterF = counterF+1;
    elseif(timeQ(j,2) == 2)    
        timeU(counterU,:) = timeQ(j,:);
        counterU = counterU+1;        
    end
end

av_lifet_F = mean(timeF(:,3))
av_lifet_U = mean(timeU(:,3))
num_events = counter

if(type == 1)
    timeF1 = timeF;
    timeU1 = timeU;
elseif(type == 2)    
    timeF2 = timeF;
    timeU2 = timeU;
elseif(type == 3)    
    timeF3 = timeF;
    timeU3 = timeU;
    timeFtot = [timeF1; timeF2; timeF3];
    timeUtot = [timeU1; timeU2; timeU3];
    av_lifet_F = mean(timeFtot(:,3))
    av_lifet_U = mean(timeUtot(:,3))
    num_events = size(timeUtot,1)
end
