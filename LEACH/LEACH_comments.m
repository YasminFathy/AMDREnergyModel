%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                      %
% SEP: A Stable Election Protocol for clustered                        %
%      heterogeneous wireless sensor networks                          %
%                                                                      %
% (c) Georgios Smaragdakis                                             %
% WING group, Computer Science Department, Boston University           %
%                                                                      %
% You can find full documentation and related information at:          %
% http://csr.bu.edu/sep                                                %
%                                                                      %  
% To report your comment or any bug please send e-mail to:             %
% gsmaragd@cs.bu.edu                                                   %
%                                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                      %
% This is the LEACH [1] code we have used.                             %
% The same code can be used for FAIR if m=1                            %
%                                                                      %
% [1] W.R.Heinzelman, A.P.Chandrakasan and H.Balakrishnan,             %
%     "An application-specific protocol architecture for wireless      % 
%      microsensor networks"                                           % 
%     IEEE Transactions on Wireless Communications, 1(4):660-670,2002  %
%                                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%
clear;
% Add: Yasmin
clc
close all
% end: Yasmin
rng('default'); % to repeat a result that you got 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add: Yasmin
% if countCHs keeps = 0 or it's equal to the number of nodes which
% means that 100% of the nodes are cluster heads. In both cases, this is
% considered as Direct Transmission.
% This in noted from the paper, entitled
% Energy-Efficient Communication Protocol forWireless Microsensor Networks
% end: Yasmin

%Field Dimensions - x and y maximum (in meters)
xm=100;
ym=100;

%x and y Coordinates of the Sink
sink.x=0.5*xm;
sink.y=0.5*ym;

% Yasmin
%sink.x=50;
%sink.y=175;

%Number of Nodes in the field
n=52

% Add: Yasmin
%if p=1 it means we have fair distributed head
% end: Yasmin


%Optimal Election Probability of a node
%to become cluster head
p=0.1;

%Energy Model (all values in Joules)
%Initial Energy 
Eo=0.5;
%Eelec=Etx=Erx
ETX=50*0.000000001;   %% Yasmin: transmission cost 
ERX=50*0.000000001;   %% Yasmin: receiving cost

% Add: Yasmin
% this is the electronic energy, the values of  Efs and Emp
% depend on distance to the receiver and the acceptable bit-error rate.
% end: Yasmin

%Transmit Amplifier types
Efs=10*0.000000000001;
Emp=0.0013*0.000000000001;
%Data Aggregation Energy
EDA=5*0.000000001;


%Values for Hetereogeneity
%Percentage of nodes than are advanced
m=0; 



%\alpha
a=0;

%maximum number of rounds
rmax=4000; %9999 is default

%%%%%%%%%%%%%%%%%%%%%%%%% END OF PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%

%Computation of do
do=sqrt(Efs/Emp);

%Creation of the random Sensor Network
figure(1);
for i=1:1:n
    % Add: Yasmin
    % S(i).xd, S(i).yd are node coordinates with respect to 
    % Field Dimensions xm and ym variables 
    % end: Yasmin
    S(i).xd=rand(1,1)*xm;
    %XR(i)=S(i).xd;
    S(i).yd=rand(1,1)*ym;
    %YR(i)=S(i).yd;
    S(i).G=0;
    
    % Add: Yasmin
    % S(i).type is the type of a node either
    % 'N' for Node or 'C' for Cluster Head
    % end: Yasmin
    
    %initially there are no cluster heads only nodes
    S(i).type='N'; 
   
    temp_rnd0=i;
    %Random Election of Normal Nodes
    if (temp_rnd0>=m*n+1) 
        S(i).E=Eo;
        S(i).ENERGY=0;
        plot(S(i).xd,S(i).yd,'o');
        hold on;
    end
    %Random Election of Advanced Nodes
    if (temp_rnd0<m*n+1)  
        S(i).E=Eo*(1+a)
        S(i).ENERGY=1;
        plot(S(i).xd,S(i).yd,'+');
        hold on;
    end
end

S(n+1).xd=sink.x;
S(n+1).yd=sink.y;
plot(S(n+1).xd,S(n+1).yd,'x');
    
        
%First Iteration
figure(1);

%counter for CHs
countCHs=0;
%counter for CHs per round
rcountCHs=0;
cluster=1;

countCHs;
rcountCHs=rcountCHs+countCHs;
flag_first_dead=0;

for r=0:1:rmax
    % Add: Yasmin
    % r is the current round
    % end: Yasmin
    r 

  %Operation for epoch
  if(mod(r, round(1/p) )==0)
    for i=1:1:n
        % Add: Yasmin
        % G is an indicator if node i has been CH or not
        % S(i).G = 0 when node i has been CH within the last 1/p rounds, 
        % otherwise S(i).G != 0 (i.e. which 
        % means that node i has not been CH in the last 1/p rounds) 
        % end: Yasmin

        S(i).G=0;  
        %S(i).cl=0;
    end
  end

hold off;

%Number of dead nodes
dead=0;
%Number of dead Advanced Nodes
dead_a=0;
%Number of dead Normal Nodes
dead_n=0;

%counter for bit transmitted to Bases Station and to Cluster Heads
packets_TO_BS=0;
packets_TO_CH=0;
%counter for bit transmitted to Bases Station and to Cluster Heads 
%per round
PACKETS_TO_CH(r+1)=0;
PACKETS_TO_BS(r+1)=0;

figure(1);

for i=1:1:n
    %checking if there is a dead node
    if (S(i).E<=0)
        plot(S(i).xd,S(i).yd,'red .'); % Yasmin: dot nodes are the dead nodes
        dead=dead+1;
        if(S(i).ENERGY==1)
            dead_a=dead_a+1;
        end
        if(S(i).ENERGY==0)
            dead_n=dead_n+1;
        end
        hold on;    
    end
    if S(i).E>0
        S(i).type='N';
        if (S(i).ENERGY==0)  
        % Add: Yasmin
        % circle nodes are alive nodes
        % end: Yasmin
        plot(S(i).xd,S(i).yd,'o'); 
        end
        if (S(i).ENERGY==1) 
        % Add: Yasmin
        % this should be cluster heads ?
        % end: Yasmin
        plot(S(i).xd,S(i).yd,'+'); 
        end
        hold on;
    end
end
plot(S(n+1).xd,S(n+1).yd,'x');


STATISTICS(r+1).DEAD=dead;
DEAD(r+1)=dead;

% Add: Yasmin
% DEAD_N => indicates dead normal nodes
% DEAD_A => indicates dead Advanced nodes
% end: Yasmin
DEAD_N(r+1)=dead_n;  
DEAD_A(r+1)=dead_a;  

%When the first node dies
if (dead==1)
    if(flag_first_dead==0)
        % Add: Yasmin
        % the round where the first node died  
        % end: Yasmin
        first_dead=r   
        flag_first_dead=1;
    end
end

countCHs=0;
cluster=1;
for i=1:1:n
   if(S(i).E>0)
   temp_rand=rand;     
   if ( (S(i).G)<=0)

 %Election of Cluster Heads
 if(temp_rand<= (p/(1-p*mod(r,round(1/p)))))
            countCHs=countCHs+1;
            packets_TO_BS=packets_TO_BS+1;
            PACKETS_TO_BS(r+1)=packets_TO_BS;
            
            S(i).type='C';
            S(i).G=round(1/p)-1;
            C(cluster).xd=S(i).xd;
            C(cluster).yd=S(i).yd;
            plot(S(i).xd,S(i).yd,'k*');
            
            distance=sqrt( (S(i).xd-(S(n+1).xd) )^2 + (S(i).yd-(S(n+1).yd) )^2 );
            C(cluster).distance=distance;
            C(cluster).id=i;
            X(cluster)=S(i).xd;
            Y(cluster)=S(i).yd;
            cluster=cluster+1;
            
            %Calculation of Energy dissipated
            distance;
            if (distance>do)
                S(i).E=S(i).E- ( (ETX+EDA)*(4000) + Emp*4000*( distance*distance*distance*distance )); 
            end
            if (distance<=do)
                S(i).E=S(i).E- ( (ETX+EDA)*(4000)  + Efs*4000*( distance * distance )); 
                
            end
        end     
    
    end
  end 
end

STATISTICS(r+1).CLUSTERHEADS=cluster-1;

% Add: Yasmin
%CLUSTERHS(r+1) => it contains # of Cluster Heads in round r+1
% end: Yasmin

CLUSTERHS(r+1)=cluster-1;   

%Election of Associated Cluster Head for Normal Nodes
for i=1:1:n
   if ( S(i).type=='N' && S(i).E>0 )
     if(cluster-1>=1)
       min_dis=sqrt( (S(i).xd-S(n+1).xd)^2 + (S(i).yd-S(n+1).yd)^2 );
       min_dis_cluster=1;
       for c=1:1:cluster-1
           temp=min(min_dis,sqrt( (S(i).xd-C(c).xd)^2 + (S(i).yd-C(c).yd)^2 ) );
           if ( temp<min_dis )
               min_dis=temp;
               min_dis_cluster=c;
           end
       end
       
       %Energy dissipated by associated Cluster Head
            min_dis;
            if (min_dis>do)
                S(i).E=S(i).E- ( ETX*(4000) + Emp*4000*( min_dis * min_dis * min_dis * min_dis)); 
            end
            if (min_dis<=do)
                S(i).E=S(i).E- ( ETX*(4000) + Efs*4000*( min_dis * min_dis)); 
            end
        %Energy dissipated
        if(min_dis>0)
          S(C(min_dis_cluster).id).E = S(C(min_dis_cluster).id).E- ( (ERX + EDA)*4000 ); 
         PACKETS_TO_CH(r+1)=n-dead-cluster+1; 
        end

       S(i).min_dis=min_dis;
       S(i).min_dis_cluster=min_dis_cluster;
           
   end
 end
end
hold on;

countCHs;
rcountCHs=rcountCHs+countCHs;



%Code for Voronoi Cells
%Unfortynately if there is a small
%number of cells, Matlab's voronoi
%procedure has some problems

[vx,vy]=voronoi(X,Y);
plot(X,Y,'r*',vx,vy,'b-');
hold on;
voronoi(X,Y);
axis([0 xm 0 ym]);


% Add: Yasmin 
% Get avg energy dissipation per round
% Calculate then plot avg. energy dissipation

energy_per_rorund = 0;
for i=1:1:n 
     if S(i).E > 0
         energy_per_rorund = energy_per_rorund + S(i).E;
     end
end 
STATISTICS(r+1).AVG_ENERGY = energy_per_rorund/n
% end: Yasmin


end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   STATISTICS    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                     %
%  DEAD  : a rmax x 1 array of number of dead nodes/round 
%  DEAD_A : a rmax x 1 array of number of dead Advanced nodes/round
%  DEAD_N : a rmax x 1 array of number of dead Normal nodes/round
%  CLUSTERHS : a rmax x 1 array of number of Cluster Heads/round
%  PACKETS_TO_BS : a rmax x 1 array of number packets send to Base Station/round
%  PACKETS_TO_CH : a rmax x 1 array of number of packets send to ClusterHeads/round
%  first_dead: the round where the first node died                   
%                                                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Add: Yasmin 
% Plot results
rmax_range = 1:1:rmax+1;
figure_name = 'plot_results/LEACH';

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
set(gca,'LooseInset',get(gca,'TightInset'));
set(gca,'fontsize',14)

saveas(fig,sprintf('%s_%d',figure_name,1));
saveas(fig,sprintf('%s_%d',figure_name,1),'pdf');

%subplot(2,2,1)
figure('Name',figure_name);
plot(rmax_range, PACKETS_TO_BS);
title('number of packets sent to base station per round');
xlabel('time steps (round)');
ylabel('no of packets sent to base station');
grid on;

set(gca,'XLim',[0 4000]);
set(gca,'XTick',(0:500:4000));

set(gca,'YLim',[0 18])
set(gca,'YTick',(0:2:18))

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
set(gca,'LooseInset',get(gca,'TightInset'));
set(gca,'fontsize',14)

saveas(fig,sprintf('%s_%d',figure_name,2));
saveas(fig,sprintf('%s_%d',figure_name,2),'pdf');

%subplot(2,2,2)
figure('Name',figure_name);
plot(rmax_range, PACKETS_TO_CH);
title('number of packets sent to cluster heads per round');
xlabel('time steps (round)');
ylabel('no of packets sent to cluster heads');
grid on;

set(gca,'XLim',[0 4000]);
set(gca,'XTick',(0:500:4000));

set(gca,'YLim',[0 60])
set(gca,'YTick',(0:10:60))

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
set(gca,'LooseInset',get(gca,'TightInset'));
set(gca,'fontsize',14)


saveas(fig,sprintf('%s_%d',figure_name,3));
saveas(fig,sprintf('%s_%d',figure_name,3),'pdf');

%subplot(2,2,3)
figure('Name',figure_name);
plot(rmax_range, (n-DEAD));
title('number of sensors still alive per round');
xlabel('time steps (round)');
ylabel('no of alive nodes');
ylim([1 55]);
grid on;

set(gca,'XLim',[0 4000]);
set(gca,'XTick',(0:500:4000));

set(gca,'YLim',[0 55])
set(gca,'YTick',(0:5:55))

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
set(gca,'LooseInset',get(gca,'TightInset'));
set(gca,'fontsize',14)

saveas(fig,sprintf('%s_%d',figure_name,4));
saveas(fig,sprintf('%s_%d',figure_name,4),'pdf');

figure('Name',figure_name);
plot([STATISTICS.AVG_ENERGY]);
title('avgerage energy dissipation for nodes per round (J)');
xlabel('time steps (round)');
ylabel('avg. energy of each round (J)');
grid on;

set(gca,'XLim',[0 4000]);
set(gca,'XTick',(0:500:4000));

set(gca,'YLim',[0 0.5])
set(gca,'YTick',(0:0.1:0.5))

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
set(gca,'LooseInset',get(gca,'TightInset'));
set(gca,'fontsize',14)

saveas(fig,sprintf('%s_%d',figure_name,5));
saveas(fig,sprintf('%s_%d',figure_name,5),'pdf');


figure('Name',figure_name);
packet_size =4000;   % in bits
througputs = PACKETS_TO_CH + PACKETS_TO_BS; % total number of packets in network
plot(1:499:rmax+1,througputs(1:499:end)*packet_size/1024,'o-')
xlabel('time steps (round)');
ylabel('throughput in Kbits');
grid on;

set(gca,'XLim',[0 4000]);
set(gca,'XTick',(0:500:4000));

set(gca,'YLim',[0 250])
set(gca,'YTick',(0:50:250))

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
set(gca,'LooseInset',get(gca,'TightInset'));
set(gca,'fontsize',14)

saveas(fig,sprintf('%s_%d',figure_name,6));
saveas(fig,sprintf('%s_%d',figure_name,6),'pdf');

% %subplot(2,2,4)
% figure(6);
% plot(rmax_range, DEAD);
% title('Number of dead sensors per round');
% xlabel('Time steps (round)');
% ylabel('No of dead sensors');


% figure(7);
% plot(CLUSTERHS, [STATISTICS.AVG_ENERGY]*n);
% title('Avgerage energy dissipation per round (J)');
% xlabel('Time steps (round)');
% ylabel('Avg. energy dissipation (J)');

%plot(1:499:10000,PACKETS_TO_BS(1:499:end)*4000/1024,'o-')
%plot(1:499:10000,PACKETS_TO_BS(1:499:end)*4000/1000,'o-')
%plot(1:499:10000,PACKETS_TO_CH(1:499:end)*4000/1000,'o-')
% end: Yasmin

% to know when last node dies
%n-DEAD
%
% to get the throughputs in the entire network 
%sum(througputs(1:499:end)*packet_size/1024);