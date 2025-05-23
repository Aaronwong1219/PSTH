 clear all
close all
LFP_Name='QSS05_2718_16';
Spike_FileName='QSS05_2718_14_1';
max_Rate=50;
N1_Win=0.05;P1_Win=0.05;% time window width to count N1 and P1
Baseline_Window=[-0.5 -0.05];  %PSTH window
  S1_Window=[0 0.3];            % PSTH window for Glu
  S2_Window=[0.58 0.88];        % PSTH window for Glu
%   S1_Window=[0 0.05];            % PSTH window for GABA
%   S2_Window=[0.58 0.63];        % PSTH window for GABA
if length(Spike_FileName)==15
    Select_Ch_ID=str2num(Spike_FileName(end-3:end-2));
else
    Select_Ch_ID=str2num(Spike_FileName(end-2));
end;
No_CH=str2num(LFP_Name(end-1:end));
SmoothingFlag=1;
Min_X=-1;Max_X=1.5;
%Sample_Fre=1000;
Sample_Fre=3051.76;
Delay=0;% time delay of trigger
Ttriger=0;% Triger time
Duration=0.08; % duration of one sound in two-tone
Sound_Interval=0.5; % duration of one sound in two-tone
BinSize=1/Sample_Fre;
Pre_Stimulus_Time=-1;
Post_Stimulus_Time=1.5;
Mean_PSTH=[];LFP=[];LFP_40=[];ersp=[];itc=[];

Neural_Data=func_read_spike_data_sorted(Spike_FileName,Pre_Stimulus_Time,Post_Stimulus_Time);
%if isnan(Spike_FileName)==0
if isempty(Neural_Data)==0
Neural_Data=func_Exclude_burst_spikes(Neural_Data,0.008);
Stimulus=length(Neural_Data);
Y_Fre(1:Stimulus)=0;
% --------------------------------------------------------------------
% plot the raster plot
    subplot(3,2,1)
    for i=1:Stimulus
        Y_Fre(i)=i;
        SpikeTime=Neural_Data(i).SpikeTime;
        if isempty(SpikeTime)==0
        hi=plot(SpikeTime,Y_Fre(i),'ks','MarkerSize',1,'MarkerFaceColor','k');
        set(hi,'MarkerFaceColor','k');hold on;
        end;
        clear SpikeTime;
    end;
    max_freq=max(Y_Fre(i));
    TickLabel_Num=[10:25:150];
    axis([Min_X,Max_X,0,max_freq]);
    set(gca,'TickDir','out','YTick',[10:25:150],'YTickLabel',TickLabel_Num,'FontSize',8);
    AxisVal1=[Ttriger Ttriger NaN Ttriger+Duration+Sound_Interval Ttriger+Duration+Sound_Interval];
    AxisVal2=[0 max_freq NaN 0 max_freq];
    h=line(AxisVal1,AxisVal2);
    set(h,'LineStyle','-','LineWidth',0.5,'color','k'); hold on;
    title(sprintf(' %s ',Spike_FileName),'FontSize',12);
    %xlabel('Time (s)','FontSize',9);
    ylabel('No trial','FontSize',8);

%+++++++++++++++++++++++++++++++++++++++++++++++
%  plot PSTH
[PSTH,TimeVal]=func_Cal_PSTH_Click(Neural_Data,BinSize,SmoothingFlag,Pre_Stimulus_Time,Post_Stimulus_Time);
Mean_PSTH=mean(PSTH);


Mean_PSTH_Baseline=mean(Mean_PSTH(((Baseline_Window(1)+0.5)*Sample_Fre)+1: (Baseline_Window(2)+0.5)*Sample_Fre));
STD_Baseline=std(Mean_PSTH(((Baseline_Window(1)+0.5)*Sample_Fre)+1: (Baseline_Window(2)+0.5)*Sample_Fre));
Mean_PSTH_S1=mean(Mean_PSTH((S1_Window(1)+0.5)*Sample_Fre:(S1_Window(2)+0.5)*Sample_Fre));
Mean_PSTH_S2=mean(Mean_PSTH((S2_Window(1)+0.5)*Sample_Fre:(S2_Window(2)+0.5)*Sample_Fre));
subplot(3,2,3);
%stairs(TimeVal,Mean_PSTH,'k');
h=[1,2,3,4,5,6,7,8,9,10,9,8,7,6,5,4,3,2,1];
Mean_PSTH3(1,:)=func_smooth(Mean_PSTH(1,:),h);
Mean_PSTH3(1,:)=func_smooth(Mean_PSTH3(1,:),h);
Mean_PSTH3(1,:)=func_smooth(Mean_PSTH3(1,:),h);
Mean_PSTH3(1,:)=func_smooth(Mean_PSTH3(1,:),h);
Mean_PSTH3(1,:)=func_smooth(Mean_PSTH3(1,:),h);
Mean_PSTH3(1,:)=func_smooth(Mean_PSTH3(1,:),h);
Mean_PSTH3(1,:)=func_smooth(Mean_PSTH3(1,:),h);
Mean_PSTH3(1,:)=func_smooth(Mean_PSTH3(1,:),h);

% Mean_PSTH3(1,:)=smoothts(Mean_PSTH(1,:),'g');
% Mean_PSTH3(1,:)=smoothts(Mean_PSTH3(1,:),'g',16);


stairs(TimeVal, Mean_PSTH3,'k');
hold on;
% max_Rate=max(Mean_PSTH)*1.2;
axis([Min_X,Max_X,0,max_Rate]);
set(gca,'TickDir','out');
line([Ttriger Ttriger NaN Ttriger+Duration+Sound_Interval Ttriger+Duration+Sound_Interval],[0 max_Rate NaN 0 max_Rate],'LineStyle','-','LineWidth',0.5,'color','k');
title(strcat('基线 ',num2str(Mean_PSTH_Baseline),' 第一声',num2str(Mean_PSTH_S1),'  第二声  ',num2str(Mean_PSTH_S2)))



% Onset_ID_S1=min(find(TimeVal>=0));
% Offset_ID_S1=min(find(TimeVal>=0.1));
% Onset_ID_S2=min(find(TimeVal>=0.58));
% Offset_ID_S2=min(find(TimeVal>=0.68));
Onset_ID_S1=min(find(TimeVal>=S1_Window(1)));
Offset_ID_S1=min(find(TimeVal>=S1_Window(2)));
Onset_ID_S2=min(find(TimeVal>=S2_Window(1)));
Offset_ID_S2=min(find(TimeVal>=S2_Window(2)));
subplot(3,2,2)
% plot(TimeVal(Onset_ID_S1:Offset_ID_S1),Mean_PSTH(Onset_ID_S1:Offset_ID_S1));hold on;
% plot(TimeVal(Onset_ID_S2:Offset_ID_S2),Mean_PSTH(Onset_ID_S2:Offset_ID_S2),'r');
plot(Mean_PSTH3(Onset_ID_S1:Offset_ID_S1));hold on;
plot(Mean_PSTH3(Onset_ID_S2:Offset_ID_S2),'r');
% Axis_Used=max([Offset_ID_S2-Onset_ID_S2,Offset_ID_S1-Onset_ID_S1]);
Axis_Used=max([Offset_ID_S2-Onset_ID_S2,Offset_ID_S1-Onset_ID_S1]);
ylabel('Firing rate');xlabel('post stim time(s)');title('blue S1;red S2');
ylim([min(Mean_PSTH3(Onset_ID_S1:Offset_ID_S1)),1.2*max(Mean_PSTH3(Onset_ID_S1:Offset_ID_S1))]);
xlim([0,Axis_Used]);
% axis([0,0.1,min(Mean_PSTH(Onset_ID_S1:Offset_ID_S1)),1.2*max(Mean_PSTH(Onset_ID_S1:Offset_ID_S1))]);
[P1,LP1]=max(Mean_PSTH3(Onset_ID_S1:Offset_ID_S1+fix(P1_Win/BinSize)));%P1:最大值，LP；最大值对应的横坐标
[P2,LP2]=max(Mean_PSTH3(Onset_ID_S2:Offset_ID_S2+fix(P1_Win/BinSize)));
IG_P=P2/P1;
% title(sprintf('blue S1;red S2;s2/s1=%.2f ',IG_P));
set(gca,'XTick',[0:25*3051.4/1000:Axis_Used],'XTickLabel',num2str([0:25:Axis_Used/3051.4*1000]'),'FontSize',7);


title(strcat('blue S1;red S2;s2/s1= ',num2str(IG_P),'%%S1 latency=',num2str(LP1/3051.76*1000),'ms S2 Latency=',num2str(LP2/3051.76*1000),'ms'));
%line([N1_Win,N1_Win NaN P1_Win,P1_Win],[Min_Amp,Max_Amp NaN Min_Amp,Max_Amp]);
% line([TimeVal(Onset_ID_S1+LP1),TimeVal(Onset_ID_S1+LP1)+0.002],[P1,P1],'Color','k');
% line([TimeVal(Onset_ID_S1+LP2),TimeVal(Onset_ID_S1+LP2)+0.002],[P2,P2],'Color','k');
line([LP1,LP1+20],[P1,P1],'Color','k');hold on
line([LP2,LP2+20],[P2,P2],'Color','k');
subplot(3,1,3)
Z_P1=(P1-Mean_PSTH_Baseline)/STD_Baseline;%Z-score
Z_P2=(P2-Mean_PSTH_Baseline)/STD_Baseline;
IG_Z=Z_P2/Z_P1;
title(strcat('P1 zscore:',num2str(Z_P1),'P2 zscore:',num2str(Z_P2),'IG zscore',num2str(IG_Z)))
% line([TimeVal(Onset_ID_S1+LP2),TimeVal(Onset_ID_S1+LP2)+0.002],[P2,P2],'Color','k');
end

Sorted_Wave_Epoch=func_read_LFP(LFP_Name);
Pre_Stim_Time=-0.5; % start time  of the analysis window relative to the stimulus onset
Post_Stim_Time=1.1;  % end time  of the analysis window
Analysis_Start_Time=0.0;
Analysis_Time_Window=0.5;
Max_Power0=0.25;Max_Power=0.15;
Min_Amp=-2;Max_Amp=2;
%%   
time2=(1:length(Sorted_Wave_Epoch))/Sample_Fre+Pre_Stim_Time+Delay;
subplot(3,2,1) 
[Trials,Len]=size(squeeze(Sorted_Wave_Epoch(Select_Ch_ID,:,:)));%squeeze:删除单维度
LFP= mean(eegfilt(squeeze(Sorted_Wave_Epoch(Select_Ch_ID,:,:)),Sample_Fre,0,300));%eegfilt(数据，采样率，locutoff,hicutoff)
plot(time2,LFP*20+60,'r'); hold on;
title(sprintf('%s, CH:%d,n=%d ',LFP_Name,Select_Ch_ID,Trials),'FontSize',7);
set(gca,'XTick',[0:0.1:Post_Stim_Time],'XTickLabel',num2str([0:0.1:Post_Stim_Time]'),'FontSize',7);
if isempty(Neural_Data)==0
    subplot(3,2,3) 
    plot(time2,LFP*20+max_Rate/1.5,'r');
else
    subplot(3,2,3) 
    plot(time2,LFP,'r');
    axis([Min_X,Max_X,Min_Amp,Max_Amp]);    
end


clear Neural_Data;


Onset_ID_S1=min(find(time2>=S1_Window(1)));
Offset_ID_S1=min(find(time2>=S1_Window(2)));
Onset_ID_S2=min(find(time2>=S2_Window(1)));
Offset_ID_S2=min(find(time2>=S2_Window(2)));
subplot(3,2,4)
plot(LFP(Onset_ID_S1:Offset_ID_S1));hold on;
plot(LFP(Onset_ID_S2:Offset_ID_S2),'r');
ylabel('Amp of LFP');xlabel('post stim time(s)');title('blue S1;red S2');
% axis([0,0.1,1.1*min(LFP(Onset_ID_S1:Offset_ID_S1)),1.2*max(LFP(Onset_ID_S1:Offset_ID_S1))]);

[N1,LN1]=min(LFP(Onset_ID_S1:Offset_ID_S1+fix(N1_Win/BinSize)));
[P1,LP1]=max(LFP(Onset_ID_S1:Offset_ID_S1+fix(P1_Win/BinSize)));
[N2,LN2]=min(LFP(Onset_ID_S2:Offset_ID_S2+fix(N1_Win/BinSize)));
[P2,LP2]=max(LFP(Onset_ID_S2:Offset_ID_S2+fix(P1_Win/BinSize)));
IG_N=N2/N1;
IG_P=P2/P1;
ABS_IG=(P2-N2)/(P1-N1);

% title(sprintf('s2/s1 N= %.2f; P=%.2f ',IG_N,IG_P));
% set(gca,'XTick',[0:25*3051.4/1000:Axis_Used],'XTickLabel',num2str([0:25:Axis_Used/3051.4*1000]'),'FontSize',7);
% xlim([0,Axis_Used]);
title(strcat('s2/s1 N= ',num2str(IG_N),'P=',num2str(IG_P),'ABS-Value=',num2str(ABS_IG),'N1 latency=',num2str(LN1/3051.76*1000),'ms P1 Latency=',num2str(LP1/3051.76*1000),'ms'));

%line([N1_Win,N1_Win NaN P1_Win,P1_Win],[Min_Amp,Max_Amp NaN Min_Amp,Max_Amp]);
line([LN1,LN1+20 NaN LP1,LP1+20],[N1,N1 NaN P1,P1],'Color','k');
line([LN2,LN2+20 NaN LP2,LP2+20],[N2,N2 NaN P2,P2],'Color','k');


% %% plot the result figure
% if isnan(Spike_FileName)==0
%     FigName=strcat('E:\407\IG\Data_Figure\',Location,'_',LFP_Name,'vs','_',Spike_FileName(12:14),'.jpg');
% else
%     FigName=strcat('E:\407\IG\Data_Figure\',Location,'_',LFP_Name,'_',num2str(Select_Ch_ID),'.jpg');
% end;
% print ('-f','-djpeg',FigName);
% close all;

