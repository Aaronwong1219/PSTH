%% analysis the continuous EEG data of one whole session of two channels 
%2021-9-17 增加了cross-frequency分析功能
%2021-9-23 完善了cross-frequency的绘图功能，增加了数据范围截断功能（时长1秒）
%2021-9-28 修复了cross-frequency的光相关错误
%2021-9-29 修复了标题道数错误的bug,增加了对一道的支持
%2021-11-3 完全调用tri文件，field，细胞，动作电位波形均支持
%2021-12-22 绘制field图时去除坏道功能实装（三图全去掉），spike未去掉相应trail
%2022-5-31 ZETA工具包 细胞raster图混匀
%2022-6-27 修复timef函数中关于winsize在Sample_Fre~=1000时的相关bug
%2022-6-27 增加了平均field的滤波功能
%%  plot the sound evoked LFPs in the order of theta power during pre-stimulus period
clear all
close all
File_path='G:\Original_Data\Data_Tri\';
%File_path='//Users/qinling/Documents/MATLAB/Data_Txt/';
%File_path='D:\Documents\O-Documents\CMU_Floor3\Data_Txt\';
% LFP_File='LTR37_1704_16.mat';% file of continuous LFP signals
% LFP_X=13;LFP_X_0=15; Sort_Ch=11;Ref_Ch=15;EMG_Ch=15;
% LKM160204
%LFP_File='LKM18_7104.mat';% file of continuous LFP signal

 LFP_File='QSS05_0710'; 
% LFP_File='QLT01_3302'; 
Spike_Selected={'12_1'};
% Spike_Selected={'1_1','2_1','3_1'};
% Spike_Selected={'1_1','2_1','3_1','5_1','6_1','16_1'};
%,'9_1','10_1','11_1','12_1','13_1','14_1','15_1','16_1'
%'1_1','2_1','3_1','4_1','5_1','6_1','7_1','8_1'
% Pre_Axis_Time=0.2; %指定轴的开始时间
% Post_Axis_Time=1.2;   %指定轴的结束时间
Pre_Axis_Time=-0.5; %指定轴的开始时间
Post_Axis_Time=2;   %指定轴的结束时间
Pre_Stim_Time=2; %最好不要改
Post_Stim_Time=3; %最好不要改
Low_Plot_Threshold=-150;%field绘图下界
High_Plot_Threshold=100;%field绘图上界
Ersp_Freq=100;%热图的纵轴（频率）最大值
MinY=-10;MaxY=100;%调整绘图轴  60
cross_frequency_analysis_tag=0;  %是否进行cross-frequency分析
Sort_Step2=0;
Time_Window=0.2;%?纳入分析时间窗的长度，从声音开始时刻算起？
Sort_Frequency=[60,80];%然后比找到的最大值？
Sort_Channel=10;%选择要排序的道数
threshold_artifact=2000; %肌电滤过阈值
global eeglab_only_ersp
eeglab_only_ersp=1;%0有热度趋势条 1是没有
Light_On=0;  %是否开光
erspmax=4;% 热图能量范围
Disorder_Control=0;%细胞放电散点图是否进行小规模重排序 1开启 0不开启
mean_filter=0;  %是否进行平均field滤波,0为关闭，1为开启
mean_low_filter=35; %平均field滤波下限
mean_high_filter=45;%平均field滤波上限
Smooth_SD=1;% PSTH 平滑指数 默认为5
global AbNormal_Point
AbNormal_Point=0;  %PSTH开头强行归零 0为不开启，1为开启该功能
%%    后面不用改
Pre_Axis_Time=-1*abs(Pre_Axis_Time);%防bug用。保证其为负值
rmpath('D:\MyDocuments\Work\download_toolbox\buzcode-master\externalPackages\eeglab14_0_0b\functions\timefreqfunc')
% cd('D:\MyDocuments\Work')

LFP_File_Seg=cat(2,LFP_File,'tri.mat');% file of segmented LFP signals of each stimulus trial
% load LFP data of mat
load(cat(2,File_path,LFP_File_Seg));
for i=1:length(Spike_Selected)
   Signal_Ch_Original(i)= str2num(Spike_Selected{i}(1:find(Spike_Selected{i}=='_')-1));
    Signal_Ch(i)=find(Signal_Ch_Original(i)==LFP_array);
end

try
Wave_Epoch=func_cut_from_continous_data2(Wave_Epoch_continue,Pre_Stim_Time,Post_Stim_Time,Stimuli_Time,Sample_Fre_Tri);
catch
    Stimuli_Time=Stimuli_Time(1:end-1);
Wave_Epoch=func_cut_from_continous_data2(Wave_Epoch_continue,Pre_Stim_Time,Post_Stim_Time,Stimuli_Time,Sample_Fre_Tri);

end

% Wave_Epoch(:,:,round(2.3796*Sample_Fre_Tri):round(2.4982*Sample_Fre_Tri)) = Wave_Epoch(:,:,round(2.3796*Sample_Fre_Tri):round(2.4982*Sample_Fre_Tri))*0.6;


% Wave_Epoch=10*Wave_Epoch;
Stim_Count=size(Wave_Epoch,2);
%Start_Trial=1;End_Trial=60;
%Start_Trial=241;End_Trial=300;
Start_Trial=1; End_Trial=Stim_Count;
% No_CH=str2num(LFP_File(12:13)); % total count of channels
Min_Color=0;Max_Color=6; % color range of spectrum
% fmin=0;fmax=30; % frequency range of FFT
max_y=0.05; %control the amplitude of LFP wave.
Sample_Fre=Sample_Fre_Tri;BinSize=1/Sample_Fre;
%% sort the LFP by light or without light 
load('D:\Work\light_click_train_control.mat')
% if Mouse_Monitor==1
%     Stimuli_Time=Stimuli_Time(2:end);
% LFP=Wave_Epoch(:,2:end,:);
% else
LFP=Wave_Epoch;
% end
light_click_train_control(121:240)=light_click_train_control;
light_click_train_control=light_click_train_control(1:size(LFP,2));
[M,IX]=sort(light_click_train_control);

light_point=find(M==1,1);
for i=1:length(IX)
    temp(:,i,:)=LFP(:,IX(i),:);
end
if length(size(temp))==2
 tempnew(1,:,:)=temp';
 temp=tempnew;
end
clear LFP;LFP=temp;clear temp;
 %%    step-2 sort the light and without light by frequency power
 switch Sort_Step2
     case 0
         
     case 1
Temp_Light=LFP(:,light_point:end,:);
Temp_Without_Light=LFP(:,1:light_point-1,:);
%把要画的几道的LFP进行顺序重拍，用chronux工具箱
% window=5;
% movingwin= [window,1];   %(in the form [window winstep] i.e length of moving window and step size.  Note that units here have  to be consistent with
%params.Fs=1000; %(sampling frequency) - optional. Default 1.
params.FsSample_Fre_Tri; %(sampling frequency) - optional. Default 1.
params.fpass=Sort_Frequency;
params.tapers=[5,3];params.tapers=dpsschk(params.tapers,Time_Window*1000,1000);
params.err=[1,0];
params.pad=0;
params.trialave=0;%0为不进行trial间平均
data=[];data=Temp_Without_Light(Sort_Channel,:,abs(Time_Window-(Pre_Stim_Time+0.45))*1000+1:(Pre_Stim_Time+0.45)*1000);data=squeeze(data);data=data';
[S_Without_Light,f,~]=mtspectrumc(data,params);  %单窗法计算各频率段能量
S_Without_Light=mean(S_Without_Light,1);
[~,B]=sort(S_Without_Light);
data=[];data=Temp_Light(Sort_Channel,:,abs(Time_Window-(Pre_Stim_Time+0.45))*1000+1:(Pre_Stim_Time+0.45)*1000);data=squeeze(data);data=data';
[S_Light,f,~]=mtspectrumc(data,params);  %单窗法计算各频率段能量
S_Light=mean(S_Light,1);
[~,B2]=sort(S_Light);
B2=B2+max(B);IX=cat(2,B,B2); %制作出来新的排序
LFP=LFP(:,IX,:);%根据chronux 工具箱结果再次调换LFP排序   
     case 2
 Temp_Light=LFP(:,light_point:end,:);
Temp_Without_Light=LFP(:,1:light_point-1,:);   
data=[];data=Temp_Without_Light(Sort_Channel,:,abs(Time_Window-(Pre_Stim_Time+0.45))*1000+1:(Pre_Stim_Time+0.45)*1000);data=squeeze(data);data=data';
data=eegfilt(data',1000,Sort_Frequency(1),Sort_Frequency(2));
for i=1:size(data,1)
    temp_data=data(i,:);
  data2=hilbert(temp_data);
  data_sort(i)=max(abs(data2));
end
[A,B]=sort(data_sort);
data=[];data=Temp_Light(Sort_Channel,:,abs(Time_Window-(Pre_Stim_Time+0.45))*1000+1:(Pre_Stim_Time+0.45)*1000);data=squeeze(data);data=data';
data=eegfilt(data',1000,Sort_Frequency(1),Sort_Frequency(2));
for i=1:size(data,1)
    temp_data=data(i,:);
  data2=hilbert(temp_data);
  data_sort(i)=max(abs(data2));
end
[A,B2]=sort(data_sort);
B2=B2+max(B);IX=cat(2,B,B2); %制作出来新的排序
LFP=LFP(:,IX,:);%根据eegfilt hilbert变换结果再次调换LFP排序
% data2=hilbert(data);
% [A,B2]=sort(data_sort);
 end
%  ID=[1:Stim_Count];      
%% plot the mean LFP together
time3=(-1*Pre_Stim_Time:1/Sample_Fre:Post_Stim_Time);
if length(LFP)<length(time3)
    time3=time3(1:length(LFP));
end
    
%% plot the LFP of each trials
Y_Scale=MaxY;
figure(1)
% if Mouse_Monitor==1
% End_Trial=End_Trial-1;
% end
MinY=Y_Scale*Start_Trial;
MaxY=Y_Scale*End_Trial;
delete_number=zeros(length(Signal_Ch),1);
for j=1:length(Signal_Ch)
    for i=Start_Trial:End_Trial
        subplot(4,length(Signal_Ch),[j j+length(Signal_Ch)])
        Wave_Trial=squeeze(LFP(Signal_Ch(j),i,:));
        if max(abs(Wave_Trial))<=threshold_artifact
        plot(time3,Wave_Trial+Y_Scale*(i-1),'b'); hold on;
%              plot(time3,Wave_Trial,'b'); hold on;
        else
            delete_number(j)=delete_number(j)+1;
        end
    end
%     axis([-1*Pre_Stim_Time,Post_Stim_Time,MinY,MaxY]); 
     axis([Pre_Axis_Time,Post_Axis_Time,MinY,MaxY]); 
     line([0.22,0.22 ],[MinY,MaxY ],'color','r');
    line([0.45,0.45 NaN 0.95,0.95],[MinY,MaxY NaN MinY,MaxY],'color','k');
    if j==1
            title(sprintf('%s--Ch:%d%s;',LFP_File,Signal_Ch_Original(j),'sort by light'));
    else
            title(sprintf('%d;',Signal_Ch_Original(j)));
    end;
    set(gca,'yticklabel','');
end
% subplot(1,length(Signal_Ch),1)
%  subplot(4,length(Signal_Ch),[1 1+length(Signal_Ch)])
% xlabel('Time');ylabel('Trial');
%%    light and without light mean value

% eegfilt(squeeze(Sorted_Wave_Epoch(CH,:,StartBin:EndBin)),Sample_Fre,locutoff,hicutoff);
% squeeze(LFP(Signal_Ch(1),light_point:end,:)),
for i=1:length(Signal_Ch)
    subplot(4,length(Signal_Ch),i+3*length(Signal_Ch))
%      plot(time3,squeeze(mean(LFP(Signal_Ch(i),1:light_point-1,:),2)),'b')
if Light_On==1
if mean_filter==0   
plot(time3,squeeze(mean(LFP(Signal_Ch(i), find(max(abs(squeeze(LFP(Signal_Ch(i),1:light_point-1,:))'))<=threshold_artifact),:),2)),'b')
else
    temp_select=find(max(abs(squeeze(LFP(Signal_Ch(i),1:light_point-1,:))'))<=threshold_artifact);
    plot(time3,mean(eegfilt(squeeze(LFP(Signal_Ch(i),temp_select,:)),Sample_Fre,mean_low_filter,mean_high_filter),1),'b')
end
    hold on
% %      plot(time3,squeeze(mean(LFP(Signal_Ch(i),light_point:end,:),2)),'r')
%           plot(time3,squeeze(LFP(Signal_Ch(i), light_point-1+find(max(abs(squeeze(LFP(Signal_Ch(i),light_point:end,:))'))<=threshold_artifact),:)),'k')


if mean_filter==0   
plot(time3,squeeze(mean(LFP(Signal_Ch(i), light_point-1+find(max(abs(squeeze(LFP(Signal_Ch(i),light_point:end,:))'))<=threshold_artifact),:),2)),'r')
else
    temp_select=light_point-1+find(max(abs(squeeze(LFP(Signal_Ch(i),light_point:end,:))'))<=threshold_artifact);
    plot(time3,mean(eegfilt(squeeze(LFP(Signal_Ch(i),temp_select,:)),Sample_Fre,mean_low_filter,mean_high_filter),1),'r')
end
else
plot(time3,squeeze(mean(LFP(Signal_Ch(i), find(max(abs(squeeze(LFP(Signal_Ch(i),1:size(LFP,2),:))'))<=threshold_artifact),:),2)),'b')


end
hold on

      line([0.22,0.22 ],[Low_Plot_Threshold,High_Plot_Threshold],'color','r');
    line([0.45,0.45 NaN 0.95,0.95],[Low_Plot_Threshold,High_Plot_Threshold NaN Low_Plot_Threshold,High_Plot_Threshold],'color','k');
     xlim([Pre_Axis_Time,Post_Axis_Time])
     ylim([Low_Plot_Threshold,High_Plot_Threshold])
     if i==1
legend('Without lLight','With Light','Location','southeast')
     end
     xlim([Pre_Axis_Time Post_Axis_Time])
     end
%% parameters of egglab
figure(1)
% tlimits=[-100 1100]; %还要改
% tlimits=[0 2000];
% Pre_Axis_Time
% Post_Axis_Time
winsize=200;%改这行造成一切后果自负！（如果想热图横向平滑些就把值改大）
cycles=[0.5 0.5];%改这行造成一切后果自负！
% Pre_Axis_Time=0;
% Post_Axis_Time=1;
for i=1:length(Signal_Ch)
    tlimits=[0 size(LFP,3)/Sample_Fre*1000]; %还要改
    frames=size(LFP,3);
%     winsize=256;
    % calculate the baseline
data=[];
    Select_Ch_ID=Signal_Ch(i);
  for trial=1:light_point-1    %计算热图的基础值用
      if max(abs(squeeze(LFP(Select_Ch_ID,trial,:))))<=threshold_artifact
         data=cat(1,data,squeeze(LFP(Select_Ch_ID,trial,:)));
      end    
  end
     if isempty(data)
         warning('threshold artifact is too small!!!! lzj')
         continue
     end
%   subplot(8,length(Signal_Ch),i+5*length(Signal_Ch))
    [~,~,powbase,times,~,~,~,~] = timef(data',frames,tlimits,Sample_Fre,cycles,'maxfreq',Ersp_Freq,'baseline',1000*abs(Pre_Stim_Time),'erspmax',6,'itcmax',1,'plotitc','off' ,'plotersp','off','winsize',floor(winsize*Sample_Fre_Tri/1000));  
%取消推迟光遗传固有延迟（200ms）
    tlimits=[0 (Post_Axis_Time-Pre_Axis_Time)*1000+winsize]; %还要改
% frames=(tlimits(2)-tlimits(1))*Sample_Fre;
%  1000*( Post_Axis_Time-Pre_Stim_Time)+100
%  1000*(Pre_Axis_Time-Pre_Stim_Time)-winsize/2+1;
%以上为确定baseline部分
 data=[];
    Select_Ch_ID=Signal_Ch(i);
  for trial=1:light_point-1  
   if max(abs(squeeze(LFP(Select_Ch_ID,trial,Sample_Fre_Tri*(Pre_Axis_Time+Pre_Stim_Time)-winsize/2*Sample_Fre_Tri/1000+1: Sample_Fre_Tri*( Post_Axis_Time+Pre_Stim_Time)+winsize/2*Sample_Fre_Tri/1000)))+1)<=threshold_artifact
   data=cat(1,data,squeeze(LFP(Select_Ch_ID,trial,round(Sample_Fre_Tri*(Pre_Axis_Time+Pre_Stim_Time)-winsize/2*Sample_Fre_Tri/1000+1): floor(Sample_Fre_Tri*( Post_Axis_Time+Pre_Stim_Time)+winsize/2*Sample_Fre_Tri/1000))));
   if trial==1
   frames=length(data);
   end
   end
  end 
    subplot(8,length(Signal_Ch),i+5*length(Signal_Ch))
    [~,~,~,~,~,~,~,~] = timef(data',frames,tlimits,Sample_Fre,cycles,'maxfreq',Ersp_Freq,'baseline',0,'erspmax',erspmax,'itcmax',1,'plotitc','off' ,'winsize',floor(winsize*Sample_Fre_Tri/1000),'powbase',powbase,'timesout',2000);  
    data=[];
    for trial=light_point:End_Trial  %以low_trail的热图基础值作为high_trail的
%      if max(abs(squeeze(LFP(Select_Ch_ID,trial,Sample_Fre*(Pre_Axis_Time+Pre_Stim_Time)-winsize/2+1: Sample_Fre*( Post_Axis_Time+Pre_Stim_Time)+100))))<=threshold_artifact
%     data=cat(1,data,squeeze(LFP(Select_Ch_ID,trial,Sample_Fre*(Pre_Axis_Time+Pre_Stim_Time)-winsize/2+1: Sample_Fre*( Post_Axis_Time+Pre_Stim_Time)+100)));
    
 if max(abs(squeeze(LFP(Select_Ch_ID,trial,Sample_Fre_Tri*(Pre_Axis_Time+Pre_Stim_Time)-winsize/2*Sample_Fre_Tri/1000+1: Sample_Fre_Tri*( Post_Axis_Time+Pre_Stim_Time)+winsize/2*Sample_Fre_Tri/1000)))+1)<=threshold_artifact
   data=cat(1,data,squeeze(LFP(Select_Ch_ID,trial,round(Sample_Fre_Tri*(Pre_Axis_Time+Pre_Stim_Time)-winsize/2*Sample_Fre_Tri/1000+1): floor(Sample_Fre_Tri*( Post_Axis_Time+Pre_Stim_Time)+winsize/2*Sample_Fre_Tri/1000))));
   end
    end
        subplot(8,length(Signal_Ch),i+4*length(Signal_Ch))
    [~,~,~,~,~,~,~,~] = timef(data',frames,tlimits,Sample_Fre,cycles,'maxfreq',Ersp_Freq,'baseline',0,'erspmax',erspmax,'itcmax',1,'plotitc','off' ,'winsize',floor(winsize*Sample_Fre_Tri/1000),'powbase',powbase,'timesout',2000);  
%     [~,~,~,~,~,~,~,~] = timef(data',frames,tlimits,Sample_Fre,cycles,'maxfreq',75,'baseline',100,'erspmax',6,'itcmax',1 ,'plotitc','off','winsize',winsize);  

    clear powbase
end
if Stim_Count<58
    suptitle(strcat('注意：声音刺激仅有',num2str(Stim_Count),'次'))
else
    try
%     suptitle(strcat('threshold_artifact',num2str(threshold_artifact),'delete',num2str(delete_number),'个trail')) %添加删除了多少个trail的信息
    catch
    end
    end
% return
%%   进行细胞相关分析
clearvars -except  Spike_Continous Disorder_Control Smooth_SD LFP_array Stimuli_Time cross_frequency_analysis_tag Mouse_Control  light_click_train_control Spike_Selected IX LFP_File Signal_Ch Totalnew Pre_Stim_Time Post_Stim_Time light_point SpikeWavefom  Spike_Continous  Pre_Axis_Time Post_Axis_Time
BinSize=5/1000;
PLX_FileName=strcat(LFP_File(1:5),LFP_File(7:10));
% PLX_FileName='LTR490904-01';
% FilePath='D:\Plexon Data\L\';
% Pre_Stim_Time=2; % start time  of the analysis window relative to the stimulus onset
% Post_Stim_Time=3;  % end time  of the analysis window
Duration=0.5;
% filename=strcat(FilePath,PLX_FileName,'.plx');%changed lzj
CH_number=16;
% Ev_CH=1;
% [G_SpikeData,N_Spk,G_Sort_Code,SOUNTstim]=func_read_PLX(filename,CH_number,Ev_CH);
% if isempty(Mouse_Control)==0
%     SOUNTstim=SOUNTstim(2:end);
% end
SOUNTstim=Stimuli_Time;
NoLight_Stim=SOUNTstim(find(light_click_train_control==0));
Light_Stim=SOUNTstim(find(light_click_train_control==1));
StimCount=length(Stimuli_Time);
%%  小范围乱序消除时间对细胞放电的影响
if Disorder_Control==1
  [~,Fake_Order]=sort(rand(1,light_point-1));
IX(1:light_point-1)=IX(Fake_Order);
  [~,Fake_Order]=sort(rand(1,length(IX)+1-light_point));
Fake_Order=Fake_Order+light_point-1;
IX(light_point:end)=IX(Fake_Order);
end
figure(2)
Fig_ID=0;
plot_control=0;
%%  repeated the raster plot for each channel
for CH_ID=1:CH_number
    Sort_Number=5;  % number of Sorted units
    if Sort_Number>0
        if CH_ID<10
            CH_Name=strcat('0',num2str(CH_ID));
        else
            CH_Name=num2str(CH_ID);
        end;
       %% repeated the raster plot for each sorted unit
       for Selected_No=1:Sort_Number

          
            for qop=1:length(Spike_Selected)
      if   strcmp(strcat(num2str(CH_ID),'_',num2str(Selected_No)),Spike_Selected{qop})==1
          plot_control=1;
          break
      end
            end  %只画选中的细胞
            if plot_control==1
                Fig_ID=Fig_ID+1;
                plot_control=0;
            else
                continue
            end
            %%
            try
            eval(['SpikeData=Spike_Continous.cell',num2str(CH_ID),'_',num2str(Selected_No),';']);
            catch
                continue
            end
            SpikeCount=length(SpikeData);
            subplot(3,length(Spike_Selected),[qop,length(Spike_Selected)+qop]);
            for i=1:StimCount
               Index1=max(find(SpikeData<=SOUNTstim(i)-Pre_Stim_Time));
               Index2=max(find(SpikeData<=SOUNTstim(i)+Post_Stim_Time));
               if isempty(Index1)==0 && isempty(Index2)==0
               hi=plot(SpikeData(Index1:Index2)-SOUNTstim(i),find(IX==i),'ks','MarkerSize',2,'MarkerFaceColor','k');hold on;
       set(hi,'MarkerFaceColor','k');
               end
               Neural_Data(find(IX==i)).SpikeTime=SpikeData(Index1:Index2)-SOUNTstim(i);
               
            end;
            axis([Pre_Axis_Time,Post_Axis_Time,0 StimCount+1]);
            AxisVal1=[0.45 0.45 NaN 0.45+Duration 0.45+Duration];
            AxisVal2=[0 StimCount NaN 0 StimCount];
            h=line(AxisVal1,AxisVal2);
            line([0.22 0.22],[0 StimCount],'Color','r')
            set(h,'LineStyle','-','LineWidth',0.5);
            title(sprintf(' %s, CH:%s(%d); ',PLX_FileName(1:end),CH_Name,Selected_No));
         
            ylabel('Repetition rate (Hz)','FontSize',8);
            
              %%  绘制细胞放电形状图
%       
%       filename1=strcat(filename(1:end-4),'-01','.plx');
% if ~exist(filename1, 'file')
%     FileName=filename;
% else
%     FileName=filename1;
% end;
%        [~, ~, ~, wave] =plx_waves_v( FileName,CH_ID,Selected_No);
%        subplot(5,2*length(Spike_Selected),2*Fig_ID)
%        plot(wave','k');
hold on
try
  eval(['wave=SpikeWavefom.wave',num2str(CH_ID),'_',num2str(Selected_No),;]);
%        plot(mean(wave,1),'r');axis off
% Pre_Axis_Time=0;
% Post_Axis_Time=1;
wave=wave';
       ratio= StimCount*0.1/(max(mean(wave,1)));
        plot(-1*Pre_Axis_Time:(Post_Axis_Time+Pre_Axis_Time)/29:Post_Axis_Time,mean(wave,1)*ratio+StimCount*0.75,'r');
      title(strcat(PLX_FileName(1:end),'-',CH_Name,'-',num2str(Selected_No),'n=',num2str(SpikeCount)))
       clear wave
catch
end
       %%
       %% 
%        if Mouse_Monitor==1
        subplot(3,length(Spike_Selected), 2*length(Spike_Selected)+qop)     
     [vecTime,vecRate]=getIFR(SpikeData,NoLight_Stim-abs(Pre_Axis_Time),abs(Pre_Axis_Time)+abs(Post_Axis_Time),Smooth_SD);   
      hold on
      getIFR(SpikeData,Light_Stim-abs(Pre_Axis_Time),abs(Pre_Axis_Time)+abs(Post_Axis_Time),Smooth_SD);   
            xlim([0 abs(Pre_Axis_Time)+abs(Post_Axis_Time)])
            %%  PSTH部分
%        Neural_Data_without_light=Neural_Data(1:light_point);
%        Neural_Data_with_light=Neural_Data(light_point+1:end);
%          subplot(3,length(Spike_Selected), 2*length(Spike_Selected)+qop)
%     [PSTH,TimeVal]=func_Cal_PSTH_Click(Neural_Data_without_light,BinSize,1,-1*Pre_Stim_Time,Post_Stim_Time);
%  Mean_PSTH=mean(PSTH);
%         Timenew=[TimeVal(1):0.0001:TimeVal(end)];
% yi= interp1(TimeVal,Mean_PSTH,Timenew,'linear');
%         stairs(TimeVal,Mean_PSTH,'b','LineWidth',2);hold on;
%    [PSTH,TimeVal]=func_Cal_PSTH_Click(Neural_Data_with_light,BinSize,1,-1*Pre_Stim_Time,Post_Stim_Time);
%  Mean_PSTH=mean(PSTH);
%         Timenew=[TimeVal(1):0.0001:TimeVal(end)];
% yi= interp1(TimeVal,Mean_PSTH,Timenew,'linear');
%     stairs(TimeVal,Mean_PSTH,'r','LineWidth',2);hold on;
%             xlim([Pre_Axis_Time,Post_Axis_Time])


       
       end
    end
end



% %% plot the spectrum averaged across the trials
% figure(2)
% %% parameters of egglab
% Low_Trial=fix(End_Trial/4);High_Trial=End_Trial-fix(End_Trial/4);
% clear data;
% data=squeeze(LFP(Signal_Ch(1),:,:));
% [Trial,Len]=size(data);
% frames=Len;tlimits=[0 Len];winsize=256;cycles=[0.5 0.5];
% for j=1:length(Signal_Ch)
%     data=[];
%     Select_Ch_ID=Signal_Ch(j);
%     for trial=1:Low_Trial  
%          data=cat(1,data,squeeze(LFP(Select_Ch_ID,trial,:)));
%     end
%     subplot(2,length(Signal_Ch),length(Signal_Ch)+j)
%     [ersp,itc,powbase,times,freqs,erspboot,itcboot,itcphase] = timef(data',frames,tlimits,Sample_Fre,cycles,'maxfreq',75,'baseline',500,'erspmax',6,'itcmax',2 );  
% 
%     for trial=High_Trial:End_Trial  
%     data=cat(1,data,squeeze(LFP(Select_Ch_ID,trial,:)));
%     end
%     subplot(2,length(Signal_Ch),j)
%     [ersp,itc,powbase,times,freqs,erspboot,itcboot,itcphase] = timef(data',frames,tlimits,Sample_Fre,cycles,'maxfreq',75,'baseline',500,'erspmax',6,'itcmax',2 );  
% end;
% 


% %FigName=strcat('E:\407\Data_Figure\',LFP_File(1:end-4),'.jpg');
% FigName=strcat('/Volumes/KINGSTON/407/Data_Figure/',LFP_File(1:end-4),'-',num2str(LFP_X),'-',num2str(Sort_Ch),'.jpg');
% print (figure(1),'-djpeg',FigName);
% close all;

%%   3 cross-frequency分析  （optional）

if cross_frequency_analysis_tag==1
else
    return
end
figure()
addpath('D:\MyDocuments\Work\download_toolbox\cross-frequency')
addpath('D:\My Documents\MATLAB\Work\lzj\download_toolbox\cross-frequency')
 FileName=LFP_File; % PPC HIP文件
 temp=find(FileName=='_');temp=temp(end);temp=str2num(FileName(temp+1:end));
%  addpath(genpath('D:\Work\dashijie_analysis\eeglab13_4_4b'))
%  LFP_array=[1:temp];%多通道电极排序，不是多通道电极填[1:16]，就可以
Bad_Ch=[];%不使用多通道电极的情况下空着不填，多通道电极有坏的就填
low_phase_range=[1,20];%横轴对应的低频范围上下界
high_amplitude_range=[21,150];%纵轴对应的高频范围上下界
%注意以上两者不能有重叠部分，否则会报错
field_select_laminar=Signal_Ch;%选择amplitude的道数来源，可以选择多个，就会出多个图
%纵轴高频数据来源的来源道数/层数
phase_select_laminar=Signal_Ch;%选择phase的道数来源
%横轴低频数据的来源道数/层数
%以上两个变量代表多通道电极从上到下的层数，不是多通道电极的时候直接填道数就可以
 Analysis_Start_Time=0;%分析时间起始时间（秒）
 Analysis_End_Time=9;%分析时间结束时间（秒）
%不能太长，太长会很慢，一般有个一两分钟就行了
% Mouse_Monitor=1;
Channels=length(LFP_array)-length(Bad_Ch);
%   targetfile=strcat('D:\Original_Data\Data_Txt\',FileName,'Seg.mat');
%  load(targetfile,'Stimuli_Time');
  targetfile=strcat('D:\Original_Data\Data_Tri\',FileName,'tri.mat');
  load(targetfile,'Wave_Epoch_continue');
Sorted_Wave_Epoch=Wave_Epoch_continue; %从连续数据中重新切seg数据
% Stimuli_Time
clear Wave_Epoch temp
% Sorted_Wave_Epoch=Sorted_Wave_Epoch(:,3,:);
for i=1:length(Stimuli_Time)
temp(:,i,:)=Sorted_Wave_Epoch(:,Stimuli_Time(i)*1000+201:Stimuli_Time(i)*1000+201+1000*Analysis_End_Time);
end
Sorted_Wave_Epoch=temp;
% if Mouse_Monitor==1
%     Sorted_Wave_Epoch=Sorted_Wave_Epoch(:,2:end,:);
% %     Stimuli_Time=Stimuli_Time(2:end);
% end
clear temp
%% parameters of egglab
 frames=4579;
 tlimits=[0 1500.4];
 winsize=256;
 %cycles=0; 
 cycles=[0.5 0.5];
locutoff=0;
hicutoff=500; 
Pre_Stim_Time=0; % start time  of the analysis window relative to the stimulus onset
Post_Stim_Time=0;  % end time  of the analysis window
Sample_Fre=1000;
BinSize=1/Sample_Fre;
time2=(1:length(Sorted_Wave_Epoch))/Sample_Fre+Pre_Stim_Time;
StartBin=min(find(time2>Analysis_Start_Time));
EndBin=min(find(time2>=time2(end)));
k=0;
%  for i=1:length(LFP_array)
%     CH= LFP_array(i);
% %      filerted_wave(i,:,:)= eegfilt(squeeze(Sorted_Wave_Epoch(CH,:,StartBin:EndBin))',Sample_Fre,locutoff,hicutoff);
%  filerted_wave(i,:,:)= squeeze(Sorted_Wave_Epoch(CH,:,StartBin:EndBin))';
%  end
 %%     光遗传参数相关
 Sorted_Wave_Epoch_Light=[];
 Sorted_Wave_Epoch_NoLight=[];
%  load('light_click_train_control.mat')
 for i=1:size(Sorted_Wave_Epoch,2)
     if light_click_train_control(i)==1
         Sorted_Wave_Epoch_Light=cat(2,Sorted_Wave_Epoch_Light,Sorted_Wave_Epoch(:,i,:));
     else
         Sorted_Wave_Epoch_NoLight=cat(2,Sorted_Wave_Epoch_NoLight,Sorted_Wave_Epoch(:,i,:));
     end
 end
%  q=squeeze(mean(filerted_wave,2));
MI_total_total_light=zeros(size(Sorted_Wave_Epoch_Light,2),length(phase_select_laminar),length(field_select_laminar),low_phase_range(2),high_amplitude_range(2));
for trial_number=1:size(Sorted_Wave_Epoch_Light,2)
    if length(size(Sorted_Wave_Epoch_Light))==2
  q=squeeze(Sorted_Wave_Epoch_Light(:,trial_number));   
  q=q';
    else
 q=squeeze(Sorted_Wave_Epoch_Light(:,trial_number,:)); 
    end
 %%    调整排序，去除坏道
  if length(Bad_Ch)>0
 for i=1:length(Bad_Ch)
 temp(i)=find(LFP_array==Bad_Ch(i));
end
temp=sort(temp,'descend');
Bad_Ch=LFP_array(temp);
clear temp
 end
 if length(Bad_Ch)>0
     clear CH
     for i=1:length(Bad_Ch)
         CH(i)= find(LFP_array==Bad_Ch(i));  
     end
     q(CH,:)= [];
 end
 
 
 %%   有光部分分析
% q=q(:,40000:43500);
%% Generate the two signals that are interacting
% addpath('D:\My Documents\MATLAB\Work\lzj\function\cross-frequency')
% for phase_laminar_choose=1:length(phase_select_laminar)
for field_laminar_choose=1:length(field_select_laminar)
    
%    if size(q,2)==1
x_field=q(field_select_laminar(field_laminar_choose),:)';%amplitude的数据来源
x_phase=q(phase_select_laminar(field_laminar_choose),:)';%phase的数据来源

 for hengzhou=low_phase_range(1):low_phase_range(2)    %设定横轴频率范围
    for zongzhou=high_amplitude_range(1):high_amplitude_range(2)%   设定纵轴频率分析范围
    f1=hengzhou;
    f2=zongzhou;
[b, a, fCenterList, nTap, fEdges] = filterBank_prop(f1, f2, 1, 3, Sample_Fre, 'fir1', f2-f1);
% [x_filtered, x_analytic] = applyFilterBankThenHT_filtfilt(b, a, x_phase, nTap);
[~, x_analytic] = applyFilterBankThenHT_filtfilt(b, a, x_phase, nTap);
x_analytic=x_analytic(1:445,:);
%MIT程序包里带的两个程序，第一个没看明白啥意思，好像是选择带通滤波的相关参数

% x_analytic(2,:)=

% amplitude = abs(x_analytic);
phase = angle(x_analytic);
[~, x_analytic] = applyFilterBankThenHT_filtfilt(b, a, x_field, nTap);
x_analytic=x_analytic(1:445,:);
amplitude = abs(x_analytic);
%希尔伯特变换相关
phase_sequence=linspace(min(phase(:,1)),max(phase(:,1)),19);
AAA=zeros(1,18); count=zeros(1,18); 
%360度角分成18份

for j=1:length(phase)
for i=1:18
    if phase(j,1)>=phase_sequence(i) && phase(j,1)<=phase_sequence(i+1)
        count(i)= count(i)+1;
        AAA(i)=AAA(i)+amplitude(j,2);
        break
    end
end
end
AAA=AAA./count;%18个区间每个区间的amplitude平均值
% bar(AAA)
H=0;
for i=1:length(AAA)
   H=H+ AAA(i)/sum(AAA)*log(AAA(i)/sum(AAA));
end
H=-H;Hmax=log(length(AAA));
MI=(Hmax-H)/Hmax;%用这个评价两频率的相关性
MI_total(hengzhou,zongzhou)=MI;
    end
end

%%  后面画最强相关性点的例子各角度振幅平均值

% temp=max(MI_total);temp=find(temp==max(temp));
% temp2=max(MI_total');temp2=find(temp2==max(temp2));
% for hengzhou=temp2:temp2
%     for zongzhou=temp:temp
%     f1=hengzhou;
%     f2=zongzhou;
% [b, a, fCenterList, nTap, fEdges] = filterBank_prop(f1, f2, 1, 3, Sample_Fre, 'fir1', f2-f1);
% [~, x_analytic] = applyFilterBankThenHT_filtfilt(b, a, x_phase, nTap);
% 
% % x_analytic(2,:)=
% 
% % amplitude = abs(x_analytic);
% phase = angle(x_analytic);
% [~, x_analytic] = applyFilterBankThenHT_filtfilt(b, a, x_field, nTap);
% amplitude = abs(x_analytic);
% phase_sequence=linspace(min(phase(:,1)),max(phase(:,1)),19);
% AAA=zeros(1,18); count=zeros(1,18); 
% figure((phase_laminar_choose-1)*length(phase_select_laminar)+field_laminar_choose)
% subplot(4,1,4)
% for j=1:length(phase)
% for i=1:18
%    if phase(j,1)>=phase_sequence(i) && phase(j,1)<=phase_sequence(i+1)
%        count(i)= count(i)+1;
%        AAA(i)=AAA(i)+amplitude(j,2);
%        break
% end
% end
% end
% AAA=AAA./count;
% BBB=cat(2,AAA,AAA);
% bar((1:length(BBB))*20,BBB)  %文献上画了一个720度的图，学习一下
% % bar(BBB)
% title(strcat('横轴',num2str(hengzhou),'Hz  vs 纵轴',num2str(zongzhou),'Hz'))
% xticks([0 180 360 540 720])
% xticklabels([180:720])
H=0;
MI_total_total_light(trial_number,field_laminar_choose,field_laminar_choose,:,:)=MI_total;

    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                            No_light_sequence
MI_total_total_Nolight=zeros(size(Sorted_Wave_Epoch_NoLight,2),length(phase_select_laminar),length(field_select_laminar),low_phase_range(2),high_amplitude_range(2));
for trial_number=1:size(Sorted_Wave_Epoch_NoLight,2)
  if length(size(Sorted_Wave_Epoch_NoLight))==2
  q=squeeze(Sorted_Wave_Epoch_NoLight(:,trial_number));   
  q=q';
    else
 q=squeeze(Sorted_Wave_Epoch_NoLight(:,trial_number,:)); 
  end
 %%    调整排序，去除坏道
  if length(Bad_Ch)>0
 for i=1:length(Bad_Ch)
 temp(i)=find(LFP_array==Bad_Ch(i));
end
temp=sort(temp,'descend');
Bad_Ch=LFP_array(temp);
clear temp
 end
 if length(Bad_Ch)>0
     clear CH
     for i=1:length(Bad_Ch)
         CH(i)= find(LFP_array==Bad_Ch(i));  
     end
     q(CH,:)= [];
 end
 
 
 %% 上面都是常规处理数据
% q=q(:,40000:43500);
%% Generate the two signals that are interacting
% addpath('D:\My Documents\MATLAB\Work\lzj\function\cross-frequency')
% for phase_laminar_choose=1:length(phase_select_laminar)
for field_laminar_choose=1:length(field_select_laminar)
x_field=q(field_select_laminar(field_laminar_choose),:)';%amplitude的数据来源
x_phase=q(phase_select_laminar(field_laminar_choose),:)';%phase的数据来源
 for hengzhou=low_phase_range(1):low_phase_range(2)    %设定横轴频率范围
    for zongzhou=high_amplitude_range(1):high_amplitude_range(2)%   设定纵轴频率分析范围
    f1=hengzhou;
    f2=zongzhou;
[b, a, fCenterList, nTap, fEdges] = filterBank_prop(f1, f2, 1, 3, Sample_Fre, 'fir1', f2-f1);
% [x_filtered, x_analytic] = applyFilterBankThenHT_filtfilt(b, a, x_phase, nTap);
[~, x_analytic] = applyFilterBankThenHT_filtfilt(b, a, x_phase, nTap);
x_analytic=x_analytic(1:445,:);

%MIT程序包里带的两个程序，第一个没看明白啥意思，好像是选择带通滤波的相关参数

% x_analytic(2,:)=

% amplitude = abs(x_analytic);
phase = angle(x_analytic);
[~, x_analytic] = applyFilterBankThenHT_filtfilt(b, a, x_field, nTap);
x_analytic=x_analytic(1:445,:);
amplitude = abs(x_analytic);
%希尔伯特变换相关
phase_sequence=linspace(min(phase(:,1)),max(phase(:,1)),19);
AAA=zeros(1,18); count=zeros(1,18); 
%360度角分成18份

for j=1:length(phase)
for i=1:18
    if phase(j,1)>=phase_sequence(i) && phase(j,1)<=phase_sequence(i+1)
        count(i)= count(i)+1;
        AAA(i)=AAA(i)+amplitude(j,2);
        break
    end
end
end
AAA=AAA./count;%18个区间每个区间的amplitude平均值
% bar(AAA)
H=0;
for i=1:length(AAA)
   H=H+ AAA(i)/sum(AAA)*log(AAA(i)/sum(AAA));
end
H=-H;Hmax=log(length(AAA));
MI=(Hmax-H)/Hmax;%用这个评价两频率的相关性
MI_total(hengzhou,zongzhou)=MI;
    end
end

%%  后面画最强相关性点的例子各角度振幅平均值

% temp=max(MI_total);temp=find(temp==max(temp));
% temp2=max(MI_total');temp2=find(temp2==max(temp2));
% for hengzhou=temp2:temp2
%     for zongzhou=temp:temp
%     f1=hengzhou;
%     f2=zongzhou;
% [b, a, fCenterList, nTap, fEdges] = filterBank_prop(f1, f2, 1, 3, Sample_Fre, 'fir1', f2-f1);
% [~, x_analytic] = applyFilterBankThenHT_filtfilt(b, a, x_phase, nTap);
% 
% % x_analytic(2,:)=
% 
% % amplitude = abs(x_analytic);
% phase = angle(x_analytic);
% [~, x_analytic] = applyFilterBankThenHT_filtfilt(b, a, x_field, nTap);
% amplitude = abs(x_analytic);
% phase_sequence=linspace(min(phase(:,1)),max(phase(:,1)),19);
% AAA=zeros(1,18); count=zeros(1,18); 
% figure((phase_laminar_choose-1)*length(phase_select_laminar)+field_laminar_choose)
% subplot(4,1,4)
% for j=1:length(phase)
% for i=1:18
%    if phase(j,1)>=phase_sequence(i) && phase(j,1)<=phase_sequence(i+1)
%        count(i)= count(i)+1;
%        AAA(i)=AAA(i)+amplitude(j,2);
%        break
% end
% end
% end
% AAA=AAA./count;
% BBB=cat(2,AAA,AAA);
% bar((1:length(BBB))*20,BBB)  %文献上画了一个720度的图，学习一下
% % bar(BBB)
% title(strcat('横轴',num2str(hengzhou),'Hz  vs 纵轴',num2str(zongzhou),'Hz'))
% xticks([0 180 360 540 720])
% xticklabels([180:720])
H=0;
MI_total_total_Nolight(trial_number,field_laminar_choose,field_laminar_choose,:,:)=MI_total;

    end
end
% end
% field_select_laminar
% figure()
%%   绘制频率间相关性热图

for i=1:size(MI_total_total_Nolight,2)
figure()
temp_NoLight=squeeze(MI_total_total_Nolight(:,i,i,:,:));
temp_NoLight=squeeze(mean(temp_NoLight,1));
temp_Light=squeeze(MI_total_total_light(:,i,i,:,:));
temp_Light=squeeze(mean(temp_Light,1));
subplot(2,2,3)
 imagesc(temp_NoLight');colorbar
title(strcat(FileName,'振幅来源第',num2str(field_select_laminar(i)),'道相位来源第',num2str(phase_select_laminar(i)),'道分析时间'...
 ,num2str(Analysis_End_Time),'无光'))
set(gca,'clim',[0 0.15],'TickDir','out', 'YDir','normal')
% set(gca,'TickDir','out', 'YDir','normal')
subplot(2,2,1) 
imagesc(temp_Light');colorbar
title(strcat(FileName,'振幅来源第',num2str(field_select_laminar(i)),'道相位来源第',num2str(phase_select_laminar(i)),'道分析时间'...
 ,num2str(Analysis_End_Time),'有光'))
set(gca,'clim',[0 0.15],'TickDir','out', 'YDir','normal')
% figure()
% %%    一次trial的field图
% subplot(3,2,3)
% plot(squeeze(mean(Sorted_Wave_Epoch_NoLight(field_select_laminar(field_laminar_choose),2,:),2)))
% 
% title('无光例子  一次trial')
% subplot(3,2,4)
% plot(squeeze(mean(Sorted_Wave_Epoch_Light(field_select_laminar(field_laminar_choose),2,:),2)))
% 
% title('有光例子 一次trial')
% %%   频率图的平均图
% % figure(1)
% subplot(3,2,5)
% % Sorted_Wave_Epoch_NoLight_Temp=squeeze(mean(Sorted_Wave_Epoch_Light_Temp(field_select_laminar,:,:),2));
for j=1:size(Sorted_Wave_Epoch_NoLight,2)
    
   wave=squeeze(Sorted_Wave_Epoch_NoLight(field_select_laminar(i),j,:))';
[fft_consequence,fft_frequency]=func_fft_plot_frequency(wave,Sample_Fre,'r');
Sorted_Wave_Epoch_NoLight_fftconsequence(:,j)=fft_consequence;

end
for j=1:size(Sorted_Wave_Epoch_Light,2)    
wave=squeeze(Sorted_Wave_Epoch_Light(field_select_laminar(i),j,:))';
[fft_consequence,fft_frequency]=func_fft_plot_frequency(wave,Sample_Fre,'r');
Sorted_Wave_Epoch_Light_fftconsequence(:,j)=fft_consequence;
end
% figure()
subplot(2,2,[2,4])
plot(fft_frequency,mean(10*log10(Sorted_Wave_Epoch_Light_fftconsequence),2),'b')
hold on
plot(fft_frequency,mean(10*log10(Sorted_Wave_Epoch_NoLight_fftconsequence),2),'r')
xlim([0 90])
% axis([0 90 15 40]);hold on
title('蓝色有光，红色无光')
xlabel('频率/Hz');
 ylabel('振幅/分贝');
end