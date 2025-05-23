
% close all
% clear all
function  [PSTHs,TimeVal]=func_Get_Tri_PSTH(Spk_File,Signal_Ch_Original,Pre_Stim_Time,Post_Stim_Time,BinSize,smooth_flag)

% Pre_Stim_Time=1; %最好不要改
% Post_Stim_Time=6; %最好不要改
% BinSize=0.005;
% smooth_flag=0;

try
    File_path='F:\Original_Data\Data_Tri\';
%     File_path='E:\Data_Tri\';
    load(cat(2,File_path,Spk_File,'tri.mat'));
catch
    File_path='F:\Original_Data\Data_Tri\';
    load(cat(2,File_path,Spk_File,'tri.mat'));
end


Stim_Count=length(Stimuli_Time);



if length(Spike_Continous)>0
    name=fieldnames(Spike_Continous);
    Spike_Times= struct2cell(Spike_Continous);

    for i=1:length(name)
        Spk_Ch{i}=name{i}(5:end);
    end
end

for i=1:length(Spk_Ch)
    if strcmp(Spk_Ch{i},Signal_Ch_Original)==1
        Spike_Selected=i;
    end;
end

%% sort the LFP by light or without light 
global Sound_Type
 if strcmp(Sound_Type,'click-train-4vs40Hz')
load('Two_sti_rand_2.mat')
light_click_train_control=Two_sti_rand(:,2);   
 else
load('D:\Work\light_click_train_control.mat');  
 end  
% load('D:\MyDocuments\Work\MatData\light_click_train_control.mat');
light_click_train_control=light_click_train_control(1:Stim_Count);
if strcmp(Sound_Type,'click-train-4vs40Hz')
    IX=1:length(light_click_train_control);
else
    [M,IX]=sort(light_click_train_control);
end
% [M,IX]=sort(light_click_train_control);

%% read and plot spike data
if isempty(Spike_Selected)==0
   
    SpikeData=Spike_Times{Spike_Selected};
    %去除5ms以内放电
            for qqw=length(SpikeData):-1:2
                if SpikeData(qqw)-SpikeData(qqw-1)<=0.005
                    SpikeData(qqw)=[];
                end
            end      
            for trial=1:Stim_Count
                try
               Index1=max(find(SpikeData<=Stimuli_Time(trial)-Pre_Stim_Time));
               Index2=max(find(SpikeData<=Stimuli_Time(trial)+Post_Stim_Time));
               Neural_Data(find(IX==trial)).SpikeTime=SpikeData(Index1:Index2)-Stimuli_Time(trial);
                catch
                    break  %防止末尾index2超出记录所有的动作电位数量？不确定
                end
            end
            %%  PSTH部分
load('D:\Work\light_click_train_control.mat');
% load('D:\MyDocuments\Work\MatData\light_click_train_control.mat');
            [M,IX]=sort(light_click_train_control);
            light_point=find(M==1,1);
            [PSTHs,TimeVal]=func_Cal_PSTH_Click(Neural_Data,BinSize,smooth_flag,-1*Pre_Stim_Time,Post_Stim_Time);
    end




 


