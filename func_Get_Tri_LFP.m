

function [LFPs,time3]=func_Get_Tri_LFP(LFP_File,Signal_Ch_Original,Pre_Stim_Time,Post_Stim_Time)


% LFP_File='QLM01_3202';
% Signal_Ch_Original=7;

% Pre_Stim_Time=1; %最好不要改
% Post_Stim_Time=6; %最好不要改

try
    File_path='G:\Original_Data\Data_Tri\';
    load(cat(2,File_path,LFP_File,'tri.mat'));
catch
   File_path='G:\Original_Data\Data_Tri\';
    load(cat(2,File_path,LFP_File,'tri.mat'));
end
try
    Sample_Fre=Sample_Fre_Tri;
catch
    Sample_Fre=1000;   
end
try
LFP=func_cut_from_continous_data(Wave_Epoch_continue,Pre_Stim_Time,Post_Stim_Time,Stimuli_Time,Sample_Fre);
catch
   Stimuli_Time=Stimuli_Time(1:end-1);
LFP=func_cut_from_continous_data(Wave_Epoch_continue,Pre_Stim_Time,Post_Stim_Time,Stimuli_Time,Sample_Fre);   
end
    
if Sample_Fre==1000   
    LFP=10*LFP;
end
Signal_Ch=find(Signal_Ch_Original==LFP_array);

%% sort the LFP by light or without light 
global Sound_Type
 if strcmp(Sound_Type,'click-train-4vs40Hz')
load('D:\MyDocuments\Work\tri_files\Two_sti_rand_2.mat')
light_click_train_control=Two_sti_rand(:,2);   
 else
load('D:\Work\light_click_train_control.mat');  
 end   
light_click_train_control=cat(1,light_click_train_control,light_click_train_control);
light_click_train_control=light_click_train_control(1:size(LFP,2));
if strcmp(Sound_Type,'click-train-4vs40Hz')
    IX=1:length(light_click_train_control);
else
    [M,IX]=sort(light_click_train_control);
%     IX=1:length(light_click_train_control);

end
light_point=length(find(light_click_train_control==0));
for i=1:length(IX)
    temp(:,i,:)=LFP(:,IX(i),:);
end
if length(size(temp))==2
 tempnew(1,:,:)=temp';
 temp=tempnew;
end
clear LFP;LFP=temp;clear temp;

%% plot the mean LFP together
time3=(-1*Pre_Stim_Time:1/Sample_Fre:Post_Stim_Time-1/Sample_Fre);
if length(LFP)<length(time3)
    time3=time3(1:length(LFP));
end
LFPs=squeeze(LFP(Signal_Ch,:,:));


 


