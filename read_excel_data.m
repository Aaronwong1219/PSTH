clear all
close all
%%
ExcelFile='TRN_recording data.xlsx';
Start_line=2; 
End_line=3334;

Layer={'C'};% control:C,optogentic:O;
% Res_Type='On';
global Sound_Type 
Sound_Type='pure_tone250';
Record_site='TRN';
Sample_Fre_Tri=3051.76;
Durg=0;% 0=before,1=after

%% load the excel file
File_path='G:\Original_Data\Data_excel\';
% File_path='D:\Original_Data\Data_Tri\';
[num, txt, raw] = xlsread(strcat(File_path,ExcelFile));
[Len,cc]=size(raw);
% LFP_mean_without=[];
% LFP_mean_with=[];
% PSTH_mean_without=[];
% PSTH_mean_with=[];
LFP_size=[];
PSTH_size=[];
LFP={};
Select_ID=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%数据读取与处理%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for Laminar_Number=1:length(Layer);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% select a group of LFP to plot
% clearvars -except Select_ID Record_site  txt Start_line Sound_Type Laminar_Number Layer raw cc End_line ExcelFile File_path Len num Post_Stim_Time Pre_Stim_Time LFP_mean_without LFP_mean_with PSTH_mean_without LFP_G PSTH_mean_with1 PSTH_mean_with LFP_size PSTH_size LFP
clearvars -except  Sample_Fre_Tri Select_ID Record_site  txt Start_line Sound_Type Laminar_Number Layer raw cc End_line ExcelFile File_path Len num Post_Stim_Time Pre_Stim_Time LFP_mean_without1 LFP_mean_without LFP_mean_with PSTH_mean_without LFP_G PSTH_mean_with1 PSTH_mean_with LFP_size PSTH_size analysis_time1 analysis_time2 LFP template

% subplot(2,2,Laminar_Number)
j=0;
for ID=Start_line:End_line
       if strcmp(raw{ID,5},Layer{Laminar_Number})==1 && strcmp(raw{ID,7},Sound_Type)==1&& strcmp(raw{ID,6},Record_site)==1
        j=j+1;
        raw_LFP_Selected(j,1:cc)=raw(ID,1:cc);
    end
end
end

j=0;
for ID=Start_line:End_line
    if strcmp(raw{ID,5},Layer{Laminar_Number})==1 && strcmp(raw{ID,7},Sound_Type)==1&& strcmp(raw{ID,6},Record_site)==1
        j=j+1;
     
        raw_spike_Selected(j,1:cc)=raw(ID,1:cc);
    end
end

%%  截取并计算PSTH
Pre_Stim_Time = 1;
Post_Stim_Time = 5;
BinSize=0.005;
smooth_flag=0;
       
[PSTH_total,TimeVal]=func_Get_pop_PSTH(raw_spike_Selected,Pre_Stim_Time,Post_Stim_Time,BinSize,smooth_flag);       
 PSTH_mean = squeeze(mean(PSTH_total,2))  ;
        
 figure()
 imagesc(TimeVal,[0 size(PSTH_total,1)],PSTH_mean,'CDataMapping','scaled');
  xlim([0 2]);
 %%
 
 save(['G:\Original_Data\Data_TRN\',strcat(Record_site,'_',Layer{1},'_',Sound_Type),'.mat'],'raw','raw_spike_Selected','raw_LFP_Selected','Sample_Fre_Tri','PSTH_total','TimeVal','-v7.3');



