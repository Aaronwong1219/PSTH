function [template_G]=func_Get_pop_template(raw_Selected,Pre_Stim_Time,Post_Stim_Time,BinSize,smooth_flag)
[Len,cc]=size(raw_Selected);
j=0;
%  BinSize=0.005;
%  smooth_flag=0;
% 
for i=1:Len
      Spk_File=raw_Selected{i,2};
      Spk_Ch=raw_Selected{i,4};
         if max(isnan(Spk_File))==0 && max(isnan(Spk_Ch))==0 
             try
             File_path='F:\Original_Data\Data_Tri\';
             load(cat(2,File_path,Spk_File,'tri.mat'));
             catch
             File_path='G:\Original_Data\Data_Tri\';
             load(cat(2,File_path,Spk_File,'tri.mat'));   
             end
             template=SpikeWaveform.(strcat('wave',num2str(Spk_Ch)));
             j=j+1;
                template_G(j,:)=template;
         end
end
