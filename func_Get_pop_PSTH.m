function [PSTH_G,TimeVal]=func_Get_pop_PSTH(raw,Pre_Stim_Time,Post_Stim_Time,BinSize,smooth_flag)
[Len,cc]=size(raw);
j=0;
%  BinSize=0.005;
%  smooth_flag=0;
% 
for i=1:Len
      Spk_File=raw{i,2};
      Spk_Ch=raw{i,4};
         if max(isnan(Spk_File))==0 && max(isnan(Spk_Ch))==0 
             try
             [PSTH,TimeVal]=func_Get_Tri_PSTH(Spk_File,Spk_Ch,Pre_Stim_Time,Post_Stim_Time,BinSize,smooth_flag);
             catch
             warning(strcat(Spk_File,'文件可能有问题，如有需要请联系李子杰'))    
             end
             j=j+1;
             try
             if size(PSTH,1)>size(PSTH_G,2) 
                 PSTH_G(:,size(PSTH_G,2)+1:size(PSTH,1),:)=NaN;
%              PSTH=PSTH(1:size(PSTH_G,2),:);
             end
            if size(PSTH,1)<size(PSTH_G,2)
                PSTH(size(PSTH,1)+1:size(PSTH_G,2),:)=NaN;
% %              PSTH_G=PSTH_G(:,1:size(PSTH,1),:);
             end  
             catch
             end
             PSTH_G(j,:,:)=PSTH;
         end
end
