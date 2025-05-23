function [LFPs,time3]=func_Get_pop_LFP(raw,Pre_Stim_Time,Post_Stim_Time)
[Len,cc]=size(raw);
j=0;
    for i=1:Len
        LFP_File=raw{i,2};
        LFP_Ch=raw{i,3};
        if i==1
            j=j+1;
           [LFP,time3]=func_Get_Tri_LFP(LFP_File,LFP_Ch,Pre_Stim_Time,Post_Stim_Time);
            LFPs(j,:,:)=LFP;
            clear LFP;
        end
        if i>1&&max(isnan(LFP_File))==0 && max(isnan(LFP_Ch))==0 && (strcmp(LFP_File,raw{i-1,2})==0 || raw{i,3}~= raw{i-1,3})
            % 增加了一个逻辑判断 文件名与道数均与上一个文件相同时才会跳过，以适应阵列电极数据
            try
             [LFP,time3]=func_Get_Tri_LFP(LFP_File,LFP_Ch,Pre_Stim_Time,Post_Stim_Time);
            catch
               warning(strcat(LFP_File,'文件可能有问题，如有需要请联系李子杰')) 
            end
             j=j+1;
              if size(LFP,1)>size(LFPs,2)
                  LFPs(:,size(LFPs,2)+1:size(LFP,1),:)=NaN;
%                   LFP=LFP(1:size(LFPs,2),:);
              end
              if size(LFP,1)<size(LFPs,2)
                  LFP(size(LFP,1)+1:size(LFPs,2),:)=NaN;
%                   LFPs=LFPs(:,1:size(LFP,1),:);
              end
             LFPs(j,:,:)=LFP;
       
        clear LFP;
        end
    end
