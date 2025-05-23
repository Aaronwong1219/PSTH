 function [Fs_sequence]=func_plot_spikedot_figure(template,spike_sample,varargin)
 %varargin 第一个变量是x_Threshold 第二个变量是y_Threshold
      %%      根据每一个细胞的第一峰半峰时间X轴，谷峰时间y轴进行绘图
    global qw
      if length(varargin)==2
          x_Threshold=varargin{1};
          y_Threshold=varargin{2};
      else
          x_Threshold=0.15;
          y_Threshold=0.25;
      end
      spike_sample=spike_sample/1000;%换算成1ms有几个采样点
      
    for i=1:size(template,1)
    temp=template(i,:);
    min_location=find(temp==min(temp));
    min_value=min(temp);
%    max_location=find(temp(min_location:end)==max(temp(min_location:end)))+min_location-1;
%     max_value=max(temp(min_location:end));
    %%
     max_location=find(temp==max(temp));
     max_value=max(temp);
%   max_location=find(temp(min_location:end)==max(temp(min_location:end)))+min_location-1;
%     max_value=max(temp(min_location:end));
    y_value=abs(max_location-min_location)/spike_sample;
    %%
    for differ_sequence=2:30
   differ_value=temp(differ_sequence)-temp(differ_sequence-1);
   if differ_value<0
       break
   end
    end
  qqw=(min_value-temp(differ_sequence))/2; %高度的一半
    qqw=temp(differ_sequence)+qqw;
  x_value=  length(find(temp<=qqw))/spike_sample;
 
%    sequence= find(qqw<0);
%    x1=sequence(1)-1;x2=sequence(1);y1=temp(sequence(1)-1);y2=temp(sequence(1));
%    y=min_value/2-temp(1);
%    k=(y2-y1)/(x2-x1);
%    b=(y1*x2-y1*x1-y2*x1+y1*x1)/(x2-x1);
%    x_before=(y-b)/k;
%    %%
%    clear x1 x2 y1 y2 y x k b
%     x1=sequence(end);x2=sequence(end)+1;y1=temp(sequence(end));y2=temp(sequence(end)+1);
%    y=min_value/2-temp(1);
%    k=(y2-y1)/(x2-x1);
%    b=(y1*x2-y1*x1-y2*x1+y1*x1)/(x2-x1);
%    x_after=(y-b)/k;
%    x_value=(x_after-x_before)/24.4414;
   %%

   hold on
  
%    text(x_value,y_value,num2str(i))
   clear temp
  
   if x_value<=x_Threshold && y_value<=y_Threshold
       Fs_sequence(i)=1;
        plot(x_value,y_value,'.b','markersize',10)
%         text(x_value,y_value,num2str(i))
       qw(i,1)=x_value;
        qw(i,2)=y_value;
   else
        plot(x_value,y_value,'.r','markersize',10)
       Fs_sequence(i)=0;  
         qw(i,1)=x_value;
        qw(i,2)=y_value(1);
   end
    end
    axis([2 15 2 25])
    xlabel('半峰宽度 （ms）')
    ylabel('峰谷间距 （ms）')
    
%     qw=0;
 end   %函数结尾