%%
load('G:\Original_Data\Data_PSGI\MRN.mat');

PSTH_mean = mean(PSTH,1);
figure(1)
plot(PSTH','b'); hold on
plot(PSTH_mean,'r','linewidth',3)


for i=1:size(PSTH,1)
    max_val=max(PSTH(i,:));
    PSTH_normolized(i,:)=PSTH(i,:)./max_val;
end
%%
 idx = sum(PSTH_normolized(:,760:1000),2) ;
%  idx = sum(PSTH_normolized(:,760:800),2) ;
 [PSTH_normolized_sorted,ID]=sort(idx,'descend');%原始数据按照过阈值的曲线下面积进行排序
normalized_PSTH_data=PSTH_normolized(ID,:);
% normalized_PSTH_data(1:30,:) = normalized_PSTH_data(1:30,:)+0.1
 figure()
 imagesc([1 size(normalized_PSTH_data,2)],[1 size(normalized_PSTH_data,1)],normalized_PSTH_data,'CDataMapping','scaled');


%%
figure()
for i=1:size(template_G,1)
    plot(template_G(i,:)); hold on
end

figure()
spike_sample=25000;
template_GG=[];
for i=1:size(template_G,1)
    template_GG(i,:)=resample(template_G(i,:),30,1);
end

[Fs_sequence]=func_plot_spikedot_figure(template_GG,spike_sample,6,7);
% 将RS 和FS 神经元波形分开画
figure()
for j=1:size(template_G,1)
    if Fs_sequence(j)== 0;
        plot(template_G(j,:),'r'); hold on
    else
        plot(template_G(j,:),'b'); hold on
    end
end




