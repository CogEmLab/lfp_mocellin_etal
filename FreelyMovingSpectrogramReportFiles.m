function rplot=FreelyMovingSpectrogramReportFiles(data,channels,SR)

%channels=[9:13];
rplot=figure;
%filename = 'N8control_freely.csv';
tend=floor(length(data)/SR);

%% 1) show PSD of channels and velocity+laser - precompute
window=2*SR;
overlap=1.8*SR;
%overlap=1.5*SR;

%precompute
ch_cnt=1;
for i=channels
    i
    [~,F,T,P] = spectrogram(data(i,:),window,overlap,2^17,SR);
    %truncate and store
    f_cut=find(F>20);
    F(f_cut)=[];
    P(f_cut,:)=[];
    TFDall(ch_cnt,:,:)=P;
    ch_cnt=ch_cnt+1;
end

%% plot it

plotshift=1;
subplot(length(channels)+plotshift,1,1);
%plot(t_adc_ds(tmin:t_adc_max),adc_ds(tmin:t_adc_max));
ds_factor=250;
light_ds=downsample(data(33,:),ds_factor);
time_ds=0:1/(SR/ds_factor):tend-1/(SR/ds_factor);
bar(time_ds,light_ds);
ylabel('light')
xlim([1 tend-1]);
cnt=2;
ch_cnt=1;
for i=channels
    subplot(length(channels)+plotshift,1,cnt);
    imagesc(T,F,squeeze(TFDall(ch_cnt,:,:)));
    colormap('jet');
    axis xy
    ylim([3 15])
    ylabel(['ch ',num2str(i)])
    xlim([1 tend-1]);
    cax=caxis;
    caxis([0 cax(2)./10])
    %caxis([0 100000])
    cnt=cnt+1;
    ch_cnt=ch_cnt+1;
end