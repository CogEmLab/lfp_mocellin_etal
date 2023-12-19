function [rplot]=FreelyMovingPSDReportNew(data_ds,t_ds,adc_ds,t_adc_ds,scale,channels,SR)
%function [rplot,frequencies,power_laser,power_nolaser]=FreelyMovingPSDReport(file,channels)
%clear all
SR_adc=SR;
%channels=[15:16];

%% options
window=3; 
overlap=2.6;
tmin=0; %index
tmax=floor(t_ds(end)); %in seconds
growby=1; %wider laser regions

threshold=2.5; 
baseline=0; %which type of normalization (0/1/2)?
min_base=110 ;%for baseline normalization
max_base=120; 

plotit=0;
rplot=figure;
MINFREQ=3;
MAXFREQ=14;
%SR = 2500;

%more ideas: statistics / epoch
power_laser_all=[];
power_nolaser_all=[];

%% fft for adc
[~,F_adc,T_adc,P_adc] = spectrogram(adc_ds,window*SR,overlap*SR,2^16,SR);
mean_F=mean(P_adc(F_adc>0.1 & F_adc<0.2,:),1);
T_laser=find(mean_F>threshold);
if(plotit)
    subplot(2,1,2);
    plot(t_adc_ds,adc_ds);
    hold on
    plot(T_adc,mean_F,'r')
    xlim([tmin tmax])
end

%% grow lasers regions
if(growby>0)
    T_laser_grown=[];
    for i=T_laser
        for j=growby:-1:1
            %minus j
            if(~any(T_laser_grown==(i-j)))
                T_laser_grown(end+1)=i-j;
            end
            % i
            if(~any(T_laser_grown==(i)))
                T_laser_grown(end+1)=i;
            end        
            %plus j
            if(~any(T_laser_grown==(i+j)))
                T_laser_grown(end+1)=i+j;
            end    
        end 
    end
    T_laser=T_laser_grown;
end

%% fft for signal

plcnt=1;
for channel=channels
channel
[~,F,T,P] = spectrogram(data_ds(channel,:),window*SR,overlap*SR,2^16,SR);
if(~isequal(T,T_adc))
    error('FFT time vectors not equal');
end
if(plotit)
    subplot(2,1,1);
    imagesc(T,F,P);
     cax=caxis;
    caxis([0 cax(2)./10])
    colormap('jet');
    axis xy
    ylim([0 20])
    xlim([tmin tmax])
end

%% get average PSD for laser / no laser

I_f = find(F>MINFREQ & F<MAXFREQ);

if(baseline) %as in book
    t_base=find(T>min_base & T<max_base); 
    freqs_baseline=P(I_f,t_base);
    normalization=mean(freqs_baseline,2);
    if(baseline==2) %as by adriano
        freqs=find(F<MAXFREQ);
        freqs_baseline=P(freqs,t_base);
        allf=mean(freqs_baseline,2);
        normalization=sum(allf);
    end
end

power_laser=[];
power_nolaser=[];

%find theta peak and peak frequency for those

for(i=1:length(T))
    if(T(i)>tmin)
        if(any(T_laser==i))
            if(baseline==1)
                power_laser(end+1,:)=10*log10((P(I_f,i))./normalization);
            elseif(baseline==2)
                power_laser(end+1,:)=(P(I_f,i))./normalization;
            else %baseline=0
                power_laser(end+1,:)=P(I_f,i);
            end
            if(plotit)
               scatter(T(i),1,[],'gx'); hold on
            end
        else
            if(baseline==1)
                power_nolaser(end+1,:)=10*log10((P(I_f,i))./normalization);
            elseif(baseline==2)
                power_nolaser(end+1,:)=(P(I_f,i))./normalization;
            else %baseline=0
                power_nolaser(end+1,:)=P(I_f,i);
            end
%             if(plotit)
%                scatter(T(i),1,[],'rx'); hold on
%             end
        end
    end
end

disp(['Laser: ', num2str(length(T_laser))])
disp(['No Laser: ', num2str(length(power_nolaser))])

%% plot 
if(true)
%figure
subplot(1,length(channels),plcnt)
h=plot(F(I_f),mean(power_laser),'r');
set(h,'LineWidth',2);

hold on

h=plot(F(I_f),mean(power_nolaser),'k');
set(h,'LineWidth',2);

plot(F(I_f),mean(power_laser)+(std(power_laser)),'r--')
%plot(F(I_f),mean(power_laser)-(std(power_laser)),'k--')
plot(F(I_f),mean(power_nolaser)+(std(power_nolaser)),'k--')
%plot(F(I_f),mean(power_nolaser)-(std(power_nolaser)),'r--')

xlim([4 12])
xlabel('frequency')
ylabel('mean power \pm std')
%title(['window=',num2str(window),', overlap=',num2str(overlap),', growby=',num2str(growby),', baseline=',num2str(baseline)]);
title(['channel ',num2str(channel)]);
legend('laser','no laser')

frequencies=F(I_f);

if(true) %fig5 wall PSD using immobile/51040ds (channel 24)
power_laser=power_laser./800000;
power_nolaser=power_nolaser./800000; 
NiceFrequencyPlot2(power_laser,power_nolaser,F,MINFREQ,MAXFREQ,'Wall','Arch','Control',10)
NiceErrorBars(power_laser,power_nolaser,F,5,8,'6-8 Hz','Control','Arch',10,6,4156); %4153
NiceErrorBars(power_laser,power_nolaser,F,8,10,'8-10 Hz','Control','Arch',10,6,4526); %4526
end

if(false)

[l_peaks,l_peaks_I]=max(power_laser,[],2);
[nl_peaks,nl_peaks_I]=max(power_nolaser,[],2);

[l_m,l_m_i]=max(l_peaks);
[nl_m,nl_m_i]=max(nl_peaks);
laser_peak_amplitude=l_m
nolaser_peak_amplitude=nl_m
laser_peak_frequency=F(I_f(l_m_i))
nolaser_peak_frequency=F(I_f(nl_m_i))

Label{1}='laser';
Label{2}='no laser';

subplot(2,2,3)
boxplot([l_peaks;nl_peaks],[zeros(1,length(l_peaks)),ones(1,length(nl_peaks))]);
title('peak theta amp')
set(gca,'xtick',[1,2],'xticklabel',Label)

subplot(2,2,4)
boxplot([F(I_f(l_peaks_I));F(I_f(nl_peaks_I))],[zeros(1,length(l_peaks)),ones(1,length(nl_peaks))]);
title('peak theta freq')
set(gca,'xtick',[1,2],'xticklabel',Label)
end
end

power_laser_all(plcnt,:)=mean(power_laser,1);
power_nolaser_all(plcnt,:)=mean(power_nolaser,1);
plcnt=plcnt+1;
end

if(false)   
NiceFrequencyPlot2(power_laser_all,power_nolaser_all,F,MINFREQ,MAXFREQ,'Wall','Control','Arch',20)
%NiceErrorBars(power_laser,power_nolaser,F,6,8,'Wall 6-8 Hz','Control','Arch',6)
%NiceErrorBars(power_laser,power_nolaser,F,8,10,'Wall 8-10 Hz','Control','Arch',6)
end
