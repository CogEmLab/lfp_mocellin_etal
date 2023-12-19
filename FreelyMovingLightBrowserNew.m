clear all
load 51040ds
SR_adc=SR;
channel=25;

%% options
window=3; 
overlap=2.6;
tmin=1; %index
tmax=floor(t_ds(end)); %in seconds
growby=1; %wider laser regions

threshold=2.5; 
baseline=0; %which type of normalization (0/1/2)?
min_base=110 ;%for baseline normalization
max_base=120; 

plotit=1;
MINFREQ=4;
MAXFREQ=14;

%more ideas: statistics / epoch

%% fft for signal
if(adc_ds~= 0)
[~,F_adc,T_adc,P_adc] = spectrogram(adc_ds,window*SR_adc,overlap*SR_adc,2^16,SR_adc);
mean_F=mean(P_adc(F_adc>0.1 & F_adc<0.2,:),1);
T_laser=find(mean_F>threshold);
    subplot(2,1,2);
    plot(t_adc_ds,adc_ds);
    hold on
    plot(T_adc,mean_F,'r')
    xlim([tmin tmax])
else
    T_laser=0;
end
[~,F,T,P] = spectrogram(data_ds(channel,:),window*SR,overlap*SR,2^16,SR);
if(adc_ds~= 0)
if(~isequal(T,T_adc))
    error('FFT time vectors not equal');
end
end

    subplot(2,1,1);
    imagesc(T,F,P);
     cax=caxis;
    caxis([0 cax(2)./10])
    colormap('jet');
    axis xy
    ylim([0 20])
    xlim([tmin tmax])
    set(gcf, 'Position', [100, 100, 1200, 1000]);
    a=gca;
    
    
if(growby>0)
    T_laser_grown=[];
    for i=T_laser
        for j=growby:-1:1
            %minus j
%             if(~any(T_laser_grown==(i-j)))
%                 T_laser_grown(end+1)=i-j;
%             end
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


while (1)
    rect = getrect(a)

    tstart=rect(1);
    tend=rect(1)+rect(3);

    
    figure
    subplot(2,1,1)
    
    %complete PSD
    try
    [Pw, Fw]=pwelch(data_ds(channel,t_ds>tstart & t_ds<tend),3000,2600,2^17,SR);
    end
    plot(Fw,Pw);
    xlim([MINFREQ MAXFREQ])
    title(['T: ',num2str(tstart) '-' num2str(tend)]);
    
    subplot(2,1,2)
    %light / nolight psd
    
    
    power_laser=[];
    power_nolaser=[];
    cnt_laser=0;
    cnt_nolaser=0;

    %find theta peak and peak frequency for those
    
    for(i=1:length(T))
        if(T(i)>tstart & T(i)<tend)
            if(any(T_laser==i))
                if(baseline==1)
                    power_laser(end+1,:)=10*log10((P(I_f,i))./normalization);
                elseif(baseline==2)
                    power_laser(end+1,:)=(P(I_f,i))./normalization;
                else %baseline=0
                    power_laser(end+1,:)=P(I_f,i);
                end
                cnt_laser=cnt_laser+1;
            else
                if(baseline==1)
                    power_nolaser(end+1,:)=10*log10((P(I_f,i))./normalization);
                elseif(baseline==2)
                    power_nolaser(end+1,:)=(P(I_f,i))./normalization;
                else %baseline=0
                    power_nolaser(end+1,:)=P(I_f,i);
                end
                cnt_nolaser=cnt_nolaser+1;
            end
        end
    end
    
    p=plot(F(I_f),mean(power_laser),'r');
    set(p,'LineWidth',2);

    hold on

    p=plot(F(I_f),mean(power_nolaser),'k');
    set(p,'LineWidth',2);

    plot(F(I_f),mean(power_laser)+(std(power_laser)),'r--')
    plot(F(I_f),mean(power_nolaser)+(std(power_nolaser)),'k--')

    xlim([MINFREQ MAXFREQ])
    legend('laser','no laser')
    title(['laser=' num2str(cnt_laser) ', nolaser=' num2str(cnt_nolaser)])
end



