function rplot=FreelyMovingSpectrogramReport(data_ds,t_ds,adc_ds,t_adc_ds,scale,channels)
%clear all
%% load ethovision position data
SR=2500;
%channels=[9:13];
rplot=figure;
%filename = 'N8control_freely.csv';

ethopos=0; %use ethovision information

if(ethopos)
ethovis = csvread(filename);
ethovis_sr=1/ethovis(2,2);

%correct video timing with effective frame rate
% effectiveFrameRate=6.47; %lower=stretch
vid_delay=0;
t_vid=vid_delay+ethovis(:,1).*(1/(effectiveFrameRate/7));
%make sure: first or second video time?
end

%% 1) show PSD of channels and velocity+laser - precompute
clear TFDall
tmin=1; %index
tmax_s=478; %in seconds
%tmax_s=88

if(ethopos)
t_etho_max=find(t_vid>tmax_s,1,'first');
end
t_data_max=find(t_ds>tmax_s,1,'first'); 
t_adc_max=find(t_adc_ds>tmax_s,1,'first'); 

window=2*SR;
overlap=1.8*SR;
%overlap=1.5*SR;

%precompute
ch_cnt=1;
for i=channels
    i
    [~,F,T,P] = spectrogram(data_ds(i,tmin*SR:tmax_s*SR),window,overlap,2^16,SR);
    %truncate and store
    f_cut=find(F>20);
    F(f_cut)=[];
    P(f_cut,:)=[];
    TFDall(ch_cnt,:,:)=P;
    ch_cnt=ch_cnt+1;
end

%% plot it
if(ethopos)
plotshift=3;
subplot(length(channels)+plotshift,1,1);
plot(t_vid(tmin:t_etho_max),ethovis(tmin:t_etho_max,5)); %plot velocity
%plot(t_adc_ds(tmin:t_adc_max),adc_ds(1,tmin:t_adc_max));
ylabel('velocity')
xlim([0 tmax_s]);

subplot(length(channels)+plotshift,1,2);
plot(t_vid(tmin:t_etho_max),sqrt(ethovis(tmin:t_etho_max,3).^2+ethovis(tmin:t_etho_max,4).^2));
ylabel('position')
xlim([0 tmax_s]);
else
    plotshift=1;
end

subplot(length(channels)+plotshift,1,1);
plot(t_adc_ds(tmin:t_adc_max),adc_ds(tmin:t_adc_max));
ylabel('light')
xlim([0 tmax_s]);
cnt=2;
ch_cnt=1;
for i=channels
    subplot(length(channels)+plotshift,1,cnt);
    imagesc(T,F,squeeze(TFDall(ch_cnt,:,:)));
    colormap('jet');
    axis xy
    ylim([3 15])
    ylabel(['ch ',num2str(i)])
    xlim([0 tmax_s]);
    cax=caxis;
    caxis([0 cax(2)./scale])
    %caxis([0 100000])
    cnt=cnt+1;
    ch_cnt=ch_cnt+1;
end

%% 2) Theta power to position mapping
% go through time steps of ethovision and compute the corresponding theta
%power
if(false)
    normalize=1;
    if(normalize)
        I_norm=find(T>280 & T<290); % baseline range!
        P_norm2=mean(TFDall(:,:,I_norm),1);
        P_norm=mean(squeeze(P_norm2),2);
    end

    power=zeros(size(P,1),size(P,2));
    diffT=diff(T(1:2));

    %average ethovision data according to spectrogram discretization
    bin_x=zeros(1,length(T));
    bin_y=zeros(1,length(T));
    bin_v=zeros(1,length(T));
    bin_r=zeros(1,length(T));

    for i=1:length(T)
        %compute average / normalize
        if(normalize)
            power(:,i)=10*log10(mean(TFDall(:,:,i),1)'./P_norm);
        else
            power(:,i)=mean(TFDall(:,:,i),1);
        end

        %average x,y,velocity data to spectrogram bin size
        I_t = find(t_vid>=(T(i)-(diffT/2)) & t_vid<(T(i)+(diffT/2)));
        if ~isempty(I_t)
           bin_x(i)=mean(ethovis(I_t,3)); 
           bin_y(i)=mean(ethovis(I_t,4));
           bin_v(i)=mean(ethovis(I_t,8));
           bin_r(i)=sqrt(bin_x(i)^2+bin_y(i)^2);
        end
    end

    %% statistics over area for no movement

    I_NoMovement=(find(bin_v<4));

    I_NoMovement_final=[];
    cnt=1;
    for (i=2:length(I_NoMovement))
        if (I_NoMovement(i)== I_NoMovement(i-1)+1)
            cnt=cnt+1;
        else
            if (cnt>5)
                I_NoMovement_final(end+1:end+cnt) = I_NoMovement(i-cnt:i-1);
            end
            cnt=0;
        end

    end


    I_center_NoMovement_final=intersect(find(bin_r<5), I_NoMovement_final);
    I_middle_NoMovement_final=intersect(find(bin_r>5 & bin_r<7.5),I_NoMovement_final) ;
    I_outer_NoMovement_final=intersect(find(bin_r>7.5), I_NoMovement_final);

    mean_center_NoMovement_final=mean(power(:,I_center_NoMovement_final),2);
    mean_middle_NoMovement_final=mean(power(:,I_middle_NoMovement_final),2);
    mean_outer_NoMovement_final=mean(power(:,I_outer_NoMovement_final),2);

    std_center_NoMovement_final=std(power(:,I_center_NoMovement_final),[],2);
    std_middle_NoMovement_final=std(power(:,I_middle_NoMovement_final),[],2);
    std_outer_NoMovement_final=std(power(:,I_outer_NoMovement_final),[],2);

    h=plot(F,mean_center_NoMovement_final,'k');
    set(h,'LineWidth',2);

    hold on

    h=plot(F,mean_middle_NoMovement_final,'r');
    set(h,'LineWidth',2);

    h=plot(F,mean_outer_NoMovement_final,'g');
    set(h,'LineWidth',2);

    plot(F,mean_center-std_center_NoMovement_final,'k--')
    plot(F,mean_middle-std_middle_NoMovement_final,'r--')
    plot(F,mean_outer-std_outer_NoMovement_final,'g--')

    plot(F,mean_center+std_center_NoMovement_final,'g--')
    plot(F,mean_middle+std_middle_NoMovement_final,'r--')
    plot(F,mean_outer+std_outer_NoMovement_final,'k--')
    xlim([0 14])

    %%
    %%% statistics over area for Movement 7-15 cm/s
    I_Movement=(find(bin_v>6 & bin_v<15));

    I_Movement_final=[];
    cnt=1;
    for (i=2:length(I_Movement))
        if (I_Movement(i)== I_Movement(i-1)+1)
            cnt=cnt+1;
        else
            if (cnt>5)
                I_Movement_final(end+1:end+cnt) = I_Movement(i-cnt:i-1);
            end
            cnt=0;
        end

    end


    I_center_Movement_final=intersect(find(bin_r<5), I_Movement_final);
    I_middle_Movement_final=intersect(find(bin_r>5 & bin_r<7.5),I_Movement_final) ;
    I_outer_Movement_final=intersect(find(bin_r>7.5), I_Movement_final);

    mean_center_Movement_final=mean(power(:,I_center_Movement_final),2);
    mean_middle_Movement_final=mean(power(:,I_middle_Movement_final),2);
    mean_outer_Movement_final=mean(power(:,I_outer_Movement_final),2);

    std_center_Movement_final=std(power(:,I_center_Movement_final),[],2);
    std_middle_Movement_final=std(power(:,I_middle_Movement_final),[],2);
    std_outer_Movement_final=std(power(:,I_outer_Movement_final),[],2);

    h=plot(F,mean_center_Movement_final,'k');
    set(h,'LineWidth',2);

    hold on

    h=plot(F,mean_middle_Movement_final,'r');
    set(h,'LineWidth',2);

    h=plot(F,mean_outer_Movement_final,'g');
    set(h,'LineWidth',2);

    plot(F,mean_center-std_center_Movement_final,'k--')
    plot(F,mean_middle-std_middle_Movement_final,'r--')
    plot(F,mean_outer-std_outer_Movement_final,'g--')

    plot(F,mean_center+std_center_Movement_final,'g--')
    plot(F,mean_middle+std_middle_Movement_final,'r--')
    plot(F,mean_outer+std_outer_Movement_final,'k--')

    xlim([0 14])
end
