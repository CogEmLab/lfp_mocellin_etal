
%%
clear all
clc
SR=1000;
%%
clear all
[t,amps,data1,aux2] = read_intan_data_leao('pinch_lightON_170531_143030_170531_150631.int');
%%
% [t2,amps2,data2,aux2] = read_intan_data_leao('hem2_2hz.int');
[t3,amps3,data3,aux3] = read_intan_data_leao('hem2_4hz.int');
[t4,amps4,data4,aux4] = read_intan_data_leao('hem2_8hz.int');
[t5,amps5,data5,aux5] = read_intan_data_leao('hem2_16hz.int');

%%
[P1 F1]=pwelch(detrend(data_1000(8,1:20000)),1000, [],2^18, 1000);
plot(F1,P1,'k');
hold on
% [P1 F1]=pwelch(detrend(data_1000(3,5000:10000)),1000, [],2^18, 1000);
% plot(F1,P1,'m');
[P2 F2]=pwelch(detrend(data_1000(8,20000:40000)),1000, [],2^18, 1000);
plot(F2,P2,'r');
[P3 F3]=pwelch(detrend(data_1000(8,40000:60000)),1000, [],2^18, 1000);
plot(F3,P3,'g');
xlim([0 20])

%%
data1=data1(:,1:25:end);

%%
[P1 F1]=pwelch(detrend(data1(16,1:20000)),1000, [],2^18, 1000);
plot(F1,P1,'k');
xlim([0 20])

hold on
[P1 F1]=pwelch(detrend(data1(16,20000:40000)),1000, [],2^18, 1000);
plot(F1,P1,'r');

[P1 F1]=pwelch(detrend(data1(16,40000:60000)),1000, [],2^18, 1000);
plot(F1,P1,'g');

xlim([0 10])
hold on
%%
[P1 F1]=pwelch(detrend(amplifier_data(22,:)),1000, [],2^18, 1000);
plot(F1,P1,'m');
[P1 F1]=pwelch(detrend(amplifier_data(24,:)),1000, [],2^18, 1000);
plot(F1,P1,'r');
xlim([0 20])
%%
[P1 F1]=pwelch(detrend(data_1000(4,1:15000)),15000, [],2^18, 1000);
plot(F1,P1,'b');
hold on
[P1 F1]=pwelch(detrend(data_1000(4,15000:35000)),20000, [],2^18, 1000);
plot(F1,P1,'r');
[P1 F1]=pwelch(detrend(data_1000(4,35000:50000)),15000, [],2^18, 1000);
plot(F1,P1,'g');
xlim([0 120]);
%%
[P2 F2]=pwelch(detrend(data_ds(9,20000:40000)),1000, [],2^16, 1000);
plot(F2,P2,'k');
hold on
[P3 F3]=pwelch(detrend(data_ds(9,40000:60000)),1000, [],2^16, 1000);
plot(F3,P3,'r');
[P4 F4]=pwelch(detrend(data_ds(9,130000:140000)),1000, [],2^18, 1000);
plot(F4,P4,'g');
xlim([3 20])


%%
[P3 F3]=pwelch(detrend(data(12,100*25000:25:130*25000)),1000, [],2^18, 1000);
plot(F3,P3,'g');
% [P4 F4]=pwelch(detrend(data(7,150*25000:25:180*25000)),1000, [],2^18, 1000);
% plot(F4,P4,'m');

xlim([0 20]);

%%
[P3 F3]=pwelch(detrend(data3(7,150000:25:200000)),[], [],2^12, 10000);
plot(F3,P3,'r');
[P4 F4]=pwelch(detrend(data4(15,150000:25:200000)),[], [],2^12, SR);
plot(F4,P4,'g');
[P5 F5]=pwelch(detrend(data5(15,150000:25:200000)),[], [],2^12, SR);
plot(F5,P5,'m');
xlim([0 50]);
hold off


%%

%%
for ch = 1:16
    SR=1000;
    subplot(4,4,ch);
[P F]=pwelch(detrend(data3(ch,1:25:50000)),[], [], 2^12, SR);
plot(F, P);
hold on
[P1 F1]=pwelch(detrend(data3(ch,150000:25:200000)),[], [],2^12, SR);

plot(F1,P1,'r');

% hold on
% [P2 F2]=pwelch(detrend(data3(ch,1:25:300000)),[], [],2^13, SR);
% plot(F2, P2, 'g');

xlim([0 70]);
hold off
end

%%
%getting the time course of theta power

theta=find(F>4 & F<12);

thetapower = mean(P(theta,:));
subplot(3,1,3)
plot(T,thetapower,'k-sq','markerfacecolor','y')
xlabel('Time (secs)')
ylabel('Theta Power')

%%
%%
[peakvalue peakindex]=max(P)

subplot(111)
plot(F,P,'k.-')
%%
% 
% <html>
% <table border=1><tr><td>one</td><td>two</td></tr></table>
% </html>
% 
hold on
plot(F(peakindex),peakvalue,'ro','markerfacecolor','r','markersize',5)

hold off

xlabel('Freq(Hz)')
ylabel('Power')
xlim([0 20])

peakfreq = F(peakindex)
text(10,0.2,['Peak Freq = ' num2str(peakfreq) 'Hz'])
%%
%Discrete TFD
clf


srate=2500;

window=2*srate;
overlap=1.8*srate;

%N1 cutted together
%begin=data_ds(18,147370:end);
%data=[begin data_ds(18,:)];

%N9 cutted together
%ch=11;
%begin=data_ds(ch,90*srate:140*srate);
%data=[begin data_ds(ch,:)];
%data(442106:end)=[];

%movement
%begin=adc_ds(2,90*srate:140*srate);
%data2=[begin adc_ds(2,:)];
%data2(442106:end)=[];

subplot(2,1,1)
[S F T P] = spectrogram(data_ds(),window,overlap,2^16,srate);
%subplot(3,1,[1 2])
imagesc(T,F,P)
colormap('jet');

ylim ([4 14])
axis xy
    cax=caxis;
    caxis([0 cax(2)./50])
    
 colorbar

ylabel('Frequency (Hz)')
xlabel('Time (secs)')

subplot(2,1,2);
plot(data2)

%%
srate=2500;
window=2*srate;
overlap=1.9*srate;

%channels=[9:16]
%channels=[17:26]
channels=[16:19]

cnt=1;

for ch = channels
    ch
    subplot(length(channels),1,cnt);    

    [S F T P] = spectrogram(double(data_ds(ch,:)),window,overlap,2^16,srate);

    f_cut=find(F>20);
    F(f_cut)=[];
    P(f_cut,:)=[];
    
    %if(cnt>2)
    %   P=P.*0.05;        
    %end
    
    TFDall(cnt,:,:)=P;
    
    imagesc(T,F,P)
    colormap('jet');
    axis xy
    
    ylim([3 14])
    ylabel(['ch ',num2str(ch)])
    
    cax=caxis;
    if(cnt<4)
        caxis([0 cax(2)./10])
    else
        caxis([0 cax(2)./10])
    end
    
    cnt=cnt+1;
end

%% average
figure
subplot(4,1,3:4)
imagesc(T,F,squeeze(mean(TFDall,1)));
colormap('jet');
axis xy

ylim([3 15])
ylabel('Frequency (Hz)')
xlabel('Time (s)')
cax=caxis;
%caxis([0 cax(2)./2.5])
caxis([0 cax(2)./10])
ylim([3 14])
xlim([5 165])
%colorbar

% traces
subplot(4,1,1)
plot(eegfilt(double(detrend(cutted(4,:))),srate,7.5,8.5));
xlim([5*srate 165*srate])
axis off
subplot(4,1,2)
plot(eegfilt(double(detrend(cutted(1,:))),srate,6.5,7.5));
xlim([5*srate 165*srate])
axis off

%compute welch

%plot(F,mean(P(:,1:583),2));
%hold on
%plot(F,mean(P(:,583:1166),2));
%plot(F,mean(P(:,1667:end),2));

%% traces
plot(eegfilt(double(detrend(cutted(4,:))),srate,8,9));


%%
%right spectrogram

SR = 1000;
%SRold=5000;

window=2*SR
overlap=1.9*SR
%data3(7,:)=eegfilt(double(data3(7,:)),SR,3.5,20);


[S,F,T,P] = spectrogram(data_ds(9,1:100000),window,overlap,2^18,SR);

imagesc(T,F,P);
colormap('jet');
axis xy
ylim([3 15])
cax=caxis;
caxis([0 cax(2)./10])

% 
% I=find(F>0);
% imagesc(T,F(I),P(I,:));
%subplot(2,1,1)

%%

P2=P;
t=find(T>56.5,1,'first');
t2=find(T>69);
P2(:,t:t+309)=P(:,672:981);

imagesc(T,F,P2);
colormap('jet');
axis xy
ylim([3 15])
cax=caxis;
caxis([0 cax(2)./20])

%%

I_f=find(F>3 & F<20);
I_t=find(T>40 & T<60);

P3=P;
P3(I_f,I_t)=P3(I_f,I_t).*2;
imagesc(T,F,P3);
colormap('jet');
axis xy
ylim([3 15])
cax=caxis;
caxis([0 cax(2)./10])
xlabel('Time (s)')
ylabel('Frequency (Hz)')
colorbar
caxis([0 9500])

%% plot
tmax_s=100;
tmin=1;
t_data_max=find(t_ds>tmax_s,1,'first'); 
t_adc_max=find(t_adc_ds>tmax_s,1,'first'); 

subplot(2,1,1);


data=data_ds(9,30*SR:50*SR);
data_f=eegfilt(double(detrend(data)),SR,3,14);
plot(t_ds(30*SR:50*SR),data_f);


xlim([0 tmax_s])

subplot(2,1,2);
imagesc(T,F,P);
colormap('jet');
axis xy
ylim([3 15])
xlim([0 tmax_s])
cax=caxis;
caxis([0 cax(2)./10])
xlabel('Time (s)')
ylabel('Frequency (Hz)')

%subplot(3,1,3);
%plot(t_adc_ds(tmin:t_adc_max),adc_ds(2,tmin:t_adc_max));
%ylabel('velocity')
%xlim([0 tmax_s]);

%ylabel('Power (\muV^2 Hz^{-1})')
%% normalize with T=0..3
t_norm=find(T<3);
norm=mean(P(:,t_norm),2);
P_norm=zeros(size(P,1),size(P,2));
for i=1:length(T)
    P_norm(:,i)=10*log10(P(:,i)./norm);
end
% plot normalized
imagesc(T,F,P_norm);
%colormap('jet');
axis xy
ylim([3 17])

%%
subplot(2,1,2)
[P1 F1]=pwelch(detrend(data1(11,1*SRold:25:5*SRold)),1000, [],2^12, 1000);
plot(F1,P1,'r');
hold on
[P2 F2]=pwelch(detrend(data1(12,10*SRold:25:15*SRold)),1000, [],2^12, 1000);
plot(F2,P2);
% [P3 F3]=pwelch(detrend(data3(7,20*SRold:25:25*SRold)),1000, [],2^12, 1000);
% plot(F3,P3,'k');
% [P4 F4]=pwelch(detrend(data3(7,40*SRold:25:60*SRold)),1000, [],2^12, 1000);
% plot(F3,P3,'g');
xlim([0 20]);
% subplot(3,1,3)
% plot(data3(7,5*SRold:25:25*SRold));




%%

%%