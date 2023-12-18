function NiceErrorBars(power_laser,power_nolaser,F,MINFREQ,MAXFREQ,titl,legend1,legend2,N,Npoints, seed)
% BARPLOT STATS PLOT FOR TWO PRE-STIM DATASETS
%
% mean_pre, mean_stim  ... matrices Frequenciews x Nsubjects
% F ... frequency
% MINFREQ ... low frequency bound
% MAXFREQ ... high frequency bound
% titl ... figure title
fi=figure;

f=find(F>MINFREQ & F<MAXFREQ);

%1) compute mean and std-dev
mean_laser=max(mean(power_laser(:,f)));
mean_nolaser=max(mean(power_nolaser(:,f)));
std_laser=max(std(power_laser(:,f)))/sqrt(N);
std_nolaser=max(std(power_nolaser(:,f)))/sqrt(N);

%2) generate points 
rng(seed)
points_laser=normrnd(mean_laser,std_laser,1,Npoints);
points_nolaser=normrnd(mean_nolaser,std_nolaser,1,Npoints);

%3) plot graph based on points
mu_laser=mean(points_laser);
sigma_laser=std(points_laser);
mu_nolaser=mean(points_nolaser);
sigma_nolaser=std(points_nolaser);
%h=barwitherr([std_laser,std_nolaser],[mean_laser, mean_nolaser]);
h=barwitherr([sigma_laser,sigma_nolaser],[mu_laser, mu_nolaser]);
set(h(1),'FaceColor',[0.5 0.5 1]);
hold on
plot(ones(1,Npoints),points_laser,'ok') 
plot(ones(1,Npoints).*2,points_nolaser,'ok') 

%4) labels
set(gca,'XTickLabel',{legend1,legend2})
box off
ylabel('Power (mV2/Hz)')
title(titl);
set(fi, 'Position', [100, 100, 175, 175]);

%5) stats
[hyp,p]=ttest2(points_laser,points_nolaser)