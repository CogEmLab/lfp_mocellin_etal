load 53010_smell_ds
%load 51040_smell_ds
%load 53055_smell_ds
%load 53056_smell_ds

ethovis(:,3)=ethovis(:,3)+7.5; %normalize
%plot(ethovis(:,3),ethovis(:,4));

%% stats
%speed in center
Icenter=find(sqrt(ethovis(:,3).^2+ethovis(:,4).^2)<5);

%speed in intermediate
Iinter=find(sqrt(ethovis(:,3).^2+ethovis(:,4).^2)>5 & sqrt(ethovis(:,3).^2+ethovis(:,4).^2)<10);

%speed in ext. inter
Iext=find(sqrt(ethovis(:,3).^2+ethovis(:,4).^2)>10 & sqrt(ethovis(:,3).^2+ethovis(:,4).^2)<15);

%speed in wall
Iwall=find(sqrt(ethovis(:,3).^2+ethovis(:,4).^2)>15);

subplot(1,4,1)
boxplot(ethovis(Icenter,5)); title('center')
subplot(1,4,2)
boxplot(ethovis(Iinter,5)); title('intermediate')
subplot(1,4,3)
boxplot(ethovis(Iext,5)); title('ext. inter')
subplot(1,4,4)
boxplot(ethovis(Iwall,5)); title('wall')

