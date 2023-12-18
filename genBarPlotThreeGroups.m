function genBarPlotThreeGroups(group1, group2, group3,ytext,xtext,group1text, group2text, group3text)

%compute mean and std
mean1=mean(group1);
std1=std(group1);
mean2=mean(group2);
std2=std(group2);
mean3=mean(group3);
std3=std(group3);


%plot bars with erros
fi=figure
h=barwitherr([std1,std2, std3],[mean1, mean2, mean3])
set(h(1),'FaceColor',[0.5 0.5 1]);

%plot points
hold on
plot(ones(1,length(group1)),group1,'ok') 
plot(ones(1,length(group2)).*2,group2,'ok') 
plot(ones(1,length(group3)).*3,group3,'ok') 

% labels
set(gca,'XTickLabel',{group1text,group2text,group3text})
box off
xlabel(xtext)
ylabel(ytext)

set(fi, 'Position', [100, 100, 175*2, 175]);