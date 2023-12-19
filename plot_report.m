plot(weak,'-o')
hold on
plot([1,4],[weak(1),weak(1)],'--')
legend('horizontal scaling','ideal')
set(gca,'XTickLabel',1:4)
set(gca,'XTickLabel',{'1 (2)','2 (4)',' 4 (8)','8 (16)'})
box off
xlabel('number of nodes (number of cores)')
ylabel('execution time (s)')

%% normalized
plot(weak./weak(1),'-o')
hold on
plot([1,4],[1,1],'--')
legend('horizontal scaling','ideal')
set(gca,'XTickLabel',1:4)
set(gca,'XTickLabel',{'1 (2)','2 (4)',' 4 (8)','8 (16)'})
box off
xlabel('number of nodes (number of cores)')
ylabel('normalized execution time')