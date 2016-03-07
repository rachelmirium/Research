[names,F1,F2,F3,F4]=textread('formants.txt','%s%d%d%d%d',10)

subplot(1,2,1);

plot(F2,F1,'r')

for i=1:10
text(F2(i),F1(i),names{i})
end;

set(gca,'xdir','reverse')
set(gca,'ydir','reverse')

xlabel('F2');
ylabel('F1');

grid on;

subplot(1,2,2);

plot(F2-F1,F1,'r')

for i=1:10
text(F2(i)-F1(i),F1(i),names{i})
end;

set(gca,'xdir','reverse')
set(gca,'ydir','reverse')

xlabel('F2-F1');
ylabel('F1');

grid on;