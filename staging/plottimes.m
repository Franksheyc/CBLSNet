function subs=plottimes(test_yy,test_y,y,x,leng,ss1,mode,num,plotflag)
%figs 画图类型
%mode 是全画还是只画一部分，mode输入为[x,y],x为类型，1为全画，2为画一部分
%num  为第num被试的图像
s0=leng(ss1);
for i=1:length(ss1)
    s(i)=sum(s0(1:i));
end
s=[0 s];

if mode(1)==1
    for i=1:length(ss1) 
        sta=s( i )+1;
        ter=s( i+1 );
        [acc,ac,f1,F1,x1,y1,test_yy1]=statisticsub(test_y,test_yy,y,x,sta,ter);
        subs(i).num=ss1(i);
        subs(i).ac=ac;
        subs(i).acc=acc;
        subs(i).f1=f1;
        subs(i).F1=F1;
        subs(i).weight=x1;
        subs(i).result=y1;
        subs(i).label=test_yy1;
        if plotflag
            figure
            plottts(test_yy,y,x,sta,ter,i);
        end
    end
else
    sta=s( num )+1;
    ter=s( num+1 );
    [acc,f1,x1]=statisticsub(test_y,test_yy,y,x,sta,ter);
    subs.num=ss1(i);
    subs.ac=ac;
    subs.acc=acc;
    subs.f1=f1;
    subs.F1=F1;
    subs.weight=x1;
    subs.result=y1;
    subs.label=test_yy1;
    if plotflag
    figure
    plottts(test_yy,y,x,sta,ter,num);
    end
end
        

        
end
function [acc,ac,f1,F1,x1,y1,test_yy1]=statisticsub(test_y,test_yy,y,x,sta,ter)
xlabel=sta:ter;
ac=length(find(test_yy(xlabel)==y(xlabel)))/length(xlabel);
[~,~,acc,~,~,f1,~,~,~] = Evaluation(y(xlabel),test_y(xlabel,:),5);
F1=sum(f1);
x1=x(xlabel,:);
y1=y(xlabel);
test_yy1=test_yy(xlabel);
end
function plottts(test_yy,y,x,sta,ter,num)
% switch figs
%     case 1
%         times=1;
%     case 2
%         weights=1;
% end
xlabel=sta:ter;
acc=length(find(test_yy(xlabel)==y(xlabel)))/length(xlabel);
xlabel=1:ter-sta+1;

hold on
subplot(3,1,1)
x1=x(sta:ter,:);
% x1=exp(x1)./sum(exp(x1)')';
% x2=x1-min(x1);
% x1=x2/sum(x2);
la1=test_yy(sta:ter,:);
y1=y(sta:ter,:);

bar(xlabel,x1,'stacked')
axis([1,ter-sta+1,min(min(x1)),max(max(x1))])
ylabel('weight')
title(['the ',num2str(num),' subject','    acc is :',num2str(acc)])
legend('N1','N2','N2','REM','WAKE','Location','best')

subplot(3,1,2)
plot(xlabel,la1)
axis([1,ter-sta+1,1,5])
ylabel('label')

subplot(3,1,3)
plot(xlabel,y1)
axis([1,ter-sta+1,1,5])
ylabel('result')


% % 画不一致的地方
% hold on
% h=find(y1~=la1);
% scatter(h,y1(h),'*')
end

