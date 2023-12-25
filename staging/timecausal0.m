function [f,f1]=timecausal0(time_related,leng,n)
%time_related 时间窗口
% leng  样本长度
% n  特征数
K=[];
for i=1:length(leng)
    K(i)=sum(leng(1:i));
end
K=[0 K];
f=[];f1=[];
for j=1:length(K)-1
    for i=1:time_related
        f=[f (K(j)+1:K(j)+time_related-(i-1)) + K(end)*(i-1)];
        f1=[f1 (K(j+1):-1:K(j+1)-(i-1)) + K(end)*(i-1+time_related+1)];
    end
end
f=[f f1];
f1=[];
for i=1:n
f1=[f1 f+(i-1)*(time_related*2+1)*K(end)];
end
end