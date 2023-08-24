function [continues,transform,stage_ratio,current_ratio,yk]=weight_features(Y,n)
s=n+1;
for i=1:n
    eval('a1(:,i)=Y(:,s-(i-1))==Y(:,s-i);');
    eval('a2(:,i)=Y(:,s+(i-1))==Y(:,s+i);');
end
% a1=Y(:,n+1)==Y(:,n);
% a2=Y(:,n+1)==Y(:,n+2);
% a=a1+a2;
% transepoch=(a+5)/10;
% for i=1:size(Y,2)-1
%     a(:,i)=Y(:,i)==Y(:,i+1);
% end
% a=sum(a')';
a=sum((a1+a2)')';
continues=a/(2*n+1);
transform=1-continues;
parfor j=1:length(Y)
    for i=1:5
        stage_ratio(j,i)= length ( find(Y(j,:)==i) )/(2*n+1) +0.1 ;
    end
    current_ratio(j)=length( find( Y(j,:)==Y(j,n+1) ) )/ (2*n+1);
    yk(j)=1/length(unique(Y(j,:)));
end
current_ratio=current_ratio';
yk=yk';
% transepoch=gather(transepoch);
% continues=gather(continues);
% transform=gather(transform);
% stage_ratio=gather(stage_ratio);
% current_ratio=gather(current_ratio);