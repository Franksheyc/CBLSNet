function [continues,transform,stage_ratio,current_ratio,yk]=weight_features(Y,n)
%����y��yy���жϴ������������ԣ������ԣ��ڵ�ǰ�����и��׶�ռ�ȣ���ǰ֡�׶��ڴ�����ռ��
% Y=gpuArray(single(Y));
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
    for i=1:5
        a=ismember(Y,i);
        stage_ratio(:,i)=sum(a')';%ÿ�������ڴ����еı���
    end
    yk=ismember(stage_ratio,0);
    yk=1/(6-sum(yk')');%���������������
    yk=yk';

    y=Y(:,n+1);
    y=label2lab(y);
    current_ratio=stage_ratio.*y/(2*n+1);%�뵱ǰ֡��ͬ�������ڴ����еı���
    current_ratio=sum(current_ratio')';
    stage_ratio=stage_ratio/(2*n+1) +0.1 ;

% parfor j=1:length(Y)    
%     for i=1:5
%         stage_ratio(j,i)= length ( find(Y(j,:)==i) )/(2*n+1) +0.1 ;
%     end
% 
%     current_ratio(j)=length( find( Y(j,:)==Y(j,n+1) ) )/ (2*n+1);
%     yk(j)=1/(length(unique(Y(j,:)))-1);
% end
% current_ratio=current_ratio';
% yk=yk';
% transepoch=gather(transepoch);
% continues=gather(continues);
% transform=gather(transform);
% stage_ratio=gather(stage_ratio);
% current_ratio=gather(current_ratio);
