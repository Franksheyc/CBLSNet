function [cr,ce]=crce(xx,x,si,bs)
  [ae,y]=sort(x,2);
  [ar,yy]=sort(xx,2);
  
% %   xx=putpos(xx,ar(:,1));
% %   x=putpos(x,ae(:,1));
% %   [ae,~]=sort(x,2);
% %   [ar,~]=sort(xx,2);
%   ce(:,1)=ae(:,5)./ae(:,4);
%   cr(:,1)=ar(:,5)./ar(:,4);
  
  ce=ae(:,5)./sum(ae,2);
  cr=ar(:,5)./sum(ar,2);
  
  by=findcc(y(:,5),si,bs);
  byy=findcc(yy(:,5),si,bs);
  cr=cr.*byy';
  ce=ce.*by';
end

% function x=putpos(x,ae)
% s=find(ae<0);
% s1=min(x(s,:)');
% s1=s1';
% x(s,:)=x(s,:)-ones(length(s),5).*s1;
% end

function b=findcc(label,n,be)
b=ones(1,n);
parfor i=n+1:length(label)-n
    a=unique(label(i-n:i+n));
    switch length(a)
        case 1
            b(i)=be(1);
        case 2
            b(i)=be(2);
        case 3
            b(i)=be(3);
        case 4
            b(i)=be(4);
        case 5
            b(i)=be(5);
    end
end
b(end+1:end+n)=ones(1,n);
end