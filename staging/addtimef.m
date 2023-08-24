function [xx,yy]=addtimef(xx,yy,ftr,train_yy)
xx1=xx;yy1=yy;
% for i=1:length(xx)
% xx(i,yy(i))=xx(i,yy(i))*ftr(i,6)*0.6+xx(i,yy(i))*ftr(i,1);
% end
% xx(yy)=xx(yy).*ftr.meanx*0.6+xx(yy).*ftr.mean(yy);

d=find(train_yy~=yy);
for i=1:length(d)
xx(d(i),yy(i))=xx(d(i),yy(i))-ftr.meanx(i);
end
% xx(yy)=xx(yy).*ftr.meanx*0.6;
% xx=xx.*ftr.mean;
[~,yy]=max(xx');
yy=yy';
if ~isempty(train_yy)
    or=length(find(train_yy==yy1));
    ne=length(find(train_yy==yy));
else
    ne=2;or=1;
end
if ne>or
  xx=xx1;yy=yy1;
end
end