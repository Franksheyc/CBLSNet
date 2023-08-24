function [train_y,study_rate,si]=trainybanlance(kk,train_x,train_yy,train_y,study_rate,xx,yy,cr1,y,m,si)
xx1=zeros(length(train_x),5);
if kk==1
    thr=5;
    [sb,sp,study_rate,si]=bili(train_yy,y,thr,study_rate,si);
    xn=thecon(train_yy,10);
%     sp=setdiff(1:5,sp);
    for hj=1:length(sp)
        sp1=sp(hj);
        xx1(train_yy==sp1,sp1)=sb(sp1)*(sb(sp1)<=thr)+(thr*(sb(sp1)>thr));
    end
    train_y=train_y+study_rate*xx1.*cr1;%+xn*study_rate*0.5;
    %             xx1(train_yy==1,1)=sb(1);
    %             xx1(train_yy==3,3)=sb(3);
    
    %             xx(train_yy==2,2)=1.2;
else
    thr=4;
    [sb,sp,study_rate,si]=bili(yy,y,thr,study_rate,si);
%     xn=thecon(train_yy,5);
%     if length(sb)>4 && min(sb)/max(sb)>8 
        xx=exp(xx)./sum(exp(xx)')';
%     end
%     m=(quantile(cr1,m)+quantile(ce1,m))/2;  
    m=(quantile(cr1,m));
    m1=find(cr1<m);
    m1=intersect(m1,find(yy~=train_yy));
%     nn=[1.5 1.2 1.5 1 1];n=0.8;
%     [xx,~]=markovxx(train_yy,xx,nn,n);
%     xx1(m1,:)=xx(m1,:).*[1.5 1.05 1.5 1 1];
    xx1(m1,:)=xx(m1,:);
%     xx1=xx;
    for hj=1:length(sp)
        h1=intersect(find(train_yy==hj),m1);
        sp1=sp(hj);
        xx1(h1,sp1)= xx1(h1,sp1)*(sb(sp1) * ( sb(sp1)<=thr) + (thr*(sb(sp1)>thr)) );
    end
    
    train_y=train_y+study_rate*xx1.*cr1;%+xn*study_rate*0.5;%s

end
%         train_y=train_y+study_rate*xx;
%         train_y=train_y+study_rate*xx1.*mapminmax(cr1(:,2));

% train_y=train_y+study_rate*xx1.*s;%cr1(:,2);
% if study_rate>0.15
% study_rate=study_rate-0.02;
% end

end

function [sb,sp,study_rate,si]=bili(train_yy,y,thr,study_rate,si)
for i=1:5
    bb=length(find(train_yy==i))/length(train_yy);
    bb1=length(find(y==i))/length(y);
%     sb(i)=sqrt(1/bb);
    if bb1==0
        sb1(i)=0;
    else
        sb1(i)=min(1/bb1,15);
    end
    b1(i)=bb1;
    if bb==0
        sb(i)=0;
    else
        sb(i)=min(1/bb,15);
    end
    b(i)=bb;
    %                    xx1(train_yy==i,i)=sb;
end
di0=abs(sqrt(b)-sqrt(b1));
dis=sum(di0.*b);
% dis=sum(di0).*bb;
if dis>=0.07
    sb=sb1;
    b=b1;
    si=20;
elseif dis<0.07 && dis>0.05
    [~,h]=sort(di0);
    h=h(3:end);
    sb(h)=sb1(h);
    b(h)=b1(h);
    si=20;
end
% %%%new
% sb=sqrt(sb);
% %%%
[~,rankb]=sort(b);
h=rankb(end-1:end);
if ~isempty(find(b>0.5,1)) || min(b)<=0.05
   minsr=0.15;
else
    minsr=0.2;
end

if isempty(find(b>=0.5, 1))
%     if thr==4
        sp=find(sb>3);
%     else
%        sp=find(sb>2.5); 
%     end

%     sp=unique([1 3 sp]);
%     spp=ismember([1,3],sp);
%     spd=[1 3];
%     if sum(spp)<4
%         sp=[sp spd(spp==0)];
%     end
else
    sp=setdiff(1:5,find(b>=0.2));
    
end

sb1=sort(sb);
sb=sqrt(sb/((sb1(1)+sb1(2))/2));
% sb=sb/((sb1(1)+sb1(2))/2);

if sb(5)>sb(4)
    sb([2 4])=sqrt(sb([2 4]));
else
    sb([2 4 5])=sqrt(sb([2 4 5]));
end
% sb(h)=sqrt(sb(h));
while min(sb(sp))*study_rate<1 && study_rate<=min(0.5,minsr)
   study_rate=study_rate+0.02; 
end
study_rate=study_rate-0.02*(0.5*max(sb>thr));
end
function xx=thecon(train_yy,n)
train_y1=[zeros(n,1);train_yy;zeros(n,1)];
trainy=label2lab(train_yy);
for i=1:length(train_yy)
    y=train_y1(i:2*n+i);
    l=length(unique(y));
    if l==1
        c(i)=0;
    elseif l==2 && y(n+1)~=y(n) && y(n+1)~=y(n+2)
        c(i)=1;
    elseif l==2 && ((y(n+1)==y(n) && y(n+1)~=y(n+2)) || (y(n+1)~=y(n) && y(n+1)==y(n+2)))
        c(i)=0.25;
    elseif l>2 && y(n+1)~=y(n) && y(n+1)~=y(n+2)
        c(i)=0.5;
    else
        c(i)=0.25;
    end
end
xx=trainy.*c';
end