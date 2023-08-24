function [C3,C4,C3t,C4t,EOG,EMG,H,leng,label]=dedata(C3,C4,C3t,C4t,EOG,EMG,H,leng,label,mode)
label(label==0)=5;
% [C3,C4,C3t,C4t,EOG,EMG,label,leng]=demark_label(C3,C4,C3t,C4t,EOG,EMG,leng,label);

C31=C3(:,1:50);C32=C3(:,51:end);
C41=C4(:,1:50);C42=C4(:,51:end);

[c31th0,c31thr]=acthr(C31,leng);
[c32th0,c32thr]=acthr(C32,leng);
[c41th0,c41thr]=acthr(C41,leng);
[c42th0,c42thr]=acthr(C42,leng);
[c3tth0,c3tthr]=acthr(C3t,leng);
[c4tth0,c4tthr]=acthr(C4t,leng);
[emgth0,emgthr]=acthr(EMG,leng);
[eogth0,eogthr]=acthr(EOG,leng);

[c31s,c31th,c3t1]=demark(C31,leng,'EEG',[1 2 3 4],c31th0,c31thr);
[c32s,c32th,c3t2]=demark(C32,leng,'EEG2',[2 5],c32th0,c32thr);
[c41s,c41th,c4t1]=demark(C41,leng,'EEG',[1 2 3 4],c41th0,c41thr);
[c42s,c42th,c4t2]=demark(C42,leng,'EEG2',[2 5],c42th0,c42thr);

[c3ts,c3tth,c3tt]=demark(C3t,leng,'EEG',[1 2 3 4],c3tth0,c3tthr);
[c4ts,c4tth,c4tt]=demark(C4t,leng,'EEG',[1 2 3 4],c4tth0,c4tthr);

[emgs,emgth,emgt]=demark(EMG,leng,'EMG',[1 2 3 4],emgth0,emgthr);
[eogs,eogth,eogt]=demark(EOG,leng,'EOG',[1 2 3 4 9 10 11 12],eogth0,eogthr);
c3s=[c31s' c32s'];c3th=[c31th c32th];c3t=[c3t1 c3t2+50];
c4s=[c41s' c42s'];c4th=[c41th c42th];c4t=[c4t1 c4t2+50];
clear C31 C32 C41 C42

[C3,C4,C3t,C4t,EOG,EMG,H,leng,label]=standdata(C3,C4,C3t,C4t,EOG,EMG,leng,label,c3s,c4s,emgs,eogs,c3th,c4th,emgth,eogth,c3t,c4t,eogt,emgt,c3ts,c3tth,c3tt,c4ts,c4tth,c4tt,H,mode);
label(label==5)=0;
end

function [thrr,thrr1]=acthr(d1,leng)
l1=1;l2=0;
for i=1:length(leng)
    l2=l2+leng(i);
    datw(i,:)=sum(d1(l1:l2,:));
    l1=l1+leng(i);
end
for i=1:size(datw,2)
    ss=quantile(datw(:,i),0.99);
    ss1=find(datw(:,i)<ss);
    md=abs(mean(datw(ss1,i)));
    sd=std(datw(ss1,i));
    thr=md+3*sd;
    thrr(i)=thr;
    thrr1(i)=md;
end

end

function [C3,C4,C3t,C4t,EOG,EMG,H,leng,label]=standdata(C3,C4,C3t,C4t,EOG,EMG,leng,label,c3s,c4s,emgs,eogs,c3th,c4th,emgth,eogth,c3t,c4t,eogt,emgt,c3ts,c3tth,c3tt,c4ts,c4tth,c4tt,H,mode)

for i=1:length(leng)
   len(i)=sum(leng(1:i));
end
len=[0 len];

if mode==1
eegs=union(c3s,c4s);
ots=union(eogs,emgs);

de=intersect(eegs,ots);
eogs1=setdiff(eogs,union(eegs,emgs));
emgs1=setdiff(emgs,union(eegs,eogs));
c3s1=setdiff(c3s,union(ots,c4s));
c4s1=setdiff(c4s,union(ots,c3s));
ont=setdiff(union(intersect(c3s,ots),intersect(c4s,ots)),de);
for i=1:length(c3s1)
    m=len(c3s1(i))+1:len(c3s1(i)+1);
    c3th1=abs(sum(C3(m,c3t)))./c3th;
    C3(m,c3t)=C3(m,c3t)./c3th1;
end
for i=1:length(c4s1)
    m=len(c4s1(i))+1:len(c4s1(i)+1);
    c4th1=abs(sum(C4(m,c4t)))./c4th;
    C4(m,c4t)=C4(m,c4t)./c4th1;
end
for i=1:length(eogs1)
    m=len(eogs1(i))+1:len(eogs1(i)+1);
    eogth1=abs(sum(EOG(m,eogt)))./eogth;
    EOG(m,eogt)=EOG(m,eogt)./eogth1;
end
for i=1:length(emgs1)
    m=len(emgs1(i))+1:len(emgs1(i)+1);
    emgth1=abs(sum(EMG(m,emgt)))./emgth;
    EMG(m,emgt)=EMG(m,emgt)./emgth1;
end
for i=1:length(de)
    m=len(de(i))+1:len(de(i)+1);
    c3th1=abs(sum(C3(m,c3t)))./c3th;
    c4th1=abs(sum(C4(m,c4t)))./c4th;
    eogth1=abs(sum(EOG(m,eogt)))./eogth;
    emgth1=abs(sum(EMG(m,emgt)))./emgth;
    c3tth1=abs(sum(C3(m,c3tt)))./c3tth;
    c4tth1=abs(sum(C4(m,c4tt)))./c4tth;
    C3(m,c3t)=C3(m,c3t)./c3th1;
    C4(m,c4t)=C4(m,c4t)./c4th1;
    EOG(m,eogt)=EOG(m,eogt)./eogth1;
    EMG(m,emgt)=EMG(m,emgt)./emgth1;
    C3t(m,c3tt)=C3t(m,c3tt)./c3tth1;
    C4t(m,c4tt)=C4t(m,c4tt)./c4tth1;
end
for i=1:length(ont)
    m=len(ont(i))+1:len(ont(i)+1);
    c3th1=abs(sum(C3(m,c3t)))./c3th;
    c4th1=abs(sum(C4(m,c4t)))./c4th;
    eogth1=abs(sum(EOG(m,eogt)))./eogth;
    emgth1=abs(sum(EMG(m,emgt)))./emgth;
    c3tth1=abs(sum(C3(m,c3tt)))./c3tth;
    c4tth1=abs(sum(C4(m,c4tt)))./c4tth;
    C3(m,c3t)=C3(m,c3t)./c3th1;
    C4(m,c4t)=C4(m,c4t)./c4th1;
    EOG(m,eogt)=EOG(m,eogt)./eogth1;
    EMG(m,emgt)=EMG(m,emgt)./emgth1;
    C3t(m,c3tt)=C3t(m,c3tt)./c3tth1;
    C4t(m,c4tt)=C4t(m,c4tt)./c4tth1;
end
else
    m0=unique([c3s c4s emgs eogs']);
    mf=[];
    for i=1:length(m0)
        mf=[mf len(m0(i))+1:len(m0(i)+1)];
    end
    m=setdiff(1:length(C3),mf);
    C3=C3(m,:);
    C4=C4(m,:);
    C3t=C3t(m,:);
    C4t=C4t(m,:);
    EOG=EOG(m,:);
    EMG=EMG(m,:);
    H=H(m,:);
    label=label(m);
    leng=leng(setdiff(1:length(leng),m0));
end

end

function [C3,C4,C3t,C4t,EOG,EMG,label,leng]=demark_label(C3,C4,C3t,C4t,EOG,EMG,leng,label)
% C31=[];C41=[];EOG1=[];EMG1=[];lab=[];leng0=[];
l1=1;l2=0;
k=0;mf=[];m0=[];
for i=1:length(leng)
    len(i)=sum(leng(1:i));
    l2=l2+leng(i);
    label1=label(l1:l2);
%     c3=C3(l1:l2,:);
%     c4=C4(l1:l2,:);
%     eog=EOG(l1:l2,:);
%     emg=EMG(l1:l2,:);
%     leng1=leng(i);
    for j=1:5
        bi(j)=length(find(label1==j))/length(label1);
    end
    if length(find(bi==0))>=3 || length(find(bi>=0.9))>=1
        k=k+1;
        m0(k)=i;
%         label1=[];
%         c3=[];
%         c4=[];
%         eog=[];
%         emg=[];
%         leng1=0;
    end
%     C31=[C31;c3];
%     C41=[C41;c4];
%     EOG1=[EOG1;eog];
%     EMG1=[EMG1;emg];
%     lab=[lab;label1];
%     leng0=[leng0 leng1];
    l1=l1+leng(i);
end
% leng0=leng0(leng0~=0);
    len=[0,len];
    for i=1:length(m0)
        mf=[mf len(m0(i))+1:len(m0(i)+1)];
    end
    m=setdiff(1:length(C3),mf);
    leng=leng(setdiff(1:length(leng),m0));
    C3=C3(m,:);
    C4=C4(m,:);
    C3t=C3t(m,:);
    C4t=C4t(m,:);
    EOG=EOG(m,:);
    EMG=EMG(m,:);
    label=label(m);
end
function [dimarks,thrr1,M]=demark(dat,leng,m,dd,thrr,thrr0)
switch m
    case 'EEG'
        n=6;
    case 'EEG2'
        n=6;
    case 'EOG'
        n=6;%2;
    case 'EMG'
        n=6;%5;
end
l1=1;l2=0;thrr1=[];M=[];
% for i=1:length(leng)
% lenn(i)=sum(leng(1:i));
% end
% lenn=[0 lenn];
% lenn=lenn(2:end);
for i=1:length(leng)
    l2=l2+leng(i);
    datw(i,:)=sum(dat(l1:l2,:));
    l1=l1+leng(i);
end
% dd=[1 2 5];%mean std eng
% if n==2
%     dd=[dd dd+16];
% end
dimark=zeros(length(leng),1);
for i=1:length(dd)
    mark=(dd(i)-1)*n+1:dd(i)*n;
    d1=abs(datw(:,mark));
    thr=thrr(mark);
    thr1=thrr0(mark);
    dma=d1>thr;
    dma1=sum(dma')';
    di=find(dma1>=n/2);
    dimark(di)=dimark(di)+1;
    thrr1=[thrr1 thr1];
    M=[M mark];
end
dimarks=find(dimark>=length(dd/2));

end