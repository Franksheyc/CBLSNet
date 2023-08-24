function [C31,C41,EOG1,EMG1,H1]=datanorms(C3,C4,EOG,EMG,H)
C31=eegnorm(C3);
C41=eegnorm(C4);
% F41=eegnorm(F4);
% O21=eegnorm(O2);
EOG1=eognorm(EOG);
EMG1=tfsnorm(EMG,'EMG');
H1=Hnorm(H);
% EMG1=[EMG1(:,1:26) EMG1(:,31:end)];
% EOG1=[EOG1(:,1:11) EOG1(:,13:29) EOG1(:,31:38) EOG1(:,40:41)];
end

function C3=eegnorm(C3)
C31=coenorm(C3(:,1:50));
C32=tfsnorm(C3(:,51:end),'EEG2');
C3=[C31,C32];
end
function EOG=eognorm(EOG)
EOG1=tfsnorm(EOG(:,1:54),'EOG');
EOG2=tfsnorm(EOG(:,55:108),'EOG');
EOG=[EOG1 EOG2 EOG(:,109:end)];
end

function dat0=coenorm(dat)
len=6;
dat0=[];
dat=dat(:,1:26);
means=dat(:,19:24);
ra1=ratio(1,sqrt(means));
ra2=ratio(2,sqrt(means));
ra3=ratio(3,sqrt(means));
ra4=ratio(4,sqrt(means));
ra5=ratio(5,sqrt(means));
ra6=ratio(6,sqrt(means));
% ra=[ra1 ra2 ra3 ra4 ra5 ra6];
ra=[ra1(:,end) ra2(:,end) ra3(:,end) ra4(:,end) ra5(:,end) ra6(:,end)];
% softmax1=@(x)exp(x)./sum(exp(x));
% % sigmoid=@(x)1./(1+exp(-x));
% ra1=softmax1(ra1(:,1:5)')';
% ra2=softmax1(ra2(:,1:5)')';
% ra3=softmax1(ra3(:,1:5)')';
% ra4=softmax1(ra4(:,1:5)')';
% ra5=softmax1(ra5(:,1:5)')';
% ra6=softmax1(ra6(:,1:5)')';
% % ra1=nor_ra(ra1(:,1:5));
% % ra2=nor_ra(ra2(:,1:5));
% % ra3=nor_ra(ra3(:,1:5));
% % ra4=nor_ra(ra4(:,1:5));
% % ra5=nor_ra(ra5(:,1:5));
% % ra6=nor_ra(ra6(:,1:5));
dat00=[ra ra1(:,1:end-1) ra2(:,1:end-1) ra3(:,1:end-1) ra4(:,1:end-1) ra5(:,1:end-1) ra6(:,1:end-1)];
% dat00=[ra ra1 ra2 ra3 ra4 ra5 ra6];
% dat(:,19:24)=log(abs(dat(:,19:24)));
for i=1:3
    dat1=zscore(dat(:,len*(i-1)+1+2*(i-4>=1):len*i+2*(i-4>=1)));
    dat0=[dat0 dat1];
end
clear dat1
i=4;
dat1=zscore(log(sqrt(dat(:,len*(i-1)+1+2*(i-4>=1):len*i+2*(i-4>=1)))));
dat10=dat1;
clear dat1
% for i=5:8
%     dat1=abs( dat(:,len*(i-1)+1+2*(i-4>=1):len*i+2*(i-4>=1)) ./ max( dat(:,len*(i-1)+1+2*(i-4>=1):len*i+2*(i-4>=1))' )' );
%     dat00=[dat00 dat1];
% end
dat01(:,1)=zscore(dat(:,25));
dat01(:,2)=zscore(dat(:,26));
dat0=[dat0 dat10 dat01 dat00];

end

function dat=tfsnorm(dat,k)
switch k
    case 'EEG'
        len=6;
    case 'EEG2'
        len=6;
    case 'EOG'
%         len=2;
        len=6;
    case 'EMG'
%         len=5;
        len=6;
end
% dat(:,len*4+1:len*5)=log(abs(dat(:,len*4+1:len*5)));
s=9;
% [m,n]=size(dat);
% dat0=zeros(m,n);
l=[1 2 3 4 5 6  7 8 9];
for i=1:length(l)
% dat1=mapminmax(dat(:,(i-1)*len+1:i*len));
if i==5
    dat1=zscore( log(sqrt( dat(:,(l(i)-1)*len+1:l(i)*len) )));
% elseif i==9
%     dat1=zscore( exp( dat(:,(l(i)-1)*len+1:l(i)*len) ));
% elseif i==7
%     dat1=dat(:,(l(i)-1)*len+1:l(i)*len);
else
    dat1=zscore(dat(:,(l(i)-1)*len+1:l(i)*len));
end
% dat0(:,len*(l(i)-1)+1:len*l(i))=dat1;
dat(:,len*(l(i)-1)+1:len*l(i))=dat1;
end
% dat0(:,len*6+1:len*7)=dat(:,len*6+1:len*7);
end

function ra=ratio(n,means)
m=setdiff(1:size(means,2),n);
for i=1:length(m)
    ra(:,i)=means(:,n)./means(:,m(i));
end
ra(:,i+1)=means(:,n)./sum(means')';
end

function H1=Hnorm(H0)
H1=[];
for i=1:5
    H=H0(:,(i-1)*18+1:i*18);
    ha=mean(H(:, ((1:6)-1)*3+1 )')';
    hm=mean(H(:, ((1:6)-1)*3+2 )')';
    hc=mean(H(:, ((1:6)-1)*3+3 )')';
    ha=zscore(ha);
    hm=zscore(hm);
    hc=zscore(hc);
    H1=[H1 ha hm hc];
end
end

function a=nor_ra(a)
a=abs(a)/10^3;
if max(max(a))
    a=a./max(a')';
end
end