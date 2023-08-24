function [C11,C12,C01,C02,M,CM1,CM2,C21,C22,di]=feature_en(C3,C4,C3t,C4t,EMG,C32,C42,EMG2)
   

[C11,c1d1]=bilif(C3t,C3t,[1 2  6],2);
[C12,c1d2]=bilif(C4t,C4t,[1 2  6],2);
[M,md]=bilif(EMG,EMG,[1 2  6],1);
[C01,c0d1]=bilif(C32(:,27:end),C32(:,27:end),[1 2  6],4);%
[C02,c0d2]=bilif(C42(:,27:end),C42(:,27:end),[1 2  6],4);
[CM1,cmd1]=bilif(C32(:,27:end),EMG2,[1  2  6],3);
[CM2,cmd2]=bilif(C42(:,27:end),EMG2,[1  2  6],3);

% C31=C3(:,1:50);
C32=C3(:,51:end);
% C41=C4(:,1:50);
C42=C4(:,51:end);
% [C21,c2d1]=bilif(C32,C32,[1:36 37:72 216:252],3);
% [C22,c2d2]=bilif(C42,C42,[1:36 37:72 216:252],3);
[C21,c2d1]=bilif(C32,C32,[1:6 7:12 37:42],3);
[C22,c2d2]=bilif(C42,C42,[1:6 7:12 37:42],3);

di=sort(unique([c1d1;c1d2;md;c0d1;c0d2;cmd1;cmd2;c2d1;c2d2]));
st=zeros(1,200);epsino=0.1;
[C11,~]=deleteoutliners_persub(C11,epsino,st);
[C12,~]=deleteoutliners_persub(C12,epsino,st);
[C01,~]=deleteoutliners_persub(C01,epsino,st);
[C02,~]=deleteoutliners_persub(C02,epsino,st);
[CM1,~]=deleteoutliners_persub(CM1,epsino,st);
[CM2,~]=deleteoutliners_persub(CM2,epsino,st);
[C21,~]=deleteoutliners_persub(C21,epsino,st);
[C22,~]=deleteoutliners_persub(C22,epsino,st);
[M,~]=deleteoutliners_persub(M,epsino,st);
C11=zscore(C11);
C12=zscore(C12);
C01=zscore(C01);
C02=zscore(C02);
C21=zscore(C21);
C22=zscore(C22);
CM1=zscore(CM1);
CM2=zscore(CM2);
M=zscore(M);
% EMG3(:,1)=[EMG(:,1);EMG(end,1)]./[EMG(1,5);EMG(:,5)];
% EMG3(:,2)=[EMG(:,16);EMG(end,16)]./[EMG(1,20);EMG(:,20)];
% 
% C43(:,1)=[C4(:,51);C4(end,51)]./[C4(1,55);C4(:,55)];
% C43(:,2)=[C4(:,69);C4(end,69)]./[C4(1,74);C4(:,74)];
% C43=C43(1:end-1,:);EMG3=EMG3(1:end-1,:);
% C43=mapminmax(C43')';EMG3=mapminmax(EMG3')';
end

function [E,di]=bilif(sig1,sig2,s,mode)
if mode==1
    n=5;%emg/emg
elseif mode==2
    n=5;%eeg/eeg
elseif mode==3
    n=0;%eeg/emg
else
    n=1;
end
[l1,~]=size(sig1);

for i=1:length(s)
    m=ones(l1,1);
    h=(s(i)-1)*(n+1)+1;
    h1=h+n;
    d=[sig1(:,h);sig1(end,h)]./[sig2(1,h1);sig2(:,h1)];
    m0=find(d(2:end-1)<0);
    m(m0)=-1;
    d=sqrt(abs(d(1:end-1))).*m;
    di=find(abs(d)>mean(abs(d))+3*std(d));
    E(:,i)=d;
end
% E=zscore(E);
end