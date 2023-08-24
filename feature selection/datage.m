function [FF,CFF,f,f0]=datage(C3,C4,c3,c4,EOG,EMG,...
    C3t,C4t,C32,C42,EMG2,EOG2,C11,C12,C01,C02,CM1,CM2,...
    C21,C22,H,mode,bu,mo)
load(filename1);
load(filename2);
% EOG=EOG(:,[1:36 37 38 40 41]);
% EOG=EOG(:,[1:32 33 34 36 37]);

cm=length(C4);
% % fe=othersignal(da);
% se=chosefea(1,[1 2 4 5 6 8],62);
% C40=[C4(:,1:12)*2 C4(:,13:18)*2 C4(:,19:26)*2 C4(:, [63:74  81:86  87:92 93:98 105:110])];%63:74 87:92      [63:74   87:92 ]
% C30=[C3(:,1:12)*2 C3(:,13:18)*2 C3(:,19:26)*2 C3(:, [63:74  81:86  87:92 93:98 105:110])];%75:80 81:86  93:98 99:104 105:110  111:116    75:80 81:86 
C40=[C4(:,1:11)*2 C4(:,12) C4(:,13:18)*2 C4(:,19:23)*2 C4(:,24) C4(:,25)*2 C4(:,26) C4(:, [63:74  81:86  87:92 93:98 105:110])];%63:74 87:92      [63:74   87:92 ]      C4(:,24) C4(:,25)*2 C4(:,26)
C30=[C3(:,1:11)*2 C3(:,12) C3(:,13:18)*2 C3(:,19:23)*2 C3(:,24) C3(:,25)*2 C3(:,26) C3(:, [63:74  81:86  87:92 93:98 105:110])];%75:80 81:86  93:98 99:104 105:110  111:116    75:80 81:86 
% C41=[C4(:,27:32) C4(:,33:37).*(rand(cm,1)/10+0.9) C4(:,38:42).*(rand(cm,1)/10+0.9) C4(:,45) C4(:,48:52).*(rand(cm,1)/10+0.9) C4(:,53:57).*(rand(cm,1)/10+0.9) C4(:,58:62).*(rand(cm,1)/10+0.9)];
% C31=[C3(:,27:32) C3(:,33:37).*(rand(cm,1)/10+0.9) C3(:,38:42).*(rand(cm,1)/10+0.9) C3(:,45) C3(:,48:52).*(rand(cm,1)/10+0.9) C3(:,53:57).*(rand(cm,1)/10+0.9) C3(:,58:62).*(rand(cm,1)/10+0.9)];

C41=[C4(:,[27:32])*2 [C4(:,[33:37]) C4(:,[38:42]) C4(:,[43:47]) C4(:,[48:52]) C4(:,[53:57]) C4(:,[58:62])]/2];
C31=[C3(:,[27:32])*2 [C3(:,[33:37]) C3(:,[38:42]) C3(:,[43:47]) C3(:,[48:52]) C3(:,[53:57]) C3(:,[58:62])]/2];
% C41=[C4(:,27:32)*4  C4(:,33:35) C4(:,36:37)/2  C4(:,38:40) C4(:,41:42)/2 C4(:,43:45) C4(:,46:47)/2  C4(:,48:50) C4(:,51:52)/2 C4(:,53:57) C4(:,58:62)]/2;%sleep edf
% C31=[C3(:,27:32)*4  C3(:,33:35) C3(:,36:37)/2  C3(:,38:40) C3(:,41:42)/2 C3(:,43:45) C4(:,46:47)/2  C3(:,48:50) C3(:,51:52)/2 C3(:,53:57) C3(:,58:62)]/2;

% c1=[c4(:,1:6) c4(:,7:12) c4(:,13:18) c4(:,19:24)].*(rand(cm,1)/10+0.7);
% c2=[c3(:,1:6) c3(:,7:12) c3(:,13:18) c3(:,19:24)].*(rand(cm,1)/10+0.7);
c1=[c4(:,7:12) c4(:,13:18) c4(:,19:24)].*(rand(cm,1)/10+0.7);% Sleep-edf
c2=[c3(:,7:12) c3(:,13:18) c3(:,19:24)].*(rand(cm,1)/10+0.7);
% c42=zscore(c41')';
% c31=zscore(c31')';
% c5=(c1+c2)/2;

C40=zscore(C40')';
% C41=zscore(C41')';
C30=zscore(C30')';
% C31=zscore(C31')';

EMG=zscore(EMG')';
EOG=zscore(EOG')';
C42=zscore(C42')';
C32=zscore(C32')';
EMG2=zscore(EMG2')';
EOG2=zscore(EOG2')';
C3t=zscore(C3t')';
C4t=zscore(C4t')';


% C11=zscore(C11')';
% CM1=zscore(CM1')';
% C21=zscore(C21')';
% C01=zscore(C01')';
% C12=zscore(C12')';
% CM2=zscore(CM2')';
% C22=zscore(C22')';
% C02=zscore(C02')';

FE1=[C40 C41 C42(:,[1:26  27:end])/2 C4t(:,[1:12 25:30 31:36 37:42 43:48])/2];%C42/2 C4t/2  (:,7:end) 25:30
FE2=[C30 C31 C32(:,[1:26  27:end])/2 C3t(:,[1:12 25:30 31:36 37:42 43:48])/2];%[ 16:18 21:23 24:26]
FO1=[EOG(:,1:54) EOG(:,[109 110]) EOG2(:,1:16)];%
FO2=[EOG(:,55:108) EOG(:,[111 112]) EOG2(:,17:32)];%
% FO1=[EOG(:,1:16) EOG(:,[33 34]) EOG2(:,1:16)];%
% FO2=[EOG(:,17:32) EOG(:,[35 36]) EOG2(:,17:32)];%
CFF1=[C01 C11 CM1/2 C21];
CFF2=[C02 C12 CM2/2 C22];
CFFO1=[C01 C11  C21];
CFFO2=[C02 C12  C22];

[H,~]=deleteoutliners_persub(H,0.1,ones(1,size(H,2)));
H=zscore(H);

% if ~isnan(H(:,3))
% HE1=H(:,1:3);
% HE2=H(:,4:6);
% HO1=H(:,7:9);
% HO2=H(:,10:12);
% HM=H(:,13:15);
% else
HE1=H(:,1:2);
HE2=H(:,4:5);
HO1=H(:,7:8);
HO2=H(:,10:11);
HM=H(:,13:14);
% end

kl=2;
for i=1:9*6/kl
    FO1p(:,i)=sum(FO1(:,(i-1)*kl+1:i*kl)')';
    FO2p(:,i)=sum(FO2(:,(i-1)*kl+1:i*kl)')';
    FMp(:,i)=sum(EMG(:,(i-1)*kl+1:i*kl)')';
end
so=chosefea(kl,[1:9],0);
FO1=[FO1p(:,so) FO1(:,57:end) ];
FO2=[FO2p(:,so) FO2(:,57:end) ];
FM=[FMp EMG2];
% FE3=(FE1+FE2)/2;
% FO3=(FO1+FO2)/2;
% CFF3=(CFF1+CFF2)/2;
% HE3=(HE1+HE2)/2;
% HO3=(HO1+HO2)/2;
clear H
% mo=2;
if mode==1
    [FE,FO,CFF,CFFO,HE,HO,c5]=chochannel(leng,FE1,FE2,FO1,FO2,CFF1,CFF2,CFFO1,CFFO2,HE1,HE2,HO1,HO2,c1,c2,mo);
%     [FE,FO,CFF,H,c5]=chochannel(leng,FE3,FE3,FO3,FO3,CFF3,CFF3,HE3,HE3,HO3,HO3,c5,c5,mo);
else
    FE=[FE1 FE2];
    FO=[FO1 FO2];
    CFF=[CFF1 CFF2];
    CFFO=[CFFO1 CFFO2];
    HE=[HE1 HE2];
    HO=[HO1 HO2];
    c5=[c1 c2];
end

if bu==1
%     FF=[FE,c55,FO,EMG,CFF];
%     FF=[FE,FO/4,FM/2,HE,HO,HM];
    FF=[FE,FO/4,FM/2,CFF,HE,HO,HM];
elseif bu==2
%     FF=[FE,c55,FO,CFF];
    FF=[FE,FO/4,CFFO,HE,HO];
else
    FF=[FE,CFFO,HE];
end
len=size(FF,2);
[f,f0]=blstimekind(len,bu,mode,kl);
end

function [f,f0]=blstimekind(len,bu,mode,kl)
% %     f1=[6 6 6 6 2 6 5*ones(1,4)    6*ones(1,9)    8 18 2*ones(1,9)    8 8 8 8 8 8    2*ones(1,18) 4   5*ones(1,9)   2*ones(1,9) 2*ones(1,9)    2*ones(1,9)];
%     f1=[6 6 6 6 2 6 5*ones(1,4)    6*ones(1,8)    8 18 2*ones(1,8)    8 8 8 8 8 8    2*ones(1,16) 4   5*ones(1,8)   2*ones(1,8) 2*ones(1,8)    2*ones(1,8) 1];%     f1=[6 6 6 6 2 6 5*ones(1,4)    6*ones(1,9)    8 17 2*ones(1,9)    8 8 8 8 8 8    2*ones(1,18) 4   5*ones(1,9)   2*ones(1,5) 1 2*ones(1,3) 2*ones(1,5) 1 2*ones(1,3)   2*ones(1,9)];
% %         f1=[6 6 6 6 2 6 5*ones(1,4)    6*ones(1,8)    8 18 2*ones(1,8)    8 8 8 8 8 8    32 4   5*ones(1,8)   16 16    16 1];%     f1=[6 6 6 6 2 6 5*ones(1,4)    6*ones(1,9)    8 17 2*ones(1,9)    8 8 8 8 8 8    2*ones(1,18) 4   5*ones(1,9)   2*ones(1,5) 1 2*ones(1,3) 2*ones(1,5) 1 2*ones(1,3)   2*ones(1,9)];
% FE=[[6 6 6 6 2  6*ones(1,6)]   [6 5 5 5 5 5 5]   [8  3 3 3 3 3 3 2*ones(1,8)] 6*ones(1,6)];
FE=[[6 6 6 6 2 6*ones(1,6)]   [6 5 5 5 5 5 5 ]   [8 9 9 1*ones(1,16)] 6*ones(1,6)];
FO=[6/kl*ones(1,9)   2*ones(1,8)];%6/kl*ones(1,9) 
FM=[6/kl*ones(1,9)  2*ones(1,8)];
HE=2;HO=2;HM=2;
CFF=[3 3 3 6 6 6];
CFFO=[3 3 6 6 6];
c5=0;
ex=[1];
if bu~=1 && bu~=2
    bu=3;
end
if bu==1 && mode==1
    f1=[FE c5 FO FM CFF HE HO HM ex ];
elseif bu==2 && mode==1
    f1=[FE c5 FO CFFO HE HO  ex];
elseif bu==3 && mode==1
    f1=[FE c5 CFFO HE ex];
elseif bu==1 && mode==2
    f1=[FE FE c5 c5 FO FO FM CFF CFF HE HE HO HO HM  ex];
elseif bu==2 && mode==2
    f1=[FE FE c5 c5 FO FO CFFO CFFO HE HE HO HO ex];
elseif bu==3 && mode==2
    f1=[FE FE c5 c5 CFFO CFFO HE HE ex];
end

% % if bu==1
% %     f1=[[6 6 6 6 2 6*ones(1,9)]   [6 5 5 5 5 5 5]   [8 18 2*ones(1,8)]  [6*ones(1,9)]  [24]    [6*ones(1,9)  2  2*ones(1,8)]  [6*ones(1,9) 2*ones(1,8)]   [3 3 3 6 6 6 ]  [3 3 3] 1];
% % elseif bu==2
% %     f1=[[6 6 6 6 2 6*ones(1,9)]   [6 5 5 5 5 5 5]   [8 18 2*ones(1,8)]  [6*ones(1,9)]  [24]    [6*ones(1,9)  2  2*ones(1,8)]  [3 3 3 6 6 6 ]   [3 3 3] 1];
% % else
% %     f1=[[6 6 6 6 2 6*ones(1,9)]   [6 5 5 5 5 5 5]   [8 18 2*ones(1,8)]  [6*ones(1,9)]  [24]    [3 3 3 6 6 6 ]   [3 3 3] 1];
% % end
% if bu==1
%     f1=[[6 6 6 6 2 6*ones(1,8)]   [6 5 5 1 5 5 5]   [8 18 2*ones(1,8)]  [8 8 8 8 8 8]  [24]    [2*ones(1,8)  2  2*ones(1,8)]  [5*ones(1,8) 2*ones(1,8)]   [3 3 3 6 6 6] 1];%  CFF   3 3 3 6 6 6    5*ones(1,8)   2*ones(1,8)
% else
%     f1=[[6 6 6 6 2 6*ones(1,8)]   [6 5 5 1 5 5 5]   [8 18 2*ones(1,8)]  [8 8 8 8 8 8]  [24]    [2*ones(1,8)  2  2*ones(1,8)]  [3 3 3 6 6 6]  1];%  CFF   3 3 3 6 6 6    5*ones(1,8)   2*ones(1,8)
% end
f=[f1 5*ones(1,(len-sum(f1))/5)];
f0=[6:10];
end

function sh=chosefea(kl,hm,b)
sh=[];
for kk=1:length(hm)
sh=[sh (hm(kk)-1)*6/kl+b+1:hm(kk)*6/kl+b];
end
end