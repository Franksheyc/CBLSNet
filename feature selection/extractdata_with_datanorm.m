% clear
clc
t='ccshs';
% ad1=['F:\sleep\data\',t,'\feature\sleep cassette'];%sleep edf
% ad1=['F:\sleep\data\',t,'\feature\dod-h'];
ad1=['F:\sleep\data\',t,'\feature'];
% ad1=['F:\sleep\data\',t,'\feature\dod'];
% ad2=['F:\sleep\data\',t,'\label\dod'];
ad2=['F:\sleep\data\',t,'\label'];
ad=['F:\sleep\open database\features\',t,'\'];
mkdir(ad);
% ad='F:\sleep\data\sof\feature';
C3=[];C4=[];F4=[];O2=[];EMG=[];EOG=[];Label=[];leng=[];a=[];
Hjorth=[];C3t=[];C4t=[];Hjorth_activity=[];Hjorth_mobility=[];Hjorth_compelicity=[];Label=[];leng=[];H=[];
stc3=zeros(1,98);stc4=zeros(1,98);stemg=zeros(1,40);steog=zeros(1,38);
stc31=zeros(1,48);stc41=zeros(1,48);
% % C31=[];C41=[];EMG1=[];EOG1=[];C3t11=[];C4t11=[];
n1=dir(ad1);
% n2=dir(ad);
n2=dir(ad2);

for i=1:length(n1)-2
    name=[ad1,'\',n1(i+2).name];
    name1=[ad2,'\',n2(i+2).name];
%     name1=[ad2,'\',n2(i+2).name];
%     load(name1)
%     [c3,c4,f4,o2,emg,eog,label1,len1]=generate_data(name);
    load(name,'C30_tfs');
    load(name,'C40_tfs');
    load(name,'C3H');
    load(name,'C4H');
    load(name,'EOGLH');
    load(name,'EOGRH');
    load(name,'EMGH');
    H1=[C3H C4H EOGLH EOGRH EMGH];
    [c3,c4,emg,eog,label,~]=generate_data1(name);
%     load(name1);
    label1=label;
    clear label
    C3t1=eegtime(C30_tfs);
    C4t1=eegtime(C40_tfs);
 
%     a = findthenans(c3,c4,emg,eog,C3t1,C4t1,Hjorth_activity1,Hjorth_mobility1,Hjorth_compelicity1);
    a = findthenans(c3,c4,emg,eog,C3t1,C4t1);
    if length(a)>length(label1)/3
        continue
    else
        b=setdiff(1:length(label1),a);
        c3=c3(b,:);
        c4=c4(b,:);
        emg=emg(b,:);
        eog=eog(b,:);
        C3t1=C3t1(b,:);
        C4t1=C4t1(b,:);
        H1=H1(b,:);
%         Hjorth_activity1=Hjorth_activity1(b,:);
%         Hjorth_mobility1=Hjorth_mobility1(b,:);
%         Hjorth_compelicity1=Hjorth_compelicity1(b,:);
        label1=label1(b);
% %         len1=length(b);
%         [c3,c4,eog,emg]=datanorms(c3,c4,eog,emg);
    end
    
    label0=union(find(label1~=5),find(label1==-1));
    label1=label1(label0);
    len1=length(label1);
    c3=c3(label0,:);
    c4=c4(label0,:);
    eog=eog(label0,:);
    emg=emg(label0,:);
    C3t1=C3t1(label0,:);
    C4t1=C4t1(label0,:);
    
    if isempty(len1)
        continue
    end
    
    C3=[C3;c3];
    C4=[C4;c4];
%     F4=[F4;f4];
%     O2=[O2;o2];
    EMG=[EMG;emg];
    EOG=[EOG;eog];
    C3t=[C3t;C3t1];
    C4t=[C4t;C4t1];
    H=[H;H1];


    Label=[Label;label1];
    leng=[leng;len1];
    leng1(i)=len1;
    he(i)=length(leng)-i;

    c3m(i,:)=mean(c3);
    c3tm(i,:)=mean(C3t1);
    c4m(i,:)=mean(c4);
    c4tm(i,:)=mean(C4t1);
    emgm(i,:)=mean(emg);
    eogm(i,:)=mean(eog);
    c3s(i,:)=std(c3);
    c3ts(i,:)=std(C3t1);
    c4s(i,:)=std(c4);
    c4ts(i,:)=std(C4t1);
    emgs(i,:)=std(emg);
    eogs(i,:)=std(eog);
        clear c3 c4 emg eog c31 c41 emg1 eog1 label1 len1 h label0
    a=[];
end 
leng=leng(leng~=0);
label=Label;
clear Label

[C3,C4,C3t,C4t,EOG,EMG,H,leng,label]=dedata(C3,C4,C3t,C4t,EOG,EMG,H,leng,label,1);

[C31,C41,~,~,~]=datanorms(C3,C4,EOG,EMG,H);
softmax1=@(x)exp(x)./sum(exp(x));
c33=C3(:,27:50);
c44=C4(:,27:50);
for i=1:4
    c3(:,(i-1)*6+1:i*6)=mapminmax(c33(:,(i-1)*6+1:i*6)')';
    c4(:,(i-1)*6+1:i*6)=mapminmax(c44(:,(i-1)*6+1:i*6)')';
end
c3(isnan(c3))=0.0001;
c4(isnan(c4))=0.0001;
C3=[C3(:,1:26) C31(:,27:62) C3(:,51:end)];
C4=[C4(:,1:26) C41(:,27:62) C4(:,51:end)];
clear C31 C41

[C32,C42,EOG2,EMG2]=feature_enhance(C3,C4,EOG,EMG);
[C11,C12,C01,C02,M,CM1,CM2,C21,C22,di]=feature_en(C3,C4,C3t,C4t,EMG,C32,C42,EMG2);

[C3]=deleteoutliners(C3,0.1);
[C4]=deleteoutliners(C4,0.1);
[EMG]=deleteoutliners(EMG,0.1);
[EOG]=deleteoutliners(EOG,0.1);
[C3t]=deleteoutliners(C3t,0.1);
[C4t]=deleteoutliners(C4t,0.1);
[C32]=deleteoutliners(C32,0.1);
[C42]=deleteoutliners(C42,0.1);
[H]=deleteoutliners(H,0.1);

C3=[C3(:,1:26) c3 C3(:,63:end)];
C4=[C4(:,1:26) c4 C4(:,63:end)];

[C3,C4,EOG,EMG,H]=datanorms(C3,C4,EOG,EMG,H);
C3t=zscore(C3t);
C4t=zscore(C4t);

C3t=single(C3t);
C4t=single(C4t);
C3=single(C3);
C4=single(C4);
EOG=single(EOG);
EMG=single(EMG);
C32=single(C32);
C42=single(C42);
EOG2=single(EOG2);
EMG2=single(EMG2);
C11=single(C11);
C12=single(C12);
C01=single(C01);
C02=single(C02);
M=single(M);
CM1=single(CM1);
CM2=single(CM2);
C21=single(C21);
C22=single(C22);
H=single(H);


mode=1;bu=2;mo=2;
[FF,CFF,f,f0]=datage(C3,C4,c3,c4,EOG,EMG,...
    C3t,C4t,C32,C42,EMG2,EOG2,C11,C12,C01,C02,CM1,CM2,...
    C21,C22,H,mode,bu,mo);
save([ad,t,'lo.mat'],'FF','label','leng')
save([ad,t,'lo_tt.mat'],'CFF')
% save([ad,t,'lo.mat'],'C3','C4','c3','c4','EOG','EMG','label','leng')
% save([ad,t,'lo_tt.mat'],'C3t','C4t','C32','C42','EMG2','EOG2','C11','C12','C01','C02','M','CM1','CM2','C21','C22','H')

function f=eegtime(C3_tfs)
% for i=1:6
%     means(:,i)=C3_tfs(i).means;
%     stds(:,i)=C3_tfs(i).stds;
%     ku(:,i)=C3_tfs(i).ku;
%     sk(:,i)=C3_tfs(i).sk;
%     eng(:,i)=log(C3_tfs(i).Eng);
%     pe(:,i)=C3_tfs(i).Pe;
%     zc(:,i)=C3_tfs(i).Z;
%     median(:,i)=C3_tfs(i).median;
means=C3_tfs.means;
stds=C3_tfs.stds;
ku=C3_tfs.ku;
sk=C3_tfs.sk;
eng=log(C3_tfs.Eng);
pe=C3_tfs.Pe;
zc=C3_tfs.Z;
median=C3_tfs.median;
pe1=C3_tfs.Pe1;
lrssv=C3_tfs.LRSSV;
    f=[means stds ku sk eng pe zc median lrssv];
% end
end

function [Hactivity,Hmobility,Hcompelicity]=shjorth(H)
Hactivity=H(:,([1:6]-1)*3+1);
Hmobility=H(:,([1:6]-1)*3+2);
Hcompelicity=H(:,([1:6]-1)*3+3);
end

% function [label,len,c3,c4,emg,eog,c31,c41,emg1,eog1,eog2]=increasesstages(label,i,c3,c4,eog,emg,c31,c41,emg1,eog1,eog2)
function [label,len,c3,c4,emg,eog,c31,c41]=increasesstages(label,i,c3,c4,eog,emg,c31,c41)
    [label,h0,len,s]=increasedata(label,i);
    if ~isempty(s)
        c3=fixdata(c3,h0,s);
        c4=fixdata(c4,h0,s);
        emg=fixdata(emg,h0,s);
        eog=fixdata(eog,h0,s);
        c31=fixdata(c31,h0,s);
        c41=fixdata(c41,h0,s);
%         emg1=fixdata(emg1,h0,s);
%         eog1=fixdata(eog1,h0,s);
%         eog2=fixdata(eog2,h0,s);
    end
end