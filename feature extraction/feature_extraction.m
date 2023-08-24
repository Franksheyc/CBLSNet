clear
clc
warning off
% objDir1='D:\database\MASS\staged\data\SS';
t='dod';
objDir1=['F:\sleep\data\',t,'\data\ccshs'];
objDir2=['F:\sleep\data\',t,'\label\ccshs'];
% objDir1=['H:\pshhs\data\',t,'1'];
% objDir2=['F:\sleep\data\',t,'\label\shhs1'];
outputDir = ['F:\sleep\data\',t,'\feature\ccshs'];
mkdir(outputDir)
% Lna1=1;Lla1=1;
k=0;
% Lna2=0;Lla2=0;mark=0;
b=7;
samplingrate=100;

%     path1=[objDir1,num2str(ij)];
%     path2=[objDir2,num2str(ij)];
n1=dir(objDir1);
n2=dir(objDir2);
% n3=dir('F:\sleep\data\shhs\feature');
mm=ceil((length(n1)-2)/1);

% for i=1:length(n3)-3
% pp(i)=str2num(n3(i+2).name(end-9:end-4));
% end
mk=1;
% hm=mm*(mk-1)+1:min(mm*mk,length(n1)-2);
% for i=hm
% pp1(i)=str2num(n1(i+2).name(end-9:end-4));
% end
% h1=length(intersect(pp,pp1));
% mk=1;
for ii=mm*(mk-1)+1:min(mm*mk,length(n1)-2)%ceil((length(n1)-2)/6)*5+1:length(n1)-2%ceil((length(n1)-2)/6)*5%%ceil((length(n1)-2)/6)*6%1:ceil((length(n1)-2)/2)%
    na=[n1(ii+2).name];
    la=[n2(ii+2).name];
    filename1=[objDir1,'\',na];
    filename2=[objDir2,'\',la];
    load(filename1);
    load(filename2);
    
    unit=size(DAT,2)/samplingrate/5;
    
    mark=find(sum(DAT,2)~=0);
    DAT=DAT(mark,:);
    label=label(mark);
    
    aa=find(label==5);
    bb=setdiff(1:length(label),aa);
    label=label(bb);
    DAT=DAT(bb,:);
    
    len=size(DAT,1);
    lla=length(label);
    lab=zeros(len,5);
    k=k+1
    %%%
    for ijk=1:len
        switch label(ijk)
            case 0
                lab(ijk,5)=1;
            case 1
                lab(ijk,1)=1;
            case 2
                lab(ijk,2)=1;
            case 3
                lab(ijk,3)=1;
            case 4
                lab(ijk,4)=1;
        end
    end
    
    C3=DAT(:,1:unit*samplingrate);
    C4=DAT(:,unit*samplingrate+1:2*unit*samplingrate);
    EOGL=DAT(:,2*unit*samplingrate+1:3*unit*samplingrate);
    EOGR=DAT(:,3*unit*samplingrate+1:4*unit*samplingrate);
    EMG =DAT(:,4*unit*samplingrate+1:end);
    % % % % % %
    
    [C3_frec,C3coe]=wavelet_features(b,C3);
    [C4_frec,C4coe]=wavelet_features(b,C4);
    C3_frec=struct2cell(C3_frec);
    C4_frec=struct2cell(C4_frec);
    for i=1:length(C3_frec)
        C3_tfs=time_features(C3_frec{i},'EEG');
        C4_tfs=time_features(C4_frec{i},'EEG');
        C3tfs{i}=C3_tfs;
        C4tfs{i}=C4_tfs;
    end
    EOGcoff=EOG_cof(EOGL,EOGR,samplingrate);
    EOGLtfs=time_features(EOGL,'EOG');
    EOGRtfs=time_features(EOGR,'EOG');
    EMGtfs=time_features(EMG,'EMG');
    u=samplingrate*(unit/6);
    C30_tfs=time_features(C3,'EEG');
    C40_tfs=time_features(C4,'EEG');
    
    C3H=[];C4H=[];EOGLH=[];EOGRH=[];EMGH=[];
    for i=1:6
        C31=C3(:,u*(i-1)+1:u*i);
        [C3a,C3m,C3c]=Hjorth_features(C31);
        C41=C4(:,u*(i-1)+1:u*i);
        [C4a,C4m,C4c]=Hjorth_features(C41);
        EOGL1=EOGL(:,u*(i-1)+1:u*i);
        [EOGLa,EOGLm,EOGLc]=Hjorth_features(EOGL1);
        EOGR1=EOGR(:,u*(i-1)+1:u*i);
        [EOGRa,EOGRm,EOGRc]=Hjorth_features(EOGR1);
        EMG1=EMG(:,u*(i-1)+1:u*i);
        [EMGa,EMGm,EMGc]=Hjorth_features(EMG1);
        C3H=[C3H C3a,C3m,C3c];
        C4H=[C4H C4a,C4m,C4c];
        EOGLH=[EOGLH EOGLa,EOGLm,EOGLc];
        EOGRH=[EOGRH EOGRa,EOGRm,EOGRc];
        EMGH=[EMGH EMGa,EMGm,EMGc];
    end   
    savename=[outputDir,na];
    save(savename,'C3coe','C4coe','C3tfs','C4tfs','C30_tfs','C40_tfs','EOGLtfs','EOGRtfs','EOGcoff','EMGtfs','C3H','C4H','EOGLH','EOGRH','EMGH','label','lab');
    clear C3coe C4coe C3tfs C4tfs C3coe1 C4coe1 C3tfs1 C4tfs1 EOGLtfs EOGRtfs EOGcoff EMGtfs C3 C4 EOGL EOGR EMG C30_tfs C40_tfs C3H C4H EOGLH EOGRH EMGH
    leng(k)=len;
    len1=length(find(label~=5));
    leng_no9(k)=len1;
end