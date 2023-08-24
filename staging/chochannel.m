function [FE,FO,CFF,CFFO,HE,HO,c5]=chochannel(leng,FE1,FE2,FO1,FO2,CFF1,CFF2,CFFO1,CFFO2,HE1,HE2,HO1,HO2,c1,c2,mo)
%mo==1,只用C4  mo==2 只用C3   mo==3混用C3，C4
len=length(leng);
s=randperm(len);
s1=randperm(len);
m=ones(length(leng),2)*mo;
%%%以下为任选一通道
if mo==3
    m=ones(length(leng),2);
    m(s(1:ceil(length(s)/2)),1)=mo-1;
%     m(s1(1:ceil(length(s)/2)),2)=mo;
end
%%%%
l1=1;l2=0;
for i=1:len
    le=leng(i);
    l2=l2+le;
    eval(['FE(l1:l2,:)=FE',num2str(m(i,1)),'(l1:l2,:);']);
    eval(['FO(l1:l2,:)=FO',num2str(m(i,1)),'(l1:l2,:);']);
    eval(['CFF(l1:l2,:)=CFF',num2str(m(i,1)),'(l1:l2,:);']);
    eval(['CFFO(l1:l2,:)=CFFO',num2str(m(i,1)),'(l1:l2,:);']);
    eval(['HE(l1:l2,:)=HE',num2str(m(i,1)),'(l1:l2,:);']);
    eval(['HO(l1:l2,:)=HO',num2str(m(i,1)),'(l1:l2,:);']);
    eval(['c5(l1:l2,:)=c',num2str(m(i,1)),'(l1:l2,:);']);
    l1=l1+le;
end