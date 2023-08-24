clear
warning off all;
format compact;
add='sleep-EDF';

Dir=['D:\database\stage scoring\open database\features\',add,'\cassette'];%cassette  telemetry
outdir=['D:\database\stage scoring\open database\result\',add];
% Dir=['F:\sleep\open database\features\',add];%cassette  telemetry
% outdir=['F:\sleep\open database\result\',add];
% Dir=['F:\sleep\open database\features\','SHHS1'];
%%%
%%%%EDF==20 edf20   EDF==78 edf78   EDF==0 others

EDF=20;
epochs = 3;
n=5;h=5;
Time_related= [8 9 10];
mode=2;% timing mode
con_st=0.75;
C = 2^-24; s = .8;%the l2 regularization parameter and the shrinkage scale of the enhancement nodes  C = 2^-30
filename1=[Dir,'\',add,'20.mat'];
filename2=[Dir,'\',add,'20_tt.mat'];
% softmax1=@(x)exp(x)./sum(exp(x)')';
load(filename1)
load(filename2)
result1=[];result=[];test_lab=[];test_label=[];

if strcmp(add,'sleep-EDF')
    if EDF==20
        leng=leng(1:39);
        d=1:sum(leng(1:39));
        FF=FF(d,:);
        CFF=CFF(d,:);
        lable=label(d);
    end
    [leng,lengo]=sleepEDFleng(leng);
end
len=size(FF,2);
label(label==0)=5;
clear Lab
Lab=label2lab(label);

if length(leng)>6000
    hd=max(6000,floor(length(leng)/2));
else
    hd=length(leng);
end
CCC0=0;
Mark=randperm(hd)+1;
for i=1:hd
    K(i)=sum(leng(1:i));
end

if EDF==20
    flods=20;
elseif EDF==78
    flods=10;
elseif EDF==0 && leng<1000
    flods=5;
elseif EDF==0 && leng>1000
    flods=2;
end

indict=ones(1,hd);
remk=rem(1:hd,flods);
for i=1:flods
    lk=find(remk==i-1);
    indict(lk)=i;
end
indict=indict';

for k=1:flods
    ss1=(indict == k);%test
    ss=~ss1;%train
    ss=find(ss==1);
    ss1=find(ss1==1);

    train_x=getpartdata(FF,ss,K);
    train_x1=getpartdata(CFF,ss,K);
    train_y=getpartdata(Lab,ss,K);
    train_yy=getpartdata(label,ss,K);
    test_x=getpartdata(FF,ss1,K);
    test_x1=getpartdata(CFF,ss1,K);
    test_y=getpartdata(Lab,ss1,K);
    test_yy=getpartdata(label,ss1,K);

    train_err=zeros(1,epochs);test_err=zeros(1,epochs);
    train_time=zeros(1,epochs);test_time=zeros(1,epochs);
    [train_x, test_x]=pre_zca(train_x,test_x);
    [train_x1, test_x1]=pre_zca(train_x1,test_x1);
    [~,leen]=size(train_x);
    wh=[];beta11=[];x_tr=zeros(length(train_x),5);x_te=zeros(length(test_x),5);

    ht=1;

    study_rate=0.15; xx=zeros(length(train_x),5); yy=zeros(length(train_x),1); cr1=ones(length(train_x),1);ce1=ones(length(test_x),1);
    train_y0=train_y;
    mode1=1;bu=2;
    [f,f0]=blstimekind1(len,bu,mode1);

    for kk=1:epochs
        if kk==1
            ten=0;

        end
        [N1,N2,fe,fe1]=trainpar(len,1,kk,ten,Time_related(1));
        if kk==1
            ml=CCC0;
        else
            ml=CCC(k,kk-1);
        end
        if kk<2
            par=[];pare=[];
        end
        [f,f0,nfp]=blstimekind(size(train_x,2),fo);

        if kk==1
            len0=leng(ss);
            for bo=1:length(ss)
                ko0(bo)=sum(len0(1:bo));
            end
            ko0=[0 ko0];
            [train_y,study_rate,si]=trainybanlance(kk,train_x,train_yy,train_y,study_rate,xx,yy,cr1,train_yy,con_st,si);
            seeds=256;
            rng(seeds)
        elseif kk==2 %|| kk==3
            [train_y,study_rate,si]=trainybanlance(kk,train_x,train_yy,train_y,study_rate,xx,yy,cr1,y,con_st,si);
        end


        [trainacc,traintime,xx,yy,par,con0,train_y,tn] = multi_scale_bls4(train_x,train_y,test_x,test_y,s,C,h,f,f0,par,[],1,N1,N2,fe,fe1,[],nfp,[]);
        leng0=lengo{ss1};
        for xi=1:length(leng0)
            K0(xi)=sum(leng0(1:xi));
        end
        K0=[0 K0];
        [testacc,testtime,x,y,par,~,~] = multi_scale_bls4(train_x,train_y,test_x,test_y,s,C,h,f,f0,par,[],2,N1,N2,fe,fe1,con0,nfp,K0);
        clear K0
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % if testacc<0.76 && kk==1
        %     pause
        % end
        subs=plottimes(test_yy,test_y,y,x,leng,ss1,[1,0],0,0);
        subjectresultf(k,kk).fbls=subs;
        clear subs

        train_length(k,kk)=tn;
        trianingacc{k,kk}=trainacc;
        [ACCo_pre1,ACC_pre1,accuracy_pre1,Precidion_pre1,Recall_pre1,F1SCORE_pre1,Kappa_pre1,mcc_pre1]=resultstatsic(y,test_y,n,k);%CCC
        ACCo_pre{k,kk}=ACCo_pre1;
        ACC_pre{k,kk}=ACC_pre1;
        accuracy_pre{k,kk}=accuracy_pre1;
        Precidion_pre{k,kk}=Precidion_pre1;
        Recall_pre{k,kk}=Recall_pre1;
        F1SCORE_pre{k,kk}=mean(F1SCORE_pre1);
        F1score_pre{k,kk}=F1SCORE_pre1;
        Kappa_pre{k,kk}=Kappa_pre1;
        Mcc_pre{k,kk}=mcc_pre1;


        Par{k,kk}=par;
        Pare{k,kk}=pare;
        if kk>1
            if mean(F1SCORE_pre1) < mean(F1SCORE_CCC1)
                x=x0;xx=xx0;y=y0;yy=yy0;
            end
        end
        x0=x;
        xx0=xx;
        y0=y;
        yy0=yy;
        test_acc(k,kk)=testacc
        si=20;bs=[1.8 1.2 1 0.9 0.7];
        lent=length(Time_related);X_tr=zeros(size(train_x,1),5);X_te=zeros(size(test_x,1),5);

        for hi=1:lent
            time_related=Time_related(hi);
            [N1,N2,fe,fe1]=trainpar(len,2,kk,0,time_related);
            [C_train_x,C_test_x,cr1,ce1,ftr,fte,ft,ft0]=time_feas(xx,x,y,yy,time_related,ss,ss1,leng,mode,train_yy,test_y,si,bs);
            [C_train_x, C_test_x]=pre_zca(C_train_x,C_test_x);

            par_2=[];
            [y_te,y_tr,x_te,x_tr,par_2,~,tnl,tns,bil]=conf_rec(C_train_x,C_test_x,y,yy,test_y,train_y,cr1,ce1,s,C,h,1,ft,ft0,par_2,N1,N2,fe1,fe,con_st,leng0);
            Bil(k,:)=bil;
            Bi(k,kk)=bil(3);
            CCrt(k,kk,hi)=length(find(y_tr==train_yy))/length(y_tr);
            CCCt(k,kk,hi)=length(find(y_te==test_yy))/length(y_te);
            train_lengthl(k,kk)=tnl;
            train_lengths(k,kk)=tns;
            clear C_test_x_e C_test_x_r C_train_x_r C_train_x_e lab_te_e lab_te_r lab_tr_e lab_tr_r
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


            Par2{k,kk,hi}=par_2;
            if length(find(y~=y_te))/length(y)>0.2
                y_te=y;
            end
            X_tr=X_tr+x_tr*time_related/sum(Time_related);
            X_te=X_te+x_te*time_related/sum(Time_related);
        end
        X_tr=X_tr./sum(X_tr')';
        X_te=X_te./sum(X_te')';
        if kk>2
            X_tr=X_tr+X_tr-xx;
            X_te=X_te+X_te-x;
        end


        [~,y_tr]=max(X_tr');
        y_tr=y_tr';
        [~,y_te]=max(X_te');
        y_te=y_te';
        if length(find(y~=y_te))/length(y)>0.2
            y_te=y;
        end

        subs=plottimes(test_yy,test_y,y_te,x,leng,ss1,[1,0],0,0);
        subjectresultt(k,kk).tbls=subs;
        clear subs
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        CCr(k,kk)=length(find(y_tr==train_yy))/length(y_tr);
        CCC(k,kk)=length(find(y_te==test_yy))/length(y_te)

        for io=1:5
            rate_e(io)=length(find(y_te==io))/length(y_te);
            rate_r(io)=length(find(test_yy==io))/length(test_yy);
        end
        Rata_e{k,kk}=rate_e;
        Rata_r{k,kk}=rate_r;
        clear C_train_x C_test_x

        [crr,cee]=crce(X_tr,X_te,0,[1 1 1 1 1]);
        [train_x,test_x,ten]=retrainset(kk,train_x,test_x,leen,X_tr,X_te,crr,cee,leng(ss),leng(ss1));

        [ACCo_CCC1,ACC_part_CCC1,accuracy_CCC1,Precidion_CCC1,Recall_CCC1,F1SCORE_CCC1,Kappa_CCC1,Mcc_1]=resultstatsic(y_te,test_y,n,k);%CCC
        ACCo_CCC{k,kk}=ACCo_CCC1;
        ACC_part_CCC{k,kk}=ACC_part_CCC1;
        accuracy_CCC{k,kk}=accuracy_CCC1;
        Precidion_CCC{k,kk}=Precidion_CCC1;
        Recall_CCC{k,kk}=Recall_CCC1;
        F1SCORE_CCC{k,kk}=mean(F1SCORE_CCC1);
        F1score_CCC{k,kk}=F1SCORE_CCC1;
        Kappa_CCC{k,kk}=Kappa_CCC1;
        Mcc_CCC{k,kk}=Mcc_1;

        X{k,kk}=x;
        Xte{k,kk}=x_te;
        Y{k,kk}=y;
        Yte{k,kk}=y_te;
        Testy{k}=test_yy;
        sss{k}=ss;
        sss1{k}=ss1;
    end
end

disp(add)
f1=zeros(1,kk);f11=zeros(1,kk);f100=zeros(1,kk);f10=zeros(1,kk);
for j=1:k
    for i=1:kk
        f11(i)=F1SCORE_pre{j,i};
        f10(i,:)=F1score_pre{j,i};
    end
    f1=f1+f11;
    f100=f100+f10;
end
f1=f1/k;
f100=f100/k;
[~,fm]=max(f1);
if size(test_acc,1)~=1
    accre=mean(test_acc);
else
    accre=test_acc;
end
[~,acm]=max(accre);
if fm~=acm
    ac1=f1*0.7+0.3*accre;
    [~,acm]=max(ac1);
end
accm=acm;

ACC=zeros(n);
acc=zeros(1,n);kappa=zeros(1);
for i=1:k
    AC1=ACC_pre{i,accm};
    ac1=accuracy_pre{i,accm};
    ka1=Kappa_pre{i,accm};
    ACC=ACC+AC1;
    acc=acc+ac1;
    kappa=kappa+ka1;
end
ac=accre(accm)
mF1=f1(accm)
F1=f100(accm,:)
ACC=ACC/k
acc=acc/k
kappa=kappa/k

[~,miniter]=min(test_acc(:,1));
% subjectresult(fit(miniter/2)).fbls;
% subjectresult(fit(miniter/2)).tbls;
plottimes(Testy{miniter},label2lab(Testy{miniter}),Y{miniter,1},X{miniter,1},leng,sss1{miniter},[1 0],0,1)

% savename=[outdir,'_with_bal_m',num2str(meth),'.mat'];
savename=[outdir,'subjectwise_',num2str(mode1),'EEG',num2str(mode1),'EOG',num2str(flods),'_fold_cv','.mat'];
save(savename,'ACC_part_CCC','accuracy_CCC','F1SCORE_CCC','F1score_CCC','Kappa_CCC','Mcc_CCC','Precidion_CCC','Recall_CCC','ACC_pre','accuracy_pre','F1SCORE_pre','F1score_pre','Kappa_pre','Mcc_pre','Precidion_pre','Recall_pre','CCC','test_acc','Xte','X','Yte','Y','sss','sss1','leng','Testy','subjectresultt','subjectresultf','seeds');

function [ACC_all_1,ACC_part_1,Accuracy_1,Precidion_1,Recall_1,F1SCORE_1,Kappa_1,Mcc_1]=resultstatsic(y,test_y,n,j)
[Acc_all_1,Acc_part_1,accuracy_1,precision_1,recall_1,F1score_1,kappa_1,mcc_1] = Evaluation(y,test_y,n);
ACC_all_1=Acc_all_1;
ACC_part_1=Acc_part_1;
Accuracy_1=accuracy_1;
Precidion_1=precision_1;
Recall_1=recall_1;
F1SCORE_1=F1score_1;
Kappa_1=kappa_1;
Mcc_1=mcc_1;
end



function [train_x,test_x,endpoint]=retrainset(kk,train_x,test_x,leen,x_tr,x_te,cr1,ce1,ss,ss1)
tl=2;
mode=2;
endpoint=tl*mode*5+5;
% corra=[1 1 1 1 1];
a1=0.15;a2=a1-0.1;
s=1:-a1:1-a1*tl;
if mode==1
    s=sort(s);
else
    s1=sort(s)+a2;
    s=[s1(1:end-1) s];
end
ls=endpoint/5;
s0=[];

s0=repmat(s,1,5);
% g1=1;g2=0;
% for i=1:length(ss)
%     g2=g2+ss(i);
%     xtr(g1:g2,:)=time_causal(zscore((x_tr(g1:g2,:).*cr1(g1:g2))')',tl,mode);
%     g1=g1+ss(i);
% end
% g1=1;g2=0;
% for i=1:length(ss1)
%    g2=g2+ss1(i);
%    xte(g1:g2,:)=time_causal(zscore((x_te(g1:g2,:).*ce1(g1:g2))')',tl,mode);
%    g1=g1+ss1(i);
% end
% % [~,ytr]=max(train_x');
% % [~,yte]=max(test_x');
% % ytr=ytr';
% % yte=yte';
% % Cr=zeros(size(x_tr,1),5);
% % Ce=zeros(size(x_te,1),5);
% % for i=1:length(ytr)
% % Cr(i,ytr(i))=cr1(i);
% % end
% % for i=1:length(yte)
% % Ce(i,yte(i))=ce1(i);
% % end
xtr=time_causal(zscore((x_tr.*cr1)')',tl,mode);
xte=time_causal(zscore((x_te.*ce1)')',tl,mode);
%
% if kk==1
% [xtr, xte]=pre_zca(xtr,xte);
% end

xtr=xtr.*s0*1.5;
xte=xte.*s0*1.5;

% % r=[];
% % for i=1:5
% %     re=(0:4)*5+i;
% %     r=[r re];
% % end
% % xtr=xtr(:,r);
% % xte=xte(:,r);


% if kk>1
%    xtr=(xtr+train_x(:,leen+1:leen+endpoint))/2;
%    xte=(xte+test_x(:,leen+1:leen+endpoint))/2;
% end
% train_x0=[train_x xtr];
% test_x0=[test_x xte];
% [train_x0, test_x0]=pre_zca(train_x0,test_x0);
% xtr=train_x0(:,size(train_x,2)+1:end);
% xte=test_x0(:,size(test_x,2)+1:end);

% if mo==1
%     train_x=[train_x xtr];
%     test_x=[test_x xte];
% else
if kk==1
    train_x(:,leen+1:leen+endpoint)=xtr;
    test_x(:,leen+1:leen+endpoint)=xte;
else
    train_x(:,leen+1:leen+endpoint)=(train_x(:,leen+1:leen+endpoint)+xtr)/2;
    test_x(:,leen+1:leen+endpoint)=(test_x(:,leen+1:leen+endpoint)+xte)/2;
end



%         train_x(:,leen+1:leen+endpoint)=xtr;test_x(:,leen+1:leen+endpoint)=xte;
%         %     train_x=[train_x x_tr];test_x=[test_x x_te];
%     else
%         %     x=CCC(kk)/CCC(kk-1);
%         %     x=2^x;
%         xtr1=x_tr./corra;
%         xte1=x_te./corra;
%         train_x(:,leen+1:leen+endpoint+5)=[xtr xtr1];test_x(:,leen+1:leen+endpoint+5)=[xte xte1];
%         %     train_x=[train_x x_tr./corra];test_x=[test_x x_te./corra];
%     end
% end
end

function L=latency(Y,ss,leng,n)
len=leng(ss);
for i=1:length(ss)
    k(i)=sum(len(1:i));
end
k=[0 k];
L=[];
for j=1:length(ss)
    ly=[];
    y=Y(k(j)+1:k(j+1));
    for i=1:n
        a=find(y==i);
        b=[0;a];
        b=diff(b)/(length(b)-1);
        x0=[b a];
        ly=[ly;x0];
    end
    ly=sortrows(ly,2);
    ly=ly(:,1);
    L=[L;ly];
end

end

function [f,f0,nfp]=blstimekind(len,f1)
f=[f1 5*ones(1,(len-sum(f1))/5)];
% if length(f)>len
%     f0=[6:10 length(f1)+1];
% else
f0=[6:14 64:length(f1)];%84
% end
nfp=(len-sum(f1))/5;
end

function [f,f0]=blstimekind1(len,bu,mode)
% %     f1=[6 6 6 6 2 6 5*ones(1,4)    6*ones(1,9)    8 18 2*ones(1,9)    8 8 8 8 8 8    2*ones(1,18) 4   5*ones(1,9)   2*ones(1,9) 2*ones(1,9)    2*ones(1,9)];%基础train_x
%     f1=[6 6 6 6 2 6 5*ones(1,4)    6*ones(1,8)    8 18 2*ones(1,8)    8 8 8 8 8 8    2*ones(1,16) 4   5*ones(1,8)   2*ones(1,8) 2*ones(1,8)    2*ones(1,8) 1];%     f1=[6 6 6 6 2 6 5*ones(1,4)    6*ones(1,9)    8 17 2*ones(1,9)    8 8 8 8 8 8    2*ones(1,18) 4   5*ones(1,9)   2*ones(1,5) 1 2*ones(1,3) 2*ones(1,5) 1 2*ones(1,3)   2*ones(1,9)];%基础train_x
% %         f1=[6 6 6 6 2 6 5*ones(1,4)    6*ones(1,8)    8 18 2*ones(1,8)    8 8 8 8 8 8    32 4   5*ones(1,8)   16 16    16 1];%     f1=[6 6 6 6 2 6 5*ones(1,4)    6*ones(1,9)    8 17 2*ones(1,9)    8 8 8 8 8 8    2*ones(1,18) 4   5*ones(1,9)   2*ones(1,5) 1 2*ones(1,3) 2*ones(1,5) 1 2*ones(1,3)   2*ones(1,9)];%基础train_x
% FE=[[6 6 6 6 2  6*ones(1,6)]   [6 5 5 5 5 5 5]   [8  3 3 3 3 3 3 2*ones(1,8)] 6*ones(1,6)];
FE=[[6 6 6 6 2 6*ones(1,6)]   [6 5 5 5 5 5 5 ]   [8 9 9 1*ones(1,16)] 6*ones(1,6)];%  有效[8  9 9 1*ones(1,16)]
FO=[6/2*ones(1,9)   2*ones(1,8)];%6/kl*ones(1,9) 
FM=[6/2*ones(1,9)  2*ones(1,8)];
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