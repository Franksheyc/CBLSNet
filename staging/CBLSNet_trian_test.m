%%%%%%%%%%%%%%%使用全部特征，time 10s/part%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
warning off all;
format compact;
add='sleep-EDF';

Dir=['D:\E\sleep\open database\features\',add,'\cassette'];%cassette  telemetry
outdir=['D:\E\sleep\open database\result\',add];
% Dir=['F:\sleep\open database\features\',add];%cassette  telemetry
% outdir=['F:\sleep\open database\result\',add];
% Dir=['F:\sleep\open database\features\','SHHS1'];
%%%
%%%%EDF==20 edf20   EDF==78 edf78   EDF==0 其他数据集
tic
EDF=20;
epochs = 3;
n=5;
h=5;%bls enhance窗
Time_related= [8 9 10];
mode=2;% timing mode
con_st=0.75;
si=20;bs=[1.8 1.2 1 0.9 0.7];%timing features[1.8 1.2 1 0.9 0.8]
bn=8;%测试集中每个SS中取多少组做训练
C = 2^-24; s = .8;%the l2 regularization parameter and the shrinkage scale of the enhancement nodes  C = 2^-30
times=2;%第一次bls循环次数
filename1=[Dir,'\',add,'20.mat'];
filename2=[Dir,'\',add,'20_tt.mat'];
%mode 取1为单独信号   mode取2为全体信号
%bu 取0不用EMG与EOG   取1用全部使用  2为只用EEG+EOG
%mo==1 取第一通道脑电，mo==2取第二通道脑电
mode1=1;bu=2;mod=2;
softmax1=@(x)exp(x)./sum(exp(x)')';
% [FF,CFF,label,leng,fo]=datage(filename1,filename2,mode1,bu);
[FF,CFF,label,leng,fo]=datage(filename1,filename2,mode1,bu,mod);

% [FF,index,phis,fo]=fea_select(FF,label,fo,0.95,1);%mRMR

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

li=[];
for i=1:length(leng)
    li=[li 0.8*(1:leng(i))/leng(i)];
end
FF=[FF li'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
len=size(FF,2);
label(label==0)=5;
% [FF,CFF,label,leng]=banwake(FF,CFF,label,leng,0.1);
clear Lab
Lab=label2lab(label);
% label=single(label);
% Lab=single(Lab);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%找出一部分样本，作为测试集%%%%%%%%%%%%%%%%%%%%
if length(leng)>6000
    hd=max(6000,floor(length(leng)/2));
else
    hd=length(leng);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% hd=length(leng)-3;%测试copy6用
CCC0=0;%CCC用来储存历次结果，第一次为0，是为了方便第一次 multibls 的使用
% % CCC0=0.84;%针对shhs2
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

% % indict=crossvalind('Kfold',floor(hd/flod)*flod,flod);
% indict=crossvalind('Kfold',hd,flods);
indict=ones(1,hd);
remk=rem(1:hd,flods);
for i=1:flods
    lk=find(remk==i-1);
    indict(lk)=i;
end
indict=indict';

% flods=1;
% indict=zeros(1,length(leng));
% ih=randperm(length(leng));
% indict(ih(1:ceil(length(leng)*0.3)))=1;
time_datap=toc;
tic
for k=1:flods
    tic
    ss1=(indict == k);%test
    ss=~ss1;%train
    ss=find(ss==1);
    ss1=find(ss1==1);

    marktr=getmark(ss,K);
    train_x=FF(marktr,:);
    train_x1=CFF(marktr,:);
    train_y=Lab(marktr,:);
    train_yy=label(marktr,:);

    markte=getmark(ss1,K);
    test_x=FF(markte,:);
    test_x1=CFF(markte,:);
    test_y=Lab(markte,:);
    test_yy=label(markte,:);
    % %     [train_x,index,phis,fo1]=fea_select(train_x,train_yy,fo,0.95,1);%mRMR
    % %     test_x=test_x(:,index);
    % %     len1=size(train_x,2);

    % 训练集样本扩充
    % %    [train_x,train_y,train_yy,leng1,dimarks]=samplebalance(train_x,train_yy,len1,leng,ss);
    % %%%%%%%%%%%%%%%fea sel%%%%%%useless%%%%%
    % [train_x,index,phis,fo]=fea_select(train_x,train_yy,fo,0.9,1);%mRMR
    % test_x=test_x(:,index);
    % %%%%%%%%%%%%%%%%%%%%%%%%%%
    train_err=zeros(1,epochs);test_err=zeros(1,epochs);
    train_time=zeros(1,epochs);test_time=zeros(1,epochs);
    [train_x, test_x]=pre_zca(train_x,test_x);%白化处理  使训练更以收敛
    [train_x1, test_x1]=pre_zca(train_x1,test_x1);
    [~,leen]=size(train_x);
    wh=[];beta11=[];x_tr=zeros(length(train_x),5);x_te=zeros(length(test_x),5);
    %     load('bls_initial_para_30_192.mat')
    ht=1;%multibls 重复次数

    study_rate=0.15; xx=zeros(length(train_x),5); yy=zeros(length(train_x),1); cr1=ones(length(train_x),1);ce1=ones(length(test_x),1);
    train_y0=train_y;

    %     train_y=conlabel(train_y,train_yy,[5]);
    for i=1:length(Time_related)
        time_related=Time_related(i);
        [f,f1]=timecausal0(time_related,leng(ss),5);
        frex{i}=f;
        f1rex{i}=f1;
        [f,f1]=timecausal0(time_related,leng(ss1),5);
        feex{i}=f;
        f1eex{i}=f1;
    end

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
        [f,f0,nfp]=blstimekind(size(train_x,2),fo);%种类参数特征

        % %         train_y(:,1)=train_y(:,1)*3;
        % %          train_y(:,3)=train_y(:,3)*3;
        %          train_y=train_y0;
        if kk==1
            len0=leng(ss);
            for bo=1:length(ss)
                ko0(bo)=sum(len0(1:bo));
            end
            ko0=[0 ko0];
            [train_y,study_rate,si]=trainybanlance(kk,train_x,train_yy,train_y,study_rate,xx,yy,cr1,train_yy,con_st,si);
            %     [train_y00,study_rate,si]=trainybanlance1(kk,train_x,train_yy,train_y,study_rate,xx,yy,cr1,train_yy,con_st,si,ss,ko0);
            %     train_y00=[];
            %     for bo=1:length(ss)
            %         trainy=train_y(ko0(bo)+1:ko0(bo+1),:);
            %         trainx=train_x(ko0(bo)+1:ko0(bo+1),:);
            %         trainyy=train_yy(ko0(bo)+1:ko0(bo+1));
            %         xx0=xx(ko0(bo)+1:ko0(bo+1),:);
            %         yy0=yy(ko0(bo)+1:ko0(bo+1),:);
            %         cr10=cr1(ko0(bo)+1:ko0(bo+1));
            %         [trainy,study_rate,si]=trainybanlance(kk,trainx,trainyy,trainy,study_rate,xx0,yy0,cr10,trainyy,con_st,si);
            %         train_y00=[train_y00;trainy];
            %     end
            seeds=56;%256   acc 86.26  f1  0.7944    kappa  0.8224
            rng(seeds)
            %     train_y=train_y00;
        elseif kk==2 %|| kk==3
            %     [train_y,study_rate,si]=trainybanlance1(kk,train_x,train_yy,train_y,study_rate,xx,yy,cr1,y,con_st,si,ss,ko0);
            [train_y,study_rate,si]=trainybanlance(kk,train_x,train_yy,train_y,study_rate,xx,yy,cr1,y,con_st,si);
        end


        % if  kk==1
        % %     sy=hist(train_yy,5);
        % %     sy=sy./sum(sy);
        % %     if sy(1)<0.05 || sy(3)<0.08
        % %         bmark=1;
        % %     else
        % %         bmark=0;
        % %     end
        % %     if bmark
        %     [train_y,study_rate,si]=trainybanlance(kk,train_x,train_yy,train_y,study_rate,xx,yy,cr1,train_yy,con_st,si);
        % %     end
        % else
        % %     if bmark
        %     [train_y,study_rate,si]=trainybanlance(kk,train_x,train_yy,train_y,study_rate,xx,yy,cr1,y,con_st,si);
        % %     end
        % end
        [trainacc,traintime,xx,yy,par,con0,train_y,tn] = multi_scale_bls4(train_x,train_y,test_x,test_y,s,C,h,f,f0,par,[],1,N1,N2,fe,fe1,[],nfp,[]);
        time_tr_f(k,kk)=toc;
        tic
        leng0=lengo{ss1};
        for xi=1:length(leng0)
            K0(xi)=sum(leng0(1:xi));
        end
        K0=[0 K0];
        [testacc,testtime,x,y,par,~,~] = multi_scale_bls4(train_x,train_y,test_x,test_y,s,C,h,f,f0,par,[],2,N1,N2,fe,fe1,con0,nfp,K0);
        clear K0
        time_te_f(k,kk)=toc;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % if testacc<0.76 && kk==1
        %     pause
        % end
        subs=plottimes(test_yy,test_y,y,x,leng,ss1,[1,0],0,0);
        subjectresultf(k,kk).fbls=subs;
        clear subs
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % [trainacc,traintime,xx,yy,par,pare,con0] = multi_scale_bls3(train_x,train_y,test_x,test_y,s,C,h,f,f0,par,pare,1,N1,N2,fe,fe1,[]);
        % [testacc,testtime,x,y,par,pare,~] = multi_scale_bls3(train_x,train_y,test_x,test_y,s,C,h,f,f0,par,pare,2,N1,N2,fe,fe1,con0);

        %         if k==1 || k==2
        %             if kk==1
        %                 n0=0.3;
        %                 nn=[2.5 1 3 1 1];
        %                 [xx,yy]=markovxx_part(ss,leng1,yy,xx,nn,n0);
        %                 [x,y]=markovxx_part(ss1,leng1,y,x,nn,n0);
        %             else
        %                 n0=n0+0.1;
        %                 nn=[2 1 2 1 1];
        %                 [xx,yy]=markovxx_part(ss,leng,train_yy,xx,nn,n0);
        %                 [x,y]=markovxx_part(ss1,leng,test_yy,x,nn,n0);
        %             end
        %         end
        train_length(k,kk)=tn;
        trianingacc{k,kk}=trainacc;
        %储存参数，方便后续测试集使用
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
        %         xx1=xx;x1=x;

        %         if length(train_x)~=len1
        %             train_x=train_x(dimarks,:);
        %             train_yy=train_yy(dimarks);
        %             train_y=train_y(dimarks,:);
        %             xx=xx(dimarks,:);
        %             yy=yy(dimarks);
        %         end


        % %         if kk==2
        % %             xtr=time_causal(xx,5,2);
        % %             train_x1=[train_x1 xtr];
        % %             xte=time_causal(x,5,2);
        % %             test_x1=[test_x1 xte];
        % %
        % %             %         else
        % %             %             xtr=time_causal(x_tr,5,2);
        % %             %             train_x=[train_x1 xtr];
        % %             %             xte=time_causal(x_te,5,2);
        % %             %             test_x=[test_x1 xte];
        % % %         end
        % % %         [TrainingAccuracy,TestingAccuracy,Trainingtime,Testingtime,x1,xx1,y1,yy1,wh,beta11] = multibls(train_x1,train_y,test_x1,test_y,s,C,ml,ht,[],[],[3 3 3 6 6 6 5*ones(1,size(xte,2)/5)],[3 7]);% C41 6 5 5 5 5
        % % %         if kk==1
        % %             fk=[3 3 3 6 6 6 5*ones(1,size(xtr,2)/5)];
        % %             fk0=[3 fk(7:end)];
        % %             [trainacc,traintime,xx1,yy1,part,con0] = multi_scale_bls4(train_x1,train_y,test_x1,test_y,s,C,h,fk,fk0,[],[],1,N1,N2,fe,fe1,[],nfp);
        % %             [testacc,testtime,x1,y1,~,~] = multi_scale_bls4(train_x1,train_y,test_x1,test_y,s,C,h,fk,fk0,part,[],2,N1,N2,fe,fe1,con0,nfp);
        % %             if length(find(y1~=y))/length(y)<0.2
        % %
        % %
        % %                 xx=(xx1+xx)/2;x=(x1+x)/2;
        % %                 [~,yy]=max(xx');yy=yy';
        % %                 [~,y]=max(x');y=y';
        % %             end
        % %         end


        % %         [~,ACC_pre1,accuracy_pre1,Precidion_pre1,Recall_pre1,F1SCORE_pre1,Kappa_pre1,mcc_pre1]=resultstatsic(y,test_y,n,k);%CCC
        % %         ACC_prec{k,kk}=ACC_pre1;
        % %         accuracy_prec{k,kk}=accuracy_pre1;
        % %         Precidion_prec{k,kk}=Precidion_pre1;
        % %         Recall_prec{k,kk}=Recall_pre1;
        % %         F1SCORE_prec{k,kk}=mean(F1SCORE_pre1);
        % %         Kappa_prec{k,kk}=Kappa_pre1;
        % %         Mcc_prec{k,kk}=mcc_pre1;
        % %

        %         [acc_stay_pre1,Acc_stay_pre1,accuracy_stay_pre1,F1score_stay_pre1,kappa_stay_pre1,ratio_stay_pre1,acc_trans_pre1,Acc1_trans_pre1,accuracy_trans_pre1,F1score_trans_pre1,kappa_trans_pre1,ratio_trans_pre1]=stage_transformation(y,test_y,test_yy);
        %         acc_stay_pre(k,kk)=acc_stay_pre1;
        %         Acc_stay_pre{k,kk}=Acc_stay_pre1;
        %         accuracy_stay_pre{k,kk}=accuracy_stay_pre1;
        %         F1score_stay_pre{k,kk}=F1score_stay_pre1;
        %         kappa_stay_pre{k,kk}=kappa_stay_pre1;
        %         ratio_stay_pre{k,kk}=ratio_stay_pre1;
        %         acc_trans_pre(k,kk)=acc_trans_pre1;
        %         Acc1_trans_pre{k,kk}=Acc1_trans_pre1;
        %         accuracy_trans_pre{k,kk}=accuracy_trans_pre1;
        %         F1score_trans_pre{k,kk}=F1score_trans_pre1;
        %         kappa_trans_pre{k,kk}=kappa_trans_pre1;
        %         ratio_trans_pre{k,kk}=ratio_trans_pre1;

        Par{k,kk}=par;
        Pare{k,kk}=pare;
        if kk>1  %%%防过拟合   这里考虑是否直接跳出
            %             if TestingAccuracy<test_acc(k,kk-1)
            if mean(F1SCORE_pre1) < mean(F1SCORE_CCC1)
                x=x0;xx=xx0;y=y0;yy=yy0;
            end
        end
        x0=x;
        xx0=xx;
        y0=y;
        yy0=yy;
        %%%%
        %         train_acc(k,kk)=TrainingAccuracy;
        test_acc(k,kk)=testacc
        %         train_time(k,kk)=Trainingtime;
        %         test_time(k,kk)=Testingtime;
        % %         [ACC_all_1,ACC_part_1,Accuracy_1,Precidion_1,Recall_1,F1SCORE_1,Kappa_1]=resultstatsic(y,test_y,n,k);
        %%%%%%%%%%%%%%%%%%%%%%第二次bls%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        lent=length(Time_related);X_tr=zeros(size(train_x,1),5);X_te=zeros(size(test_x,1),5);

        %         for i=1:5
        %             a(i)=length(find(yy==i))/length(yy);
        %             b(i)=length(find(y==i))/length(y);
        %         end
        %         varmark=mean(abs(a-b));
        %         if varmark>0.1 %|| std(b)>
        %             con_st=con_st+max(varmark,0.15);
        %             if con_st>0.95
        %                 con_st=0.95;
        %             end
        %         end


        for hi=1:lent
            tic
            ffr=frex{hi};ffe=feex{hi};
            ff1r=f1rex{hi};ff1e=f1eex{hi};
            time_related=Time_related(hi);%-(kk-1);%每次逐步减小区间
            [N1,N2,fe,fe1]=trainpar(len,2,kk,0,time_related);
            [C_train_x,C_test_x,cr1,ce1,ftr,fte,ft,ft0]=time_feas(xx,x,y,yy,time_related,ss,ss1,leng,mode,train_yy,test_y,si,bs,ffr,ff1r,ffe,ff1e);
            %             C_train_x=[C_train_x train_x(:,end)];
            %             C_test_x=[C_test_x test_x(:,end)];

            [C_train_x, C_test_x]=pre_zca(C_train_x,C_test_x);
            %             lyy=latency(train_yy,ss,leng,5);
            %             ly=latency(test_yy,ss1,leng,5);
            %             C_train_x=[C_train_x lyy];
            %             C_test_x=[C_test_x ly];


            %                 if kk==1
            par_2=[];
            %                 end
            %                 elseif kk==3
            %                     par_2=Par2{k-2,kk,hi};
            %                 else
            %                     par_2=Par2{k,kk-1,hi};
            %                 end
            % %               [trainacc,traintime,x_tr,y_tr,par_2,con0,~,~] = multi_scale_bls4(C_train_x,train_y,C_test_x,test_y,s,C,h,ft,ft0,par_2,[],1,N1,N2,fe,fe1,[],nfp,[]);
            % %               leng0=lengo{ss1};
            % %               for xi=1:length(leng0)
            % %                   K0(xi)=sum(leng0(1:xi));
            % %               end
            % %               K0=[0 K0];
            % %               [testacc,testtime,x_te,y_te,~,~,~] = multi_scale_bls4(C_train_x,train_y,C_test_x,test_y,s,C,h,ft,ft0,par_2,[],2,N1,N2,fe,fe1,con0,nfp,K0);
            % %               clear K0
            %                 tic
            %                 [~,~,~,~,~,~,~,~,de_s,de_l,dr_s,dr_l,~,~,bil]=gain_time_set(C_test_x,C_train_x,y,yy,test_y,train_y,cr1,ce1,con_st);
            %                 ssrl(k)=length(dr_l)/length(yy);
            %                 ssel(k)=length(de_l)/length(y);
            %                 ssre(k)=length(dr_s)/length(yy);
            %                 sses(k)=length(de_s)/length(y);
            %                 smark(k)=abs(ssl(k)-sse(k));
            [~,y_tr,~,x_tr,par_2,~,tnl,tns,bil,de_s,de_l,dr_s,dr_l]=conf_rec_sepre(C_train_x,C_test_x,y,yy,test_y,train_y,cr1,ce1,s,C,h,1,ft,ft0,par_2,N1,N2,fe1,fe,con_st,leng0,1);
            gxs{k,kk}=de_s;%每次的分组情况
            gxl{k,kk}=de_l;
            gxxs{k,kk}=dr_s;
            gxxl{k,kk}=dr_l;
            %                 toc
            Bil(k,:)=bil;
            Bi(k,kk)=bil(3);
            %输入用train_y0结果有提升，但是F1下降较多，不推荐
            CCrt(k,kk,hi)=length(find(y_tr==train_yy))/length(y_tr);

            train_lengthl(k,kk)=tnl;
            train_lengths(k,kk)=tns;
            Par2{k,kk,hi}=par_2;
            %             cr0=CCR{k,kk};

            %             X_tr0{hi}=x_tr;
            %             X_te0{hi}=x_te;
            X_tr=X_tr+x_tr*time_related/sum(Time_related);

            time_tr_t(k,kk,hi)=toc;




            tic

            [y_te,~,x_te,~,~,~,~,~,~,~,~,~,~]=conf_rec_sepre(C_train_x,C_test_x,y,yy,test_y,train_y,cr1,ce1,s,C,h,1,ft,ft0,par_2,N1,N2,fe1,fe,con_st,leng0,2);
            CCCt(k,kk,hi)=length(find(y_te==test_yy))/length(y_te);
            X_te=X_te+x_te*time_related/sum(Time_related);
            if length(find(y~=y_te))/length(y)>0.2
                y_te=y;
            end
            time_te_t(k,kk,hi)=toc;
            clear C_test_x_e C_test_x_r C_train_x_r C_train_x_e lab_te_e lab_te_r lab_tr_e lab_tr_r
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        end
        %         [~,rk]=sort(CCrt(k,kk,:),'descend');
        %         X_tr=X_tr0{rk(2)}+X_tr0{rk(1)};

        %         con_st=0.75;
       tic
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
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % if length(find(y_te==test_yy))/length(y_te)-testacc<0.02 && kk==1
        %    pause
        % end
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
        %         corra=[0.08 0.12 0.08 0.12 0.3];%dierwei 0.15
        %         corra=[0.05 0.12 0.06 0.12 0.5];
        corra=ones(1,5);
        %         corra=1./[2 1 1.2 1.2 1.2];
        mo=1;%1 为每次把结果续在之前的矩阵之后    2 为每次把结果放到相同位置，无需每次替换
        %         x_tr=softmax1(x_tr')';
        %         x_te=softmax1(x_te')';
        [crr,cee]=crce(X_tr,X_te,0,[1 1 1 1 1]);
        [train_x,test_x,ten]=retrainset(kk,train_x,test_x,leen,X_tr,X_te,crr,cee,leng(ss),leng(ss1));
        % %         [y_te,y_tr,x_te,x_tr]=conf_rec(train_x,test_x,y_te,y_tr,test_y,train_y,cr1,ce1,s,C,CCC0,ht,1);
        % %         [train_x,test_x]=retrainset(kk,train_x,test_x,leen,x_tr,x_te,corra,cr1,ce1);
        %         train_x=Train_x;test_x=Test_x;
        %         clear Train_x
        %         clear Test_x

        %         x_tr0=x_tr;
        %         x_te0=x_te;
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
        %         [acc_stay_C1,Acc_stay_C1,accuracy_stay_C1,F1score_stay_C1,kappa_stay_C1,ratio_stay_C1,acc_trans_C1,Acc1_trans_C1,accuracy_trans_C1,F1score_trans_C1,kappa_trans_C1,ratio_trans_C1]=stage_transformation(y,test_y,test_yy);
        %         acc_stay_pre(k,kk)=acc_stay_C1;
        %         Acc_stay_C{k,kk}=Acc_stay_C1;
        %         accuracy_stay_C{k,kk}=accuracy_stay_C1;
        %         F1score_stay_C{k,kk}=F1score_stay_C1;
        %         kappa_stay_C{k,kk}=kappa_stay_C1;
        %         ratio_stay_C{k,kk}=ratio_stay_C1;
        %         acc_trans_C(k,kk)=acc_trans_C1;
        %         Acc1_trans_C{k,kk}=Acc1_trans_C1;
        %         accuracy_trans_C{k,kk}=accuracy_trans_C1;
        %         F1score_trans_C{k,kk}=F1score_trans_C1;
        %         kappa_trans_C{k,kk}=kappa_trans_C1;
        %         ratio_trans_C{k,kk}=ratio_trans_C1;


        %         sname=['\epoch',num2str(k),'turn',num2str(kk)];
        %         save(sname,'xx','x','yy','y','test_yy','train_yy','x_tr','y_tr','x_te','y_te','ss','ss1','Mark','hd','indict');
        X{k,kk}=x;
        Xte{k,kk}=x_te;
        Y{k,kk}=y;
        Yte{k,kk}=y_te;
        Testy{k}=test_yy;
        sss{k}=ss;
        sss1{k}=ss1;
        time_reset(k,kk)=toc;
    end
end
runningtiem=toc
sublabel=[];suby=[];
for i=1:length(subjectresultf)
    sublabel1=subjectresultf(i,3).fbls.label;
    suby1=subjectresultf(i,3).fbls.result;
    sublabel=[sublabel;sublabel1];
    suby=[suby;suby1];
end
acc=length(find(sublabel==suby))/length(suby)
[~,confusion_matrix,accuracy,precision,recall,F1score,kappa,~] = Evaluation(suby,label2lab(sublabel),5)
mean(F1score)
% disp(add)
% f1=zeros(1,kk);f11=zeros(1,kk);f100=zeros(1,kk);f10=zeros(1,kk);
% for j=1:k
%     for i=1:kk
%         f11(i)=F1SCORE_pre{j,i};
%         f10(i,:)=F1score_pre{j,i};
%     end
%     f1=f1+f11;
%     f100=f100+f10;
% end
% f1=f1/k;
% f100=f100/k;
% [~,fm]=max(f1);
% if size(test_acc,1)~=1
%     accre=mean(test_acc);
% else
%     accre=test_acc;
% end
% [~,acm]=max(accre);
% if fm~=acm
%    ac1=f1*0.7+0.3*accre;
%    [~,acm]=max(ac1);
% end
% accm=acm;
%
% ACC=zeros(n);
% acc=zeros(1,n);kappa=zeros(1);
% for i=1:k
%     AC1=ACC_pre{i,accm};
%     ac1=accuracy_pre{i,accm};
%     ka1=Kappa_pre{i,accm};
%     ACC=ACC+AC1;
%     acc=acc+ac1;
%     kappa=kappa+ka1;
% end
% ac=accre(accm)
% mF1=f1(accm)
% F1=f100(accm,:)
% ACC=ACC/k
% acc=acc/k
% kappa=kappa/k
%
% [~,miniter]=min(test_acc(:,1));
% % subjectresult(fit(miniter/2)).fbls;%在结果中使用
% % subjectresult(fit(miniter/2)).tbls;
% plottimes(Testy{miniter},label2lab(Testy{miniter}),Y{miniter,1},X{miniter,1},leng,sss1{miniter},[1 0],0,1)

% savename=[outdir,'_with_bal_m',num2str(meth),'.mat'];%meth为使用方法类型，可在程序中设置
savename=[outdir,'subjectwise_',num2str(mode1),'EEG',num2str(mode1),'EOG',num2str(flods),'_fold_cv','.mat'];%meth为使用方法类型，可在程序中设置
save(savename,'ACC_part_CCC','accuracy_CCC','F1SCORE_CCC','F1score_CCC','Kappa_CCC','Mcc_CCC','Precidion_CCC','Recall_CCC','ACC_pre','accuracy_pre','F1SCORE_pre','F1score_pre','Kappa_pre','Mcc_pre','Precidion_pre','Recall_pre','CCC','test_acc',...
    'Xte','X','Yte','Y','sss','sss1','leng','Testy','subjectresultt','subjectresultf','seeds','Par','Par2','gxxs','gxxl','gxs','gxl');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%结果统计%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % end

%%%%%%%%%%%%%%%%%%%%%%用于全正特征
function spq=rate_entropy(XX)
b=size(XX,2);
for i=1:b
    p(:,i)=XX(:,i)./sum(XX,2);
end
% q=1-p;
% spq=(p./q).*log(p./q);
spq=p.*log(p);
end
function [c,m0]=find_error(y,test_yy,x)
% test_yy(test_yy==0)=5;
m0=find(y~=test_yy);
k=0;
for i=1:5
    for j=i+1:5
        k=k+1;
        c(:,k)=abs(x(:,i)./x(:,j));
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%用混淆矩阵（lab）扩展  构造时间序列

%%%%%%%%%%%%%%label2lab

%%%%%%%%%%%%%%%%%%%%%%%%%%%%事件序列结果

%%%%%%%%%%%%%%%%%%%找到m

% % % % % function lo_ri=loss_right(y,test_y)
% % % % %
% % % % % find(y~=testy);
% % % % % lo_ri=length
% % % % % end
% % % % % function lo_con=loss_continuous(y,test_yy)
% % % % % find(y~=test_yy);
% % % % % lo_con=length
% % % % % end
%%%%%%%%%%%%%%%%%%%%%%%%产生时序训练特征

%%%%%%%%%%%%用规则修改结果
function yy_2=changerule(yy_2,x_2)
p0=length(yy_2);
[~,b0]=max(x_2(:,1:4)');
for i=1:p0
    if x_2(i,5)>sum(x_2(i,1:4))
        yy_2(i)=5;
    else
        yy_2(i)=b0(i);
    end
end
end

%%%%%%%%%这个程序针对只对置信度较小的样本进行重新训练%%%%%%%%%%%%%%%%
% % % % function [x_2,yy_2]=ranky(y_s1,de_s,x_s1,y_l0,x,de_l)
% % % %     y_s1=[y_s1 de_s x_s1];%%%%%%%%%%%%%%第二次的y_s
% % % %     y_l1=[y_l0 x(de_l,:)];
% % % %     yy_2=sortrows([y_l1;y_s1],2);%%有对高置信度操作时，y_l0改为y_l1
% % % %     x_2=yy_2(:,3:7);
% % % %     yy_2=yy_2(:,1);
% % % % %     yy_21=sortrows([y_l0;y_s1(:,1:2)],2);
% % % % %     yy_21=yy_21(:,1);
% % % % end






function [ACC_all_1,ACC_part_1,Accuracy_1,Precidion_1,Recall_1,F1SCORE_1,Kappa_1,Mcc_1]=resultstatsic(y,test_y,n,j)
[Acc_all_1,Acc_part_1,accuracy_1,precision_1,recall_1,F1score_1,kappa_1,mcc_1] = Evaluation(y,test_y,n);
ACC_all_1=Acc_all_1;%ACC是总体百分比
ACC_part_1=Acc_part_1;%ACC每一项是百分比
Accuracy_1=accuracy_1;
Precidion_1=precision_1;
Recall_1=recall_1;
F1SCORE_1=F1score_1;
Kappa_1=kappa_1;
Mcc_1=mcc_1;
end

function fe=othersignal(da)
da1=da(da<quantile(da,0.9));
me1=mean(da1);
me=mean(da);
if me/me1>6
    da=da1;
    da(da>quantile(da,0.9))=quantile(da,0.9);
end
da=da-mean(da);
s=diff(da>0);
l1=length(find(s==1));
l2=length(find(s==-1));
zc=l1+l2;
st=std(da1);
mea=median(da1);
fe=[me mea st zc];
end



% function [y_te,y_tr,x_te,x_tr]=conf_rec(C_train_x,C_test_x,y,yy,test_y,train_y,cr1,ce1,s,C,CCC0,ht,tm,test_yy,f,f0)
% [C_test_x_e,C_test_x_r,C_train_x_r,C_train_x_e,lab_te_e,lab_te_r,lab_tr_e,lab_tr_r,de_s,de_l,dr_s,dr_l]=gain_time_set(C_test_x,C_train_x,y,yy,test_y,train_y,cr1,ce1,0.75);
% if tm==1
%     f1=f(1:5);
%     [~,Test_accuracy_l,~,~,x_l,xx_l,y_l,yy_l] = multibls_t(C_train_x_r(:,1:sum(f1)),lab_tr_r,C_test_x_r(:,1:sum(f1)),lab_te_r,s,C,9,CCC0,ht,f1,f0);
% else
%     [~,Test_accuracy_l,~,~,x_l,xx_l,y_l,yy_l] = multibls_t(C_train_x_r,lab_tr_r,C_test_x_r,lab_te_r,s,C,9,CCC0,ht,f,f0);
% end
% %     lab_tr_e(:,1)=lab_tr_e(:,1)*2;
% %     lab_tr_e(:,3)=lab_tr_e(:,3)*2;
%     [~,Test_accuracy_s,~,~,x_s,xx_s,y_s,yy_s] = multibls_t(C_train_x_e,lab_tr_e,C_test_x_e,lab_te_e,s,C,9,CCC0,ht,f,f0);%1:50 为C_xe
% %%%%%%%%%%%%%按被试排序
% y_te=rank_y(y_s,de_s,y_l,de_l);
% y_tr=rank_y(yy_s,dr_s,yy_l,dr_l);
% x_te=rank_y(x_s,de_s,x_l,de_l);
% x_tr=rank_y(xx_s,dr_s,xx_l,dr_l);
% length(find(y_te==test_yy))/length(y_te)
% end

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
% for i=1:ls
%     s0=[s0 s(i)*ones(1,5)];
% end
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
f=[f1 5*ones(1,(len-sum(f1))/5)];%增加的 c
% if length(f)>len
%     f0=[6:10 length(f1)+1];
% else
f0=[6:14 64:length(f1)];%84
% end
nfp=(len-sum(f1))/5;
end


function [Label,leng,X]=increasesstages(label,stage,x,mark,leng,si)
len1=[];X=[];Label=[];l1=1;l2=0;
for i=1:length(si)
    le=leng(si(i));
    l2=l2+le;
    label1=label(l1:l2);
    x1=x(l1:l2,:);
    [label1,h0,len,s]=increasedata(label1,stage);
    if ~isempty(s)
        x1=fixdata(x1,h0,s,mark);
    end
    Label=[Label;label1];
    len1=[len1 len];
    X=[X;x1];
    clear x1 len label1
    l1=l1+le;
end
leng(si)=len1;
end

function [train_x,train_y,train_yy,leng1,dimarks]=samplebalance(train_x,train_yy,len1,leng,ss)
mark=setdiff(1:size(train_x,2),[80:100]);
%     mark=1:size(train_x,2);
train_x=[train_x zeros(len1,1)];
leng1=leng;
m1=[1 3 5];
m2=[0.05 0.15 0.15];
k=0;
for i=1:length(m1)
    if length(find(train_yy==m1(i)))/length(train_yy)<m2(i)
        k=k+1;
        mk(k)=m1(i);
    end
end
if ~isempty(mk)
    for i=1:length(mk)
        [train_yy,leng1,train_x]=increasesstages(train_yy,mk(i),train_x,mark,leng1,ss);
    end
end
dimark=find(train_x(:,end)==1);
dimarks=setdiff(1:length(train_x),dimark);
train_x=train_x(:,1:end-1);
train_y=label2lab(train_yy);
end

function [X1,X2,LABEL,leng1]=banwake(X,X0,Label,leng,b)
l1=1;l2=0;X1=[];X2=[];LABEL=[];leng1=[];
for i=1:length(leng)
    le=leng(i);
    l2=l2+le;
    x=X(l1:l2,:);
    x1=X0(l1:l2,:);
    label=Label(l1:l2);
    n1=length(find(label==2))/length(label);
    n2=length(find(label==5))/length(label);
    if n2>n1 || n2>0.5
        m1=find(label~=5,1);
        m2=find(label~=5,1,'last');
        %         m=find(label==5);
        %         s=randperm(length(m));
        %         s1=s(1:floor(length(s)/5));
        s1=1:ceil(m1*b);
        s20=length(label)-m2;
        s2=floor(m2+s20*(1-b)):length(label);
        dimark=setdiff(1:size(x,1),union(s1,s2));
        x=x(dimark,:);
        x1=x1(dimark,:);
        label=label(dimark);
    end
    len1=length(label);
    X1=[X1;x];
    X2=[X2;x1];
    LABEL=[LABEL;label];
    leng1=[leng1;len1];
    l1=l1+le;
end
end

function [train_y00,study_rate,si]=trainybanlance1(kk,train_x,train_yy,train_y,study_rate,xx,yy,cr1,train_y0,con_st,si,ss,ko0)
train_y00=[];
for bo=1:length(ss)
    trainy=train_y(ko0(bo)+1:ko0(bo+1),:);
    trainx=train_x(ko0(bo)+1:ko0(bo+1),:);
    trainyy=train_yy(ko0(bo)+1:ko0(bo+1));
    if size(train_y0,1)==size(train_yy,1)
        trainy0=train_y0(ko0(bo)+1:ko0(bo+1));
    else
        trainy0=train_y0;
    end
    xx0=xx(ko0(bo)+1:ko0(bo+1),:);
    yy0=yy(ko0(bo)+1:ko0(bo+1),:);
    cr10=cr1(ko0(bo)+1:ko0(bo+1));
    [trainy,study_rate,si]=trainybanlance(kk,trainx,trainyy,trainy,study_rate,xx0,yy0,cr10,trainy0,con_st,si);
    train_y00=[train_y00;trainy];
end
end