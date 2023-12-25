function [y_te,y_tr,x_te,x_tr,PAR,train_y,tnl,tns,bil,de_s,de_l,dr_s,dr_l]=conf_rec_sepre(C_train_x,C_test_x,y,yy,test_y,train_y,cr1,ce1,s,C,h,tm,f,f0,PAR,N1,N2,fe,fe1,con_st,leng,mode)
[C_test_x_e,C_test_x_r,C_train_x_r,C_train_x_e,lab_te_e,lab_te_r,lab_tr_e,lab_tr_r,de_s,de_l,dr_s,dr_l,~,~,bil]=gain_time_set(C_test_x,C_train_x,y,yy,test_y,train_y,cr1,ce1,con_st);
[par_l,pare_l,par_s,pare_s]=repar(PAR);
[kl,ks]=partk(leng,de_l,de_s);
if tm==1
    f1=f;%(1:5);
    if mode==1
    [trainaccl,~,xx_l,yy_l,par_l,con0,lab_tr_r,tnl] = multi_scale_bls4(C_train_x_r(:,1:sum(f1)),lab_tr_r,C_test_x_r(:,1:sum(f1)),lab_te_r,s,C,h,f1,f0,par_l,pare_l,1,N1,N2,fe,fe1,[],0,[]);
    end
    if mode==2
        con0=[0 1];
    [~,~,x_l,y_l,~,~,~] = multi_scale_bls4(C_train_x_r(:,1:sum(f1)),lab_tr_r,C_test_x_r(:,1:sum(f1)),lab_te_r,s,C,h,f1,f0,par_l,pare_l,2,N1,N2,fe,fe1,con0,0,kl);
    end
%     [~,~,xx_l,yy_l,par_l,pare_l,con0] = multi_scale_bls3(C_train_x_r(:,1:sum(f1)),lab_tr_r,C_test_x_r(:,1:sum(f1)),lab_te_r,s,C,h,f1,f0,par_l,pare_l,1,N1,N2,fe,fe1,[]);
%     [~,~,x_l,y_l,~,~,~] = multi_scale_bls3(C_train_x_r(:,1:sum(f1)),lab_tr_r,C_test_x_r(:,1:sum(f1)),lab_te_r,s,C,h,f1,f0,par_l,pare_l,2,N1,N2,fe,fe1,con0);
%     [~,Test_accuracy_l,~,~,x_l,xx_l,y_l,yy_l] = multibls_t(C_train_x_r(:,1:sum(f1)),lab_tr_r,C_test_x_r(:,1:sum(f1)),lab_te_r,s,C,9,CCC0,ht,f1,f0);
else
%     lab_tr_r(:,1)=lab_tr_r(:,1)*2;
%     lab_tr_r(:,3)=lab_tr_r(:,3)*2;
 if mode==1
[trainaccl,~,xx_l,yy_l,par_l,con0,lab_tr_r,tnl] = multi_scale_bls4(C_train_x_r,lab_tr_r,C_test_x_r,lab_te_r,s,C,h,f,f0,par_l,pare_l,1,N1,N2,fe,fe1,[],0,[]);
 end
 if mode==2
     con0=[0 1];
[~,~,x_l,y_l,~,~,~] = multi_scale_bls4(C_train_x_r,lab_tr_r,C_test_x_r,lab_te_r,s,C,h,f,f0,par_l,pare_l,2,N1,N2,fe,fe1,con0,0,kl);
 end
%     [~,~,xx_l,yy_l,par_l,pare_l,con0] = multi_scale_bls3(C_train_x_r,lab_tr_r,C_test_x_r,lab_te_r,s,C,h,f,f0,par_l,pare_l,1,N1,N2,fe,fe1,[]);
%     [~,~,x_l,y_l,~,~,~] = multi_scale_bls3(C_train_x_r,lab_tr_r,C_test_x_r,lab_te_r,s,C,h,f,f0,par_l,pare_l,2,N1,N2,fe,fe1,con0);
%     [~,Test_accuracy_l,~,~,x_l,xx_l,y_l,yy_l] = multibls_t(C_train_x_r,lab_tr_r,C_test_x_r,lab_te_r,s,C,9,CCC0,ht,f,f0);
end
%     lab_tr_e(:,1)=lab_tr_e(:,1)*2;
%     lab_tr_e(:,3)=lab_tr_e(:,3)*2;
    
%     [~,lare]=max(lab_te_r');
%     length(find(lare==test_yy(de_l)))/length(de_l)
 if mode==1
[trainaccs,~,xx_s,yy_s,par_s,con0,lab_tr_e,tns] = multi_scale_bls4(C_train_x_e,lab_tr_e,C_test_x_e,lab_te_e,s,C,h,f,f0,par_s,pare_s,1,N1,N2,fe,fe1,[],0,[]);
 end
 if mode==2
     con0=[0 1];
[~,~,x_s,y_s,~,~,~] = multi_scale_bls4(C_train_x_e,lab_tr_e,C_test_x_e,lab_te_e,s,C,h,f,f0,par_s,pare_s,2,N1,N2,fe,fe1,con0,0,ks);
 end
%     [~,~,xx_s,yy_s,par_s,pare_s,con0] = multi_scale_bls3(C_train_x_e,lab_tr_e,C_test_x_e,lab_te_e,s,C,h,f,f0,par_s,pare_s,1,N1,N2,fe,fe1,[]);
%     [~,~,x_s,y_s,~,~,~] = multi_scale_bls3(C_train_x_e,lab_tr_e,C_test_x_e,lab_te_e,s,C,h,f,f0,par_s,pare_s,2,N1,N2,fe,fe1,con0);

% 
% [~,lare]=max(lab_tr_r');
% [~,train_yy]=max(train_y');
% length(find(lare'==train_yy(dr_l)))/length(dr_l)

%%%%%%%%%%%%%∞¥±ª ‘≈≈–Ú
if mode==2
y_te=rank_y(y_s,de_s,y_l,de_l);
x_te=rank_y(x_s,de_s,x_l,de_l);
y_tr=[];
x_tr=[];
tnl=[];
tns=[];
end
if mode==1
y_tr=rank_y(yy_s,dr_s,yy_l,dr_l);
x_tr=rank_y(xx_s,dr_s,xx_l,dr_l);
y_te=[];
x_te=[];
train_y=rank_y(lab_tr_e,dr_s,lab_tr_r,dr_l);
PAR=gepar(par_l,pare_l,par_s,pare_s);
end
end

function PAR=gepar(par_l,pare_l,par_s,pare_s)
PAR.hs=par_l;
PAR.hse=pare_l;
PAR.hd=par_s;
PAR.hde=pare_s;
end

function [par_l,pare_l,par_s,pare_s]=repar(PAR)
if isempty(PAR)
    par_l=[];
    pare_l=[];
    par_s=[];
    pare_s=[];
else
    par_l=PAR.hs;
    pare_l=PAR.hse;
    par_s=PAR.hd;
    pare_s=PAR.hde;
end
end
function [kl,ks]=partk(leng,de_l,de_s)
for i=1:length(leng)
    k(i)=sum(leng(1:i));
end
k=[0,k];
for i=1:length(leng)
    lengl0(i)=length(intersect(find(de_l<=k(i+1)),find(de_l>k(i))));
    lengs0(i)=length(intersect(find(de_s<=k(i+1)),find(de_s>k(i))));
end
for i=1:length(leng)
    kl(i)=sum(lengl0(1:i));
    ks(i)=sum(lengs0(1:i));
end
kl=[0 kl];ks=[0 ks];
end