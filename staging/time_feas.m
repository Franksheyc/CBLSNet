function [C_train_x,C_test_x,cr1,ce1,ftr,fte,f,f0]=time_feas(xx,x,y,yy,time_related,ss,ss1,leng,mod,train_yy,train_y,si,bs)
% test_yy=[];m=0;
% train_y=label2lab(train_yy);
% [~,~,~,~,~,f1,~,~] = Evaluation(yy,train_y,5);
% ac1=sum(f1);ac2=ac1+1/length(yy);ac21=ac1+1/length(yy);
% h=0;
% tw=2*time_related+1;
% while ac2>ac1
%     h=h+1;
% [C_xr1,train_xm1,cr1,C_xe1,test_xm1,ce1,sumweight_tr1,sumweight_te1,C_yr1,C_ye1,Cr2,Ce2,wxr1,Cxxr1,wxe1,Cxxe1,x_b,xx_b,crf1,cef1,yrf1,yef1]=gaintimefeature(xx,x,y,yy,time_related,ss,ss1,leng,mode,si,bs);
[C_xr1,train_xm1,cr1,C_xe1,test_xm1,ce1,sumweight_tr1,sumweight_te1,C_yr1,C_ye1,Cr2,Ce2,yrf1,yef1]=gaintimefeature(xx,x,y,yy,time_related,ss,ss1,leng,mod,train_yy,si,bs);
% [Cwxr1,Cwxr2,Cwxr3,crx,ftr]=new_cx(C_xr1,time_related,mode,yy);
% [Cwxe1,Cwxe2,Cwxe3,cex,fte]=new_cx(C_xe1,time_related,mode,y);
pyr=(yrf1(:,end-1)+yrf1(:,end))/2;
pye=(yef1(:,end-1)+yef1(:,end))/2;

[Cwxr2,crx,ftr]=new_cx(C_xr1,time_related,mod,yy);
[Cwxe2,cex,fte]=new_cx(C_xe1,time_related,mod,y);
[xx,yy]=addtimef(xx,yy,ftr,train_yy);

% n=5;
% xx=conti_prob(xx,n,pyr,yy);
% x=conti_prob(x,n,pye,y);
% [x,y]=addtimef(x,y,fte,test_yy);
% % [~,~,~,~,~,f1,~,~] = Evaluation(yy,train_y,5);
% % ac2=sum(f1);
% % if ac2<ac21
% %     ac2=0;
% % else
% %     m=m+1;
% % %     [C_xr1,train_xm1,cr1,C_xe1,test_xm1,ce1,sumweight_tr1,sumweight_te1,C_yr1,C_ye1,Cr2,Ce2,wxr1,Cxxr1,wxe1,Cxxe1,x_b,xx_b,crf1,cef1,yrf1,yef1]=gaintimefeature(xx,x,y,yy,time_related,ss,ss1,leng,mode);
% % %     [Cwxr2,crx,ftr,per]=new_cx(C_xr1,time_related,mode,yy);
% % %     [Cwxe2,cex,fte,pee]=new_cx(C_xe1,time_related,mode,y);
% % end
% % ac21=ac2;
% % 
% % if h>2
% %    continue 
% % end
% % 
% % end
%         C_train_x = [C_xr1 mean(C_xr1')' std(C_xr1')' skewness(C_xr1')' kurtosis(C_xr1')' Cwxr2 train_xm1 Cr2 sumweight_tr1 C_yr1 crx crf1 yrf1];%用训练结果  Cr2 crf1 yrf1    sumweight_tr Time_r  Cxxr1 wxr1 xx_b
%         C_test_x =[C_xe1 mean(C_xe1')' std(C_xe1')' skewness(C_xe1')' kurtosis(C_xe1')' Cwxe2 test_xm1 Ce2 sumweight_te1 C_ye1 cex cef1 yef1];%              Ce2 cef1 yef1    sumweight_te Time_e  Cxxe1 wxe1 x_b
% if m
    C_train_x = [C_xr1 ftr.pe ftr.kurto ftr.pex zscore(Cwxr2) Cr2 sumweight_tr1  C_yr1/5 pyr yrf1(:,1:7)];%用训练结果  Cr2 crf1 yrf1    sumweight_tr Time_r  Cxxr1 wxr1 xx_b    ftr(:,[6 8 10])
    C_test_x =[C_xe1 fte.pe fte.kurto fte.pex zscore(Cwxe2) Ce2 sumweight_te1 C_ye1/5 pye yef1(:,1:7)];%
%     C_train_x = [C_xr1 ftr.pe ftr.kurto ftr.pex zscore(Cwxr2)  sumweight_tr1 C_yr1/5  pyr  yrf1(:,3:7)];%用训练结果  Cr2 crf1 yrf1    sumweight_tr Time_r  Cxxr1 wxr1 xx_b    ftr(:,[6 8 10])
%     C_test_x =[C_xe1 fte.pe fte.kurto fte.pex zscore(Cwxe2)  sumweight_te1 C_ye1/5 pye yef1(:,3:7)];%
    [f,f0]=blstimekindt(time_related,1);
% else
%     C_train_x = [C_xr1  ftr.pe ftr.kurto ftr.pex  zscore(Cwxr2) Cr2 sumweight_tr1 C_yr1/5 pyr yrf1(:,1:7)];%用训练结果  Cr2 crf1 yrf1    sumweight_tr Time_r  Cxxr1 wxr1 xx_b
%     C_test_x = [C_xe1  fte.pe fte.kurto fte.pex  zscore(Cwxe2) Ce2 sumweight_te1 C_ye1/5 pye yef1(:,1:7)];%
% %     C_train_x = [C_xr1  ftr.pe ftr.kurto ftr.pex  zscore(Cwxr2) sumweight_tr1 C_yr1/5  pyr yrf1(:,3:7)];%用训练结果  Cr2 crf1 yrf1    sumweight_tr Time_r  Cxxr1 wxr1 xx_b
% %     C_test_x = [C_xe1  fte.pe fte.kurto fte.pex  zscore(Cwxe2)  sumweight_te1 C_ye1/5 pye yef1(:,3:7)];%
%     [f,f0]=blstimekindt(time_related,0);
% end

% mp=find(train_yy==mode(train_yy));
% mp1=setdiff(1:length(train_yy),mp);
% yrr=yrf1(:,3:7);
% for i=1:length(mp1)
% ms=train_yy(mp1(i));
% sj=-1*ones(1,size(train_y,2));
% sj(ms)=1;
% % if train_yy(time_related+1)
% train_y(mp1(i),:)=train_y(mp1(i),:)+0.05*yrr(mp1(i),:).*sj;
% % end
% end
end

function [f,f0]=blstimekindt(time_related,m)
h=time_related*2+1;
f=[ [h*ones(1,5)] [5 5 1]  [5*ones(1,h)] h 10 h 1 2 5];
f0=[1:5];
end

function [C_xr1,train_xm1,cr1,C_xe1,test_xm1,ce1,sumweight_tr1,sumweight_te1,C_yr1,C_ye1,Cr2,Ce2,yrf1,yef1]=gaintimefeature(xx,x,y,yy,time_related,ss,ss1,leng,mode,train_yy,si,bs)
% function [C_xr1,train_xm1,cr1,C_xe1,test_xm1,ce1,sumweight_tr1,sumweight_te1,C_yr1,C_ye1,Cr2,Ce2,wxr1,Cxxr1,wxe1,Cxxe1,x_b,xx_b,crf1,cef1,yrf1,yef1]=gaintimefeature(xx,x,y,yy,time_related,ss,ss1,leng,mode,si,bs)
        C_xr1=[];
        train_xm1=[];
        cr1=[];
        C_xe1=[];
        test_xm1=[];
        ce1=[];
        sumweight_tr1=[];
        sumweight_te1=[];
        C_yr1=[];
        C_ye1=[];
%         wxr1=[];
%         Cxxr1=[];
%         wxe1=[];
%         Cxxe1=[];
%         crf1=[];
%         cef1=[];
        yrf1=[];
        yef1=[];
        tic
        for i=1:length(time_related)
            %     [C_xr,train_xm,cr,C_yr,Cr1,Cr2,C_xe,test_xm,ce,C_ye,Ce1,Ce2,sumweight_tr,sumweight_te]=gaintimefeatures_1(xx,x,y,yy,time_related,ss,ss1,leng);
%             [C_xr,train_xm,cr,Cr1,Cr2,C_xe,test_xm,ce,Ce1,Ce2,sumweight_tr,sumweight_te,C_yr,C_ye]=gaintimefeatures_1(xx,x,y,yy,time_related(i),ss,ss1,leng,mode);
%             [C_xr,train_xm,cr,Cr2,C_xe,test_xm,ce,Ce2,sumweight_tr,sumweight_te,C_yr,C_ye,wxr,Cxxr,wxe,Cxxe,x_b,xx_b,crf,cef,yrf,yef]=gaintimefeatures_1(xx,x,y,yy,time_related(i),ss,ss1,leng,mode,si,bs);
[C_xr,train_xm,cr,Cr2,C_xe,test_xm,ce,Ce2,sumweight_tr,sumweight_te,C_yr,C_ye,yrf,yef]=gaintimefeatures_1(xx,x,y,yy,time_related(i),ss,ss1,leng,mode,train_yy,si,bs);
            C_xr1=[C_xr1 C_xr];
            train_xm1=[train_xm1 train_xm];
            cr1=[cr1 cr];
            C_xe1=[C_xe1 C_xe];
            test_xm1=[test_xm1 test_xm];
            ce1=[ce1 ce];
            sumweight_tr1=[sumweight_tr1 sumweight_tr];
            sumweight_te1=[sumweight_te1 sumweight_te];
            C_yr1=[C_yr1 C_yr];
            C_ye1=[C_ye1 C_ye];
%             wxr1=[wxr1 wxr];
%             Cxxr1=[Cxxr1 Cxxr];
%             wxe1=[wxe1 wxe];
%             Cxxe1=[Cxxe1 Cxxe];
%             crf1=[crf1 crf];
%             cef1=[cef1 cef];
            yrf1=[yrf1 yrf];
            yef1=[yef1 yef];
        end
%         time_remake=toc;
end