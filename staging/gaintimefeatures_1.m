function [C_xr,train_xm,cr,Cr,C_xe,test_xm,ce,Ce,sumweight_tr,sumweight_te,C_yyr,C_yye,yrf,yef]=gaintimefeatures_1(xx,x,y,yy,time_related,ss,ss1,leng,mode,train_yy,si,bs)
% function [C_xr,train_xm,cr,Cr,C_xe,test_xm,ce,Ce,sumweight_tr,sumweight_te,C_yyr,C_yye,wxr,Cxxr,wxe,Cxxe,x_b,xx_b,crf,cef,yrf,yef]=gaintimefeatures_1(xx,x,y,yy,time_related,ss,ss1,leng,mode,si,bs)
%   x=softmax_1(x);
%   xx=softmax_1(xx);
%    x_b=x.*(2*randn(size(x,1),size(x,2))-1)+x;
%    xx_b=xx.*(2*rand(size(xx,1),size(xx,2))-1)+xx;
% %  numrelcal=1;
% %  accr=length(find(yy==train_yy));
% %  accr1=0;
% % while numrelcal<2 || accr<accr1
   [cr,ce]=crce(xx,x,si,bs);
   
  [ae,be]=sort(x,2);
  [ar,br]=sort(xx,2);

%   hr=find(abs(ar(:,5))>2);
%   he=find(abs(ae(:,5))>2);
%   ar(hr,:)=ar(hr,:)./ar(hr,5);
%   ae(he,:)=ae(he,:)./ae(he,5);
  
%   ce(:,1)=ae(:,5)./ae(:,4);%%%%%%%置信度之比
%   cr(:,1)=ar(:,5)./ar(:,4);
%   ce(:,2)=ae(:,5)./sum(ae,2);%%%%%%%%%%最大置信度占总体的比例
%   cr(:,2)=ar(:,5)./sum(ar,2);
%   cr(:,2)=cr(:,2)./max(cr);
%     [cr(:,2),ps]=mapminmax(cr(:,2),0,1);
%     ce(:,2)=mapminmax(ce(:,2),ps);

  
%   ce(:,2)=ce(:,2)./max(ce);
    %%%%%%置信度之比归一化

%%%%%%%%%%%%%%%%%%%%%%%%%%按被试做
    y_tr=label2lab(yy);%
    y_te=label2lab(y);%
%     x_pr=1-xx;
%     x_pe=1-x;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%排序
%     train_xm1=sleep_embeding(br(:,5));%%%%把标签按embedding方式输出
%     test_xm1=sleep_embeding(be(:,5));
%     train_xm2=sleep_embeding(br(:,4));
%     test_xm2=sleep_embeding(be(:,4));
    train_xm1=label2lab(br(:,5));
    test_xm1=label2lab(be(:,5));
    train_xm2=label2lab(br(:,4));
    test_xm2=label2lab(be(:,4));
    train_xm=zscore([train_xm1 train_xm2]);%1st 2nd position
    test_xm=zscore([test_xm1 test_xm2]);

sumweight_tr=[ar(:,5).*train_xm1 ar(:,4).*train_xm2];
sumweight_te=[ae(:,5).*test_xm1 ae(:,4).*test_xm2];
    
    gg1=1;gg2=0;
    singlelen=time_related*mode+1;
    for i=1:length(ss)
        gg=leng((ss(i)));
        gg2=gg2+gg;
        for ik=1:size(y_tr,2)
%         C_yr(gg1:grg2,(ik-1)*singlelen+1:singlelen*ik)=time_causal(y_tr(gg1:gg2,ik),time_related);%
        C_xr(gg1:gg2,(ik-1)*singlelen+1:singlelen*ik)=time_causal(xx(gg1:gg2,ik),time_related,mode);
%         C_xpr(gg1:gg2,(ik-1)*singlelen+1:singlelen*ik)=time_causal(x_pr(gg1:gg2,ik),time_related);
        end
% C_xr(gg1:gg2,:)=time_causal(xx(gg1:gg2,:),time_related,mode);
        C_yyr(gg1:gg2,1:singlelen)=time_causal(yy(gg1:gg2),time_related,mode);
%         Cr1(gg1:gg2,1:singlelen)=time_causal(cr(gg1:gg2,1),time_related,mode);
        Cr(gg1:gg2,1:singlelen)=time_causal(cr(gg1:gg2),time_related,mode);
%         for ik=1:size(train_xm,2)
%             Train_xm(gg1:gg2,(ik-1)*singlelen+1:singlelen*ik)=time_causal(train_xm(gg1:gg2,ik),time_related);
%         end
        gg1=gg1+gg;
    end
    gg1=1;gg2=0;
    for i=1:length(ss1)
        gg=leng((ss1(i)));
        gg2=gg2+gg;
        for ik=1:size(y_te,2)
%         C_ye(gg1:gg2,(ik-1)*singlelen+1:singlelen*ik)=time_causal(y_te(gg1:gg2,ik),time_related);
        C_xe(gg1:gg2,(ik-1)*singlelen+1:singlelen*ik)=time_causal(x(gg1:gg2,ik),time_related,mode);
%         C_xpe(gg1:gg2,(ik-1)*singlelen+1:singlelen*ik)=time_causal(x_pe(gg1:gg2,ik),time_related);
        end
% C_xe(gg1:gg2,:)=time_causal(x(gg1:gg2,:),time_related,mode);
        C_yye(gg1:gg2,1:singlelen)=time_causal(y(gg1:gg2),time_related,mode);
%         Ce1(gg1:gg2,1:singlelen)=time_causal(ce(gg1:gg2,1),time_related,mode);
        Ce(gg1:gg2,1:singlelen)=time_causal(ce(gg1:gg2),time_related,mode);
%         for ik=1:size(test_xm,2)
%             Test_xm(gg1:gg2,(ik-1)*singlelen+1:singlelen*ik)=time_causal(test_xm(gg1:gg2,ik),time_related);
%         end
        gg1=gg1+gg;
    end
    
%     nh=markss(mode,time_related,x,y,singlelen);
%     C_xe(nh)=eps;
%     nh1=markss(mode,time_related,xx,yy,singlelen);
%     C_xr(nh1)=eps;
    
%     [wxr,Cxxr]=weight_ml(C_xr,time_related);
%     [wxe,Cxxe]=weight_ml(C_xe,time_related);
    
%     meancr=mean(Cr')';
%     stdcr=std(Cr')';
%     meance=mean(Ce')';
%     stdce=std(Ce')';
%     crf=[meancr stdcr];
%     cef=[meance stdce];
    [continues_r,transform_r,stage_ratio_r,current_ratio_r,domain_kindn_r]=weight_features(C_yyr,time_related);%当前帧的连续性，条变形，各种帧的占比，当前帧的占比，共几种帧
    [continues_e,transform_e,stage_ratio_e,current_ratio_e,domain_kindn_e]=weight_features(C_yye,time_related);
    yrf=[continues_r,transform_r,stage_ratio_r,current_ratio_r,domain_kindn_r];
    yef=[continues_e,transform_e,stage_ratio_e,current_ratio_e,domain_kindn_e];
% %     n=5;
% %     pyr=(yrf(:,end-1)+yrf(:,end))/2;
% %     pye=(yef(:,end-1)+yef(:,end))/2;
% %     xx1=conti_prob(xx,n,pyr,yy);
% %     x1=conti_prob(x,n,pye,y);
% %     [~,yy1]=max(xx1');
% %     yy1=yy1';
% %     [~,y1]=max(x1');
% %     y1=y1';
% %     accr1=length(find(yy1==train_yy));
% %     if accr1>accr
% %         y=y1;
% %         yy=yy1;
% %         x=x1;
% %         xx=xx1;
% %     end    
% %      numrelcal = numrelcal+1;
% % end
end

function nh=markss(mode,time_related,x,y,singlelen)
if mode == 1
    m=1:time_related;
else
    m=[1:time_related time_related+2:2*time_related+1];
end


for j=1:size(x,2)
    h=[];
    k=setdiff(1:size(x,2),j);
    for i=1:size(x,2)-1
        h=[h (k(i)-1)*singlelen+m];
    end
    H(j,:)=h;
    h=[];
end

nh=zeros(size(y,1),size(H,2));
for i=1:size(y,1)
   nh(i,:)=H(y(i),:);
end
end