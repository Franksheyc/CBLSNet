function [activity,mobility,complexity]=Hjorth_features(dat)
dat(dat==0)=0.001;
[len,wid]=size(dat);
m=mean(dat')';
me=ones(len,wid).*m;
sd=std(dat')';
sd1=std(diff([zeros(1,wid);dat])')';
sd2=std(diff([zeros(2,wid);dat],2)')';
activity=1/len*(sum((dat-me)'))';

mobility=sd1./sd;
complexity=(sd1./sd)./(sd2./sd1);

% aa=diff([ zeros(size(dat,1),1) dat]')';
% aa(aa==0)=eps;
% mobility=sqrt(var(aa'./dat'))';
% aaa=diff([zeros(size(aa,1),1) aa]')';
% complexity=sqrt( sqrt(var(aaa'./aa'))./sqrt(var(aa'./dat')) )';
end
