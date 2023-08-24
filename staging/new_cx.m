function [wx2,xx,f]=new_cx(x,n,mode,y)
a=[0.6 0.6 0.6 0.6 0.6 0.6 0.6 0.6 0.6 0.6 0.6 0.6 0.7 0.8 0.9 0.95 1 0.85 0.7 0.6 0.5 0.4 0.4 0.4 0.4 0.4 0.4 0.4 0.4 0.4 0.4 0.4 0.4 ];
n1=(length(a)-1)/2;
aa=a(n1-n:n1+n);
t=repmat(aa,1,5);
x(x==0)=eps;
x=x.*t;
m0=size(x,2);
u=mode*n+1;
m=m0/u;

for i=1:m
    xx=x(:,(i-1)*u+1:i*u);
    meanxr(:,i)=mean(xx')';
    stdxr(:,i)=std(xx')';
    skew(:,i)=skewness(xx')';
    kurto(:,i)=kurtosis(xx')';
    xx(xx==0)=eps;
    pe(:,i)=sum(xx'.*log(abs(xx')))';
    clear xx
end
for i=1:length(y)
meanxr1(i)=meanxr(i,y(i));
stdxr1(i)=stdxr(i,y(i));
skew1(i)=skew(i,y(i));
kurto1(i)=kurto(i,y(i));
pe1(i)=pe(i,y(i));
end
% meanxr=meanxr1';stdxr=stdxr1';skew=skew1';kurto=kurto1';pe=pe1';
% % meanxr=mean(x')';
% % stdxr=std(x')';
% % skew=skewness(x')';
% % kurto=kurtosis(x')';
% % x1=x;
% % x1(x1==0)=1;
% % pe=sum(x1'.*log(abs(x1')))';
f.mean=meanxr;
f.std=stdxr;
f.skew=mapminmax(skew')';
f.kurto=mapminmax(kurto')';
f.pe=mapminmax(pe')';
f.meanx=meanxr1';
f.stdx=stdxr1';
f.skewx=mapminmax(skew1)';
f.kurtox=mapminmax(kurto1)';
f.pex=mapminmax(pe1)';
% f=[meanxr stdxr skew kurto pe meanxr1' stdxr1' skew1' kurto1' pe1'];
% f=mapminmax(f')';
clear x1


xx=[];
for i=1:u
    x0=x(:,([1:m]-1)*u+i);
    xx=[xx x0];
end
% for i=1:u
%     for j=1:m
%     x0(j,:)=x(:,(j-1)*u+i);
%     end
%     xx=[xx x0];
% end
clear x0

% x0=xx(:,m*n+1:m*(n+1));
for i=1:u
    if i<n+1
        x1=xx(:,(i-1)*m+1:(i)*m);
        x0=xx(:,(i)*m+1:(i+1)*m);
    elseif i>n+1
        x0=xx(:,(i-1)*m+1:(i)*m);
        x1=xx(:,(i-2)*m+1:(i-1)*m);
    else
        x1=xx(:,(i-1)*m+1:(i)*m);
        
    end
    % %     if find(x1==0,1)
    % %        x1(x1==0)=1;
    % %     end
    %     wx(:,(i-1)*m+1:i*m)=x0.*log( abs(x0./x1) ) .* (x1 ./abs(x1));%ml
    if i~=n+1
%         wx1(:,(i-1)*m+1:i*m)=(x1+x0)/2;
        wx2(:,(i-1)*m+1:i*m)=x1+abs(x1-x0);
%         wx3(:,(i-1)*m+1:i*m)=x1-x0;
    else
%         wx1(:,(i-1)*m+1:i*m)=x1;
        wx2(:,(i-1)*m+1:i*m)=x1;
%         wx3(:,(i-1)*m+1:i*m)=x1;
    end
end
end
