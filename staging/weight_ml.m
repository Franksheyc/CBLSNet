function [wx,xx]=weight_ml(x,n)
m0=size(x,2);
u=2*n+1;
m=m0/u;
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

x0=xx(:,m*n+1:m*(n+1));
for i=1:2*n
    h=1*(i>n);
    x1=xx(:,(i-1+h)*m+1:(i+h)*m);
    if find(x1==0,1)
       x1(x1==0)=1; 
    end
%     wx(:,(i-1)*m+1:i*m)=x0.*log( abs(x0./x1) ) .* (x1 ./abs(x1));%ml
    wx(:,(i-1)*m+1:i*m)=x0+zscore(  sqrt( abs(x0./x1) )'  )';
end