function x=getpartdata(X,s,k)
%%%%%%% k
%%%%%%% X data
%%%%%%% s the requierd data number
%%%%%%% k the position of X
%%%%%%% k must longer than s 1
s=s+1;
if size(X,1)==1
    X=X';
end
    k=[0 k];

 
l1=1; l2=0; L=[]; L1=[];  
for i=1:length(s)   
    len=k(s(i))-k(s(i)-1);   
    l2=l2+len;    
    x(l1:l2,:)=X(k(s(i)-1)+1:k(s(i)),:);    
    if size(x,2)>63
       L=[L (1:len)/len*0.8];
       [~,lk]=sort(abs(X(k(s(i)-1)+1:k(s(i)),67)));%63 2up 3down  64  4down   65wending
       L1=[L1 lk'/len*0.8];
    end
    l1=l1+len;
end
% x=[x L' L1'];
x=[x L'];

