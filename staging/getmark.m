function L=getmark(s,k)
s=s+1;
k=[0 k];
 
l1=1; l2=0; L=[]; 
for i=1:length(s)   
    len=k(s(i))-k(s(i)-1);   
    l2=l2+len;    
    L=[L k(s(i)-1)+1:k(s(i))];

%     k(s(i)-1)+1:k(s(i)
%     if size(x,2)>63
%        L=[L (1:len)/len*0.8];
%        [~,lk]=sort(abs(X(k(s(i)-1)+1:k(s(i)),67)));%63 2up 3down  64  4down   65wending
%        L1=[L1 lk'/len*0.8];
%     end
    l1=l1+len;
end


end