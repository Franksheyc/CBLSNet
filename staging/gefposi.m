function [F,f00]=gefposi(req,f,nfp,f0)

for i=1:length(f)
    f1(i)=sum(f(1:i));
end
F1=[0 f1(1:end-1)]';
F0=[F1+1 f1'];
clear F F1
for i=1:size(F0,1)
   F{i}=F0(i,1):F0(i,2); 
end

F1=F(1:length(f)-nfp);
F2=F(length(f)-nfp+1:length(f));
if ~isempty(req)
for j=1:length(req)
% F1{j+length(f)}=findposi(req{j},f1,fe);
F3{j}=req{j};
end
else
    F3=[];
end
clear F
F=[F1 F3 F2];
f00=[f0 f0(end)+1:f0(end)+length(F3)];
end