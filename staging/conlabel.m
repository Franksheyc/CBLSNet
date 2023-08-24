function [con,dconj,dcon]=conlabel(label,dikind)
s=0;
kind=setdiff(1:5,dikind);
mark=[];m=0;
for i=1:length(kind)
label1=find(label==kind(i));
mark=[mark;label1];
mark=sort(mark);
end

[con,dconj,dcon]=checkcon(label);
con=intersect(con,mark);
dconj=intersect(dconj,mark);
dcon=intersect(dcon,mark);
% for i=1:length(label)-1
%     if label(i+1)==label(i)
%         m=m+1;
%         s(m)=i;
%     end
% end
% s1=intersect(s,mark);
% s=setdiff(1:length(label),s1);

