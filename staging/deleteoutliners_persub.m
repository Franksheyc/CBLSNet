function [B,st]=deleteoutliners_persub(A,epsino,st)
st=st(1:size(A,2));
B=[];
for i=1:size(A,2)
    a=A(:,i);
    if length(find(isinf(a)==1))
       a(isinf(a))=eps;
    end
    if length(find(isnan(a)==1))
       a(isnan(a))=max(a); 
    end
    if (max(a)>mean(a)+3*std(a) || min(a)<mean(a)-3*std(a)) && max(a)-min(a)>2*(mean(a)+3*std(a))
        s=find(a<0, 1);
        if ~isempty(s)
            su=prctile(a,99);
            sd=prctile(a,1);
            a1=a(intersect(find(a<su),find(a>sd)));
        else
            su=prctile(a,99);
            a1=a(a<su);
        end
%         stander1=(mean(a1')+3*std(a1'))';
%         stander2=(mean(a1')-3*std(a1'))';
        stander1=mean(a1)+3*std(a1);
        stander2=mean(a1)-3*std(a1);
        if sum(st)~=0
        stander1=min(stander1,st(i));
        stander2=max(stander2,-st(i));
        end
        outliner1=a>stander1*(1+epsino);
        outliner2=a<stander2*(1+epsino);
        a(outliner1)=stander1;
        a(outliner2)=stander2;
    else
        stander1=mean(a)+3*std(a);
        stander2=mean(a)-3*std(a);
    end
    B=[B a];
    if length(find(a>0))==length(a)
        st1(i)=stander1;
    else
        st1(i)=(abs(stander1)+abs(stander2))/2;
    end
end
if sum(st)~=0
    for i=1:length(st)
        if abs(st(i)-st1(i))/abs(st(i))>=10
            st0(i)=0;
        else
            st0(i)=st1(i);
        end
    end
else
    st0=st1;
end
st=(st+st0)/2;