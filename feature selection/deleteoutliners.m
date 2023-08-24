function B=deleteoutliners(A,epsino)
B=[];
for i=1:size(A,2)
    a=A(:,i);
    if (max(a)>mean(a)+3*std(a) || min(a)<mean(a)-3*std(a)) && std(a)>5
        s=find(a<0, 1);
        if ~isempty(s)
            su=prctile(a,99);
            sd=prctile(a,1);
            a1=a(intersect(find(a<su),find(a>sd)));
        else
            su=prctile(a,99);
            a1=a(a<su);
        end
        
        stander1=(mean(a1')+4*std(a1'))';
        stander2=(mean(a1')-4*std(a1'))';
        outliner1=a>stander1*(1+epsino);
        outliner2=a<stander2*(1+epsino);
        a(outliner1)=stander1;
        a(outliner2)=stander2;
    end
    B=[B a];
end