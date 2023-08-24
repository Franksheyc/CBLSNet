function x=conti_prob(x,n,pyr,yy)
[~,pyrr]=findpeaks(-pyr);
n=3;
% x=[ones(n,5);x;ones(n,5)];
% yy=[ones(n,1);yy;ones(n,1)];
len=length(yy);
for i=1:length(pyrr)
    mh=adjoin_label(yy,pyrr,i,len,n);
%     p1=mode(yy(pyrr(i)-n):yy(pyrr(i)+n));
%     mh=p1==yy(pyrr(i));
    if mh
        x(pyrr(i),:)=jumping(x(pyrr(i),:),yy(pyrr(i)));     
    end
end

pyrr1=setdiff(1:length(pyr),pyrr);
for i=1:length(pyrr1)
%     p2=mode(yy(pyrr1(i)-n):yy(pyrr1(i)+n));
%     mh2=p2==yy(pyrr1(i));
    mh2=adjoin_label(yy,pyrr1,i,len,n);
%     mh2=find(pyr(pyrr1)~=yy(pyrr1));
    if mh2
        x(pyrr1(i),:)=continuous(x(pyrr1(i),:),yy(pyrr1(i)));
    end
end
% x=x(n+1:end-n,:);
end

function mh=adjoin_label(yy,pyrr,i,len,n)
    p1=mode( yy( max(pyrr(i)-n,1) : min(pyrr(i)+n,len) )  );
    mh=p1==yy(pyrr(i));
end

function x=jumping(x,k)
switch k
    case 2
        x(3)=x(3)*2;
    case 3
        x(2)=x(2)*2;
    case 4
        x(5)=x(5)*2;
end
end

function x=continuous(x,k)

    x(k)=x(k)*0.5;

end