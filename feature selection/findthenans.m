function a = findthenans(C3,C4,EMG,EOG,c31,c41,h1,h2,h3)
a=[];
if nargin==4
    a1=findthenan(C3);
    a2=findthenan(C4);
    a3=findthenan(EOG);
    a4=findthenan(EMG);
    a=[a1;a2;a3;a4];
    
elseif nargin==6
    a1=findthenan(C3);
    a2=findthenan(C4);
    a3=findthenan(EOG);
    a4=findthenan(EMG);
    a5=findthenan(c31);
    a6=findthenan(c41);
    a=[a1;a2;a3;a4;a5;a6];
else
    a1=findthenan(C3);
    a2=findthenan(C4);
    a3=findthenan(EOG);
    a4=findthenan(EMG);
    a5=findthenan(c31);
    a6=findthenan(c41);
    a7=findthenan(h1);
    a8=findthenan(h2);
    a9=findthenan(h3);
    
    a=[a1;a2;a3;a4;a5;a6;a7;a8;a9];
end
a=unique(a);
end

function a=findthenan(C3)
for i=1:size(C3,2)
    m=find (isnan(C3(:,i))==1);
    n=find (isinf(C3(:,i))==1);
    s1{i}=unique(union(m,n));
end
e=[];
for j=1:length(s1)
    e=[e;s1{j}];
end
a=unique(e);
end
