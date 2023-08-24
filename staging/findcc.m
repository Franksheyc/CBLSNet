function b=findcc(label,n)
b=ones(1,n);
for i=n+1:length(label)-n
    a=unique(label(i-n:i+n));
    switch length(a)
        case 1
            b(i)=1.8;
        case 2
            b(i)=1.2;
        case 3
            b(i)=1;
        case 4
            b(i)=0.9;
        case 5
            b(i)=0.7;
    end
end
b(end+1:end+n)=ones(1,n);
end