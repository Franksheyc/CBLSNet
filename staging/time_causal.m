function yy=time_causal(lab,time_related,mode)
[a,b]=size(lab);
if b>a
    lab=lab';
end
[a,b]=size(lab);
% yy=zeros(length(lab),(time_related+1)*c);
% yy(:,time_related*c+1:(time_related+1)*c)=lab;
%    a=time_related+1;
%    b=length(lab);
%    for i=1:time_related
%        for j=1:time_related
%        yy(i,c*(j-1)+1:c*j)=0;
%        end
%        yy(a+1-i:end,c*(i-1)+1:c*i)=lab(1:end-a+i,:);
%    end
%    for i=1:time_related
%        for j=1:time_related
%        yy(b+1-i,c*(a+j-1)+1:c*(a+j))=0;
%        end
%        yy(1:b-a+i,c*(a+i-1)+1:c*(a+i))=lab(a-i+1:end,:);
%    end
%    for i=1:time_related/2
%        yy1=yy(:,c*(i-1)+1+a*c:c*i+a*c);
%        yy(:,c*(i-1)+1+a*c:c*i+a*c)=yy(:,a*c+c*(a-i-1)+1:c*(a-i)+a*c);
%        yy(:,a*c+c*(a-i-1)+1:c*(a-i)+a*c)=yy1;
%    end
   m=2*time_related+1;
   yy=zeros(a+m-1,m*b);
   y=zeros(a+m-1,1);
   for j=1:b
      for i=1:m
          y(m-i+1:a+m-i,i)=lab(:,j);
      end
      yy(:,(j-1)*m+1:j*m)=y;
   end
   yy=yy(time_related+1:time_related+length(lab),:);
   if mode==1
       yy0=[];
       for i=1:b
           yy1=yy(:,(i-1)*m+1:(i-1)*m+time_related+1);
           yy0=[yy0 yy1];
       end
       yy=yy0;
   end
end