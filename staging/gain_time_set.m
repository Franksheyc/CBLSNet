function [C_test_x_s,C_test_x_l,C_train_x_l,C_train_x_s,lab_te_s,lab_te_l,lab_tr_s,lab_tr_l,de_s,de_l,dr_s,dr_l,y_l,yy_l,mk]=gain_time_set(C_test_x,C_train_x,y,yy,test_y,train_y,cr,ce,m)
%%%%%%%%%  l is the right part of data
%%%%%%%%%  s is the wrong part of data
% %     de_l=find(test_yy==y);de_s=find(test_yy~=y);
% %     dr_l=find(train_yy==yy);dr_s=find(train_yy~=yy);
by=findcc(y,20);
byy=findcc(yy,20);
cr=cr.*byy';
ce=ce.*by';

%     m1=m/2;

k=1;mark0=0;
while k==1
    k=2;
        m1=(quantile(cr,m)+quantile(ce,m))/2;
        de_l=find(ce>=m1);de_s=find(ce<m1);
        dr_l=find(cr>=m1);dr_s=find(cr<m1);
        sl=length(dr_l)/length(dr_s);
        se=length(de_l)/length(de_s);
    % %     mark=abs(se-sl);
    %     con1=intersect(find(sl-se>-0.0015),find(sl-se<0.2));
    %     con2=intersect(intersect(find(se<0.5),find(se>0.2)),find(sl>0.25));
    %     mark=intersect(con1,con2);
%     [de_l,de_s,dr_l,dr_s]=partsample(m,yy,y,cr,ce);
    
%     sl=length(dr_l)/length(yy);
%     se=length(de_l)/length(y);
    %         mark=sl-se;
    a=histtt(yy(dr_l),5);
    b=histtt(y(de_l),5);
%     a=hist(yy(dr_l),5);
%     b=hist(y(de_l),5);
    a=a/sum(a);
    b=b/sum(b);
    mark=sum(abs(a-b))/5;
    mk=[sl se mark a b];
    mark0=mark;
        if mark>0.1
            mi=1;
            if mark0>mark
                mi=-1;
            end
            m=m+0.05*mi;
            if m>0.8 || m<0.5
                break
            end
            k=1;
        end
        mark0=mark;
end
%     [dr_l,dr_s]=gedre(yy);
%     [de_l,de_s]=gedre(y);

C_test_x_l=C_test_x(de_l,:);
C_test_x_s=C_test_x(de_s,:);
C_train_x_l=C_train_x(dr_l,:);
C_train_x_s=C_train_x(dr_s,:);

yy_l=yy(dr_l,:);
%     yy_l0=[yy_l0 dr_l];
y_l=y(de_l,:);
%     y_l0=[y_l0 de_l];


lab_tr_l=train_y(dr_l,:);
lab_tr_s=train_y(dr_s,:);
lab_te_s=test_y(de_s,:);
lab_te_l=test_y(de_l,:);
end
%%%%%%%%%%%%%%
% function b=findcc(label)
%     b=1;
% for i=2:length(label)-1
%     if label(i)==label(i-1) && label(i)==label(i+1)
%         b(i)=1.5;
%     elseif label(i)==label(i-1) && label(i)~=label(i+1)
%         b(i)=1.1;
%     elseif label(i)~=label(i-1) && label(i)==label(i+1)
%         b(i)=1;
%     else
%         b(i)=0.9;
%     end
%
% end
%     b(end+1)=1;
% end

function b=findcc(label,n)
% if length(unique(label(1:n)))==1
%     b=ones(1,n)*1.8;
% else
    b=ones(1,n);
% end
parfor i=n+1:length(label)-n
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
% if length(unique(label(end-n+1:end)))==1
%     b(end+1:end+n)=ones(1,n)*1.8;
% else
    b(end+1:end+n)=ones(1,n);
% end
end

function [de_l,de_s,dr_l,dr_s]=partsample1(m,yy,y,cr,ce)
de_l=[];de_s=[];dr_l=[];dr_s=[];
m1=(quantile(cr,m)+quantile(ce,m))/2;
de_l=find(ce>=m1);de_s=find(ce<m1);
dr_l=find(cr>=m1);dr_s=find(cr<m1);
[a]=hist(de_l);[b]=hist(dr_l);

end

function [de_l,de_s,dr_l,dr_s]=partsample(m,yy,y,cr,ce)
de_l=[];de_s=[];dr_l=[];dr_s=[];
for i=1:5
    myy=find(yy==i);
    my=find(y==i);
    crr=cr(myy);
    cee=ce(my);
    m1(i)=(quantile(crr(m))+quantile(cee(m)))/2;
    e_l=my(find(cee>=m1(i)));e_s=my(find(cee<m1(i)));
    r_l=myy(find(crr>=m1(i)));r_s=myy(find(crr<m1(i)));
    de_l=[de_l;e_l];
    de_s=[de_s;e_s];
    dr_l=[dr_l;r_l];
    dr_s=[dr_s;r_s];
end
de_l=sort(de_l);
de_s=sort(de_s);
dr_l=sort(dr_l);
dr_s=sort(dr_s);
end
function [myy]=histtt(yy,n)
for i=1:n
    myy(i)=length(find(yy==i))/length(yy);
end
end