function [x_2,y_2]=ranky(y_s,de_s,x_s,y_l,x_l,de_l)
    y_s1=[y_s de_s x_s];
    y_l1=[y_l de_l x_l];  
    y_2=sortrows([y_l1;y_s1],2);
    x_2=y_2(:,3:7);
    y_2=y_2(:,1);
%     yy_21=sortrows([y_l0;y_s1(:,1:2)],2);
%     yy_21=yy_21(:,1);
end
