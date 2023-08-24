% function y_2=rank_y(y_s,de_s,y_l,de_l,y_m,de_m)
function y_2=rank_y(y_s,de_s,y_l,de_l)
    y_s1=[y_s de_s];
    y_l1=[y_l de_l];
%     y_m1=[y_m de_m];
%     y_2=sortrows([y_l1;y_m1;y_s1],2);
    y_2=sortrows([y_l1;y_s1],size(y_s1,2));
    y_2=y_2(:,1:end-1);
%     yy_21=sortrows([y_l0;y_s1(:,1:2)],2);
%     yy_21=yy_21(:,1);
end
