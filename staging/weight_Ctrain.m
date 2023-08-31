function [C_xr2,C_xe2]=weight_Ctrain(C_xr1,C_xe1,time_related,mod)

mi=time_related+1;
xc1=ones(1,2*time_related+1);
if mod==1
%     xc1=[0.2 0.3 0.4 0.5 0.8 1 0.8 0.5 0.4 0.3 0.2];
xc1=[0.6 0.6 0.6 0.6 0.6 0.6 0.6 0.6 0.6 0.6 0.6 0.6 0.7 0.8 0.9 0.95 1 0.85 0.7 0.6 0.5 0.4 0.4 0.4 0.4 0.4 0.4 0.4 0.4 0.4 0.4 0.4 0.4];
% xc1=[0.4 0.5 0.6 0.7 0.8 1 0.8 0.7 0.6 0.5 0.4];
else
%     xc1=[0.2 0.3 0.3 0.5 0.8 1 0.3 0.2 0.1 0.1 0.05];
xc1=[0.6 0.6 0.6 0.6 0.6 0.6 0.6 0.6 0.6 0.6 0.6 0.6 0.7 0.8 0.9 0.9 1 0.8 0.7 0.6 0.5 0.4 0.4 0.4 0.4 0.4 0.4 0.4 0.4 0.4 0.4 0.4 0.4];
% xc1=[0.4 0.5 0.6 0.7 0.8 1 0.8 0.6 0.5 0.3 0.2];
end
xc1=xc1(mi-time_related:mi+time_related);
xc=[xc1 xc1 xc1 xc1 xc1];
C_xr2=C_xr1;
C_xe2=C_xe1;
C_xr2(:,1:(time_related*2+1)*5)=C_xr1(:,1:(time_related*2+1)*5).*xc;
C_xe2(:,1:(time_related*2+1)*5)=C_xe1(:,1:(time_related*2+1)*5).*xc;
end