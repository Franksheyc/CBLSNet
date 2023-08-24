function [N1,N2,fe,fe1]=trainpar(b,tur,kk,ten,time_related)
n=6;
b1=ceil(b/n);
unit=time_related*2+1;
if tur==1
%     N1=[ceil(b1/12)  ceil(b1/6)  ceil(b1/4) ceil(b1/3)  ceil(b1/12) ];
    N1=[1 ceil(b1/12)  ceil(b1/6)  ceil(b1/6) ceil(b1/6)  ceil(b1/3)];
%     N1=[1 ceil(b1/12)  ceil(b1/6)  ceil(b1/6) ceil(b1/6)  ceil(b1/6) ceil(b1/6)];
    N2=[1 2];
    if kk>2
%         b1=b1*n;
          fe={[1:6 7:12];[33:62];[63:68 69:74];[b+1:b+ten]};
%         7:12];[62:67];[68:74];[220:228];[b+1:b+ten]};%[220:228]  acc85.86
%         fe1={time_related+1+unit*([1:5]-1)};
        fe1={[6 17 28 39 50]};
    else
%         fe={[1:6 7:12];[62:67];[68:74];[220:228]};
        fe={[1:6 7:12];[33:62];[63:68 69:74]};
        fe1={};
    end
%     fe={};
%     fe1={};
else
%     b=ceil(b/n);
    N1=[1 ceil(b1/12)  ceil(b1/6)  ceil(b1/6) ceil(b1/6)  ceil(b1/3)];
    N2=[4];%[1 1 2]
    fe={1:5*(time_related*2+1)};
    fe1={1:unit*5};
% fe={};
% fe1={};
end