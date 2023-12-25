multi
ty1=0;ty2=0;ty3=0;tt1=0;tt2=0;
for i=1:20
for j=1:3
ty1=ty1+time_reset(i,j);
ty2=ty2+time_te_f(i,j);
ty3=ty3+time_tr_f(i,j);
end
end
for i=1:20
for j=1:2
for k=1:3
tt1=tt1+time_te_t(i,j,k);
tt2=tt2+time_tr_t(i,j,k);
end
end
end
time_tr=ty1+ty3+tt2;
time_te=ty2+tt1;
time_tr/20
time_te/20
time_te
time_tr


single
ty1=0;ty2=0;ty3=0;tt1=0;tt2=0;
for i=1:20
for j=1:3
ty1=ty1+time_reset(i,j);
ty2=ty2+time_te_f(i,j);
ty3=ty3+time_tr_f(i,j);
end
end
for i=1:20
for j=1:2
ty2=ty2+time_te_t(i,j);
ty3=ty3+time_tr_t(i,j);
end
end
time_tr=ty1+ty3+tt2;
time_te=ty2+tt1;
time_tr/20
time_te/20
time_te
time_tr