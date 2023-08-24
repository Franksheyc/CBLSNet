function EOGcof=EOG_cof(dat,dat1,samplingrate)
dat(dat==0)=eps;
dat1(dat1==0)=eps;
path1 = [0.25 1.2] / samplingrate * 2;%SW
path2 = [1.2 35] / samplingrate * 2;
[p1_B,p1_A] = butter( 4 , path1 , 'bandpass' ); 
[p2_B,p2_A] = butter( 4 , path2 , 'bandpass' ); 
for i=1:size(dat,1)
    EOG1_low = filtfilt( p1_B , p1_A , dat(i,:)' );
    EOG1_high = filtfilt( p2_B , p2_A , dat(i,:)' );
    EOG2_low = filtfilt( p1_B , p1_A , dat1(i,:)' );
    EOG2_high = filtfilt( p2_B , p2_A , dat1(i,:)' );
    a=angle(hilbert(EOG1_low));
    b=angle(hilbert(EOG2_low));
    
    a1=angle(hilbert(EOG1_high));
    b1=angle(hilbert(EOG2_high));
    
    angcof_l2l=min(min(corrcoef(a,b)));
    angcof_h2h=min(min((corrcoef(a1,b1))));
    angcof_l2h=min(min(corrcoef(a,b1)));
    angcof_h2l=min(min((corrcoef(a1,b))));
%     
%     dat0=dat(i,:)+dat1(i,:);
%     b=abs(dat0)/abs(dat(i,:));    
    angcof=[angcof_l2l angcof_h2h angcof_l2h angcof_h2l];
    EOGcof(i,:)=angcof;
end
end
