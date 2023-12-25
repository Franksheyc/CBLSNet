function timefeatures=time_features1(DAT,type)
%%%%提取信号的时域特征  比原本更省时间
switch type
    case 'EEG'
        window=size(DAT,2)/6;
        overlap=0;
        max_n=8192;
    case 'EOG'
        window=size(DAT,2)/6;
        overlap=0;
        max_n=8192;
    case 'EMG'
        window=size(DAT,2)/6;
        overlap=0;
        max_n=8192;
end

[t,s]=size(DAT);
if s>=window*2
    timepiece=(s-overlap)/(window-overlap);
else
    timepiece=1;
end
% zc=zeros(t,1);pe1=zeros(t,1);
    esmeans=[];
    esstds=[];
    esku=[];
    essk=[];
    esEng=[];
    esPe=[];
    esZ=[];
    esmedian=[];
    esPe1=[];
    esLRSSV=[];
for i=1:timepiece    
    dat=DAT(:,(i-1)*window-(i-1)*overlap+1:i*window-(i-1)*overlap);
    means=mean(dat')';
    stds=std(dat')';
    Eng=sqrt(sum(dat'.^2))';
    sk=skewness(dat')';
    ku=kurtosis(dat')';
    med = median(dat')';%每一段中值
    lrssv=log10(sqrt(sum((diff(dat').^2))))';
    pw=abs(dat)'./sum(abs(dat)');
    pe = -sum(pw.*log(pw))';
    z0=dat>0;
    z=diff(z0')';
    parfor j=1:t
        zc(j,:)=length(find(z(j,:)==0))/window;
        pe10 = PermutationEntropy( dat(j,:) , floor(window/8), 4, floor(window/2));
        pe1(j,:)=mean(pe10')'; 
    end
    
    
    esmeans=[esmeans means];
    esstds=[esstds stds];
    esku=[esku ku];
    essk=[essk sk];
    esEng=[esEng Eng];
    esPe=[esPe pe];
    esZ=[esZ zc];
    esmedian=[esmedian med];
    esPe1=[esPe1 pe1];
    esLRSSV=[esLRSSV lrssv];

end
    timefeatures.means=esmeans;
    timefeatures.stds=esstds;
    timefeatures.ku=esku;
    timefeatures.sk=essk;
    timefeatures.Eng=esEng;
    timefeatures.Pe=esPe;
    timefeatures.Z=esZ;
    timefeatures.median=esmedian;
    timefeatures.Pe1=esPe1;
    timefeatures.LRSSV=esLRSSV;