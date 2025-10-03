function    [widths,vwidth]= get_width(anti_diag)
gwin=gausswin(51); gwin=gwin./sum(gwin);
ix= size(anti_diag,2)/2;
vwidth=zeros(size(anti_diag,1),1);
for i=1:size(anti_diag,1)
    tmp= conv2(anti_diag(i,:),gwin','same');
    tmp=tmp-mean(tmp(1:100));
    if tmp(ix)/2<0
        tmp=tmp.*-1;
    end
    index= find(tmp(1:ix-1)<tmp(ix)/2,1,'last');
    try
        w1=interp1(tmp(index:ix),index:ix,tmp(ix)/2);
    catch me
        w1=nan;
    end
    index= ix+find(tmp(ix+1:end)<tmp(ix)/2,1,'first');
    try
        w2=interp1(tmp(ix:index),ix:index,tmp(ix)/2);
    catch me
        w2=nan;
    end
    widths(i)=(w2-w1)/2;
    if abs(tmp(ix))>2*std(tmp)
        vwidth(i)=1;
    end
end
end