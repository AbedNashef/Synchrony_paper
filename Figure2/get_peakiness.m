function    [peakiness,vpeak]= get_peakiness(anti_diag,ax)
ix= size(anti_diag,2)/2;
bix= [ix-30 ix+30];
gwin=gausswin(51); gwin=gwin./sum(gwin);
for i=1:size(anti_diag,1)
    tmp= conv2(anti_diag(i,:),gwin','same');
    tmp=(tmp-mean(tmp))./std(tmp);
    v=1;
    % if tmp(ix)/2<0
    %     tmp=tmp.*-1;
    %     v=-1;
    % end
    peakiness(i)= tmp(ix)/(abs(tmp(ix))+abs(mean(tmp(bix))));%(tmp(ix)-(mean(tmp(bix))))/(tmp(ix)+(mean((tmp(bix)))));%
    peakiness(i)=peakiness(i)*v;
    if abs(tmp(ix))>2*std(tmp)
        vpeak(i)=1;
    else
        vpeak(i)=0;
    end
    vpeak(i)=1-length(find(tmp<tmp(ix)))/length(tmp);
end
end