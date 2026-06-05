function ccg=get_ccg(mat)
ax=-10:10;
ccg=[];
for i=1:size(mat,1)
    for ii=i+1:size(mat,1)
        x1= mat(i,:);
        x2=mat(ii,:);
        nx2=get_null_fr(x2,1);
        ix=find(x1==1);
        ix=ix(find(ix>abs(ax(1)) & ix<size(mat,2)-ax(end)));
        index=repmat(ix(:),1,length(ax))+repmat(ax,length(ix),1);
        tmp=mean(x2(index)-nx2(index),1);
        ccg=[ccg; tmp(find(ax==0))];

    end
end
% ccg=nanmean(ccg);
end