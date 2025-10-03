function [anti_diag,all_diag]= get_anti_diag(mat,ax,time0)
ix= find(round(ax*1000)==round(time0*1000));
xix=1:length(ax)+time0*1000;
yix= length(ax)+time0*1000:-1:1;
index= find(xix>0);
xix=xix(index);
yix=yix(index);
if size(mat,3)>1
    for i=1:size(mat,1)
        mat(i,:,:)=imgaussfilt(squeeze(mat(i,:,:)),2);
        for ii=1:length(xix)
            all_diag(i,ii)=diag(squeeze(mat(i,xix(ii),yix(ii))));
            anti_diag(ii)=ax(end)-xix(ii)./1000;

        end
    end
else % just one matrix
    mat=imgaussfilt(squeeze(mat),2);
    for ii=1:length(xix)
        all_diag(ii)=diag(squeeze(mat(xix(ii),yix(ii))));
        anti_diag(ii)=ax(end)-xix(ii)./1000;

    end
end
end