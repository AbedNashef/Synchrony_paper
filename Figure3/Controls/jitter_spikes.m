function jittered_mat=jitter_spikes(mat,varargin)
jittered_mat=zeros(size(mat));
if ~isempty(varargin)
    n=varargin{1};
    for i=1:size(mat,1)
        ix= find(mat(i,:));
        jitter_sign= randn(1,length(ix));
        jitter_sign=jitter_sign./abs(jitter_sign);
        dt= randi([0 n],[1 length(ix)]) .* jitter_sign;
        ix=ix+dt;
        ix2remove=find(ix<1 | ix>size(mat,2));
        ix(ix2remove)=ix(ix2remove)+-1.*(dt(ix2remove));
        jittered_mat(i,:)=0;
        jittered_mat(i,ix)=1;
    end
else
    for i=1:size(mat,1)
        ix= [0 find(mat(i,:)) size(mat,2)];
        
        for ii=2:length(ix)-1
            this_ix= randi([ceil(ix(ii)-0.5*(ix(ii)-ix(ii-1))) floor(ix(ii)+0.5*(ix(ii+1)-ix(ii)))],1);
            jittered_mat(i,this_ix)=1;
        end
    end
end
end