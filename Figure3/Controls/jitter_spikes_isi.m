function X=jitter_spikes_isi(X)
new_x=zeros(size(X));
for i=1:size(X,1)
    this_trial=X(i,:);
    INDEX= [find(this_trial,1,'first') diff(find(this_trial))];
    [~,ix]= sort(rand(1,length(INDEX)));
    INDEX=INDEX(ix);

    new_x(i,cumsum(INDEX))=1;
end
X=new_x;
end