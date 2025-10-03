function X=montecarlo_jitter_spikes(X,bin_size)
jitter_tail= mod(size(X,2),bin_size);
addon= X(:,end-jitter_tail+1:end);
X=X(:,1:end-jitter_tail);
for i=1:size(X,1)
    this_trial=X(i,:);
    win_trial= reshape(this_trial,bin_size,[]);
    new_trial=zeros(size(win_trial));
    for ii=1:size(win_trial,2)
        ix= find(win_trial(:,ii));
        if ~isempty(ix)
            for j=1:length(ix)
                rix= randi(size(new_trial,1),1);
                new_trial(rix,ii)=1;
            end
        end
    end
    new_X(i,:)=new_trial(:)';
    
end

X=[new_X addon];
end