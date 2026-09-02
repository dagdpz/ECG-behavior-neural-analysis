function p= ecg_bna_population_significant_n_bins(in)

in=abs(in);
bins1=16:26;
bins2=37:47;
for f=1:size(in,1)
    n_before=sum(squeeze(in(f,bins1,:)),1);
    n_after=sum(squeeze(in(f,bins2,:)),1);
%     %pu=ranksum(n_after,n_before); %unpaired
     p(f)=signrank(n_after-n_before); %paired
    
end


end