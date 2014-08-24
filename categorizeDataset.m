function [ cat_dataset ] = categorizeDataset(dataset,wt)

[ngenes,nsamples]=size(dataset);
n=100;
fc_dataset=zeros(ngenes,nsamples);

% transform into fold-chage values
for i = 1:nsamples
    fc_dataset(:,i)=log2(dataset(:,i)./mean(wt,2));
end

% get mean and std from randomly selected 100 samples
tmp=[];
for i = randi(nsamples,1,n)
    tmp=[tmp fc_dataset(:,i)];
end
m=mean(tmp,2); s=std(tmp,0,2);

% categorize
weight=1;
cat_dataset=zeros(ngenes,nsamples);
for i = 1:ngenes
    for j = 1:nsamples
        if fc_dataset(i,j) <= m(i)-s(i)*weight
            cat_dataset(i,j)=1;
        elseif m(i)-s(i)*weight < fc_dataset(i,j) && fc_dataset(i,j) < m(i)+s(i)*weight
            cat_dataset(i,j)=2;
        elseif m(i)+s(i)*weight <= fc_dataset(i,j)
            cat_dataset(i,j)=3;
        else
            fprintf('error (%f %f %f) %f %f %f\n',m(i),s(i),weight,m(i)-s(i)*weight,m(i)+s(i)*weight,fc_dataset(i,j));
        end
    end 
end

end

