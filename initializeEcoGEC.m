conds=dlmread('ecogec_v1.conditions_v3.txt');
partition=dlmread('ecogec_v1.partition.txt');
total_dataset=dlmread('ecogec_v1.dataset.adj.txt');

[ngenes,nsamples]=size(total_dataset);

%grn=loadGRN('GRN_RDB+Costello.txt',ngenes); 

WT=0;
DATA=1;
TESTABLE=2;

GPL199=0;
GPL3154=1;
RNASeq=2;

new_dataset=zeros(ngenes,nsamples);

GPL199_wt=total_dataset(:,partition(:,2)==GPL199 & partition(:,1)==WT);
GPL3154_wt=total_dataset(:,partition(:,2)==GPL3154 & partition(:,1)==WT);
RNASeq_wt=total_dataset(:,partition(:,2)==RNASeq & partition(:,1)==WT);

new_dataset(:,partition(:,2)==GPL199 & partition(:,1)~=WT)=categorizeDataset(total_dataset(:,partition(:,2)==GPL199 & partition(:,1)~=WT),GPL199_wt);
new_dataset(:,partition(:,2)==GPL3154 & partition(:,1)~=WT)=categorizeDataset(total_dataset(:,partition(:,2)==GPL3154 & partition(:,1)~=WT),GPL3154_wt);
new_dataset(:,partition(:,2)==RNASeq & partition(:,1)~=WT)=categorizeDataset(total_dataset(:,partition(:,2)==RNASeq & partition(:,1)~=WT),RNASeq_wt);

%GPL_wt=total_dataset(:,(partition(:,2)==GPL199 | partition(:,2)==GPL3154 | partition(:,2)==RNASeq) & partition(:,1)==WT);
%new_dataset(:,(partition(:,2)==GPL3154 | partition(:,2)==GPL199 | partition(:,2)==RNASeq) & partition(:,1)~=WT)=categorizeDataset(total_dataset(:,(partition(:,2)==GPL3154 | partition(:,2)==GPL199 |  partition(:,2)==RNASeq) & partition(:,1)~=WT),GPL_wt);


new_dataset=new_dataset(:,partition(:,1)~=WT);
new_conds=conds(conds(:,1)~=WT,:);

