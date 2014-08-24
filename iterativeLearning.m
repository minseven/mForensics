function [pheno] = iterativeLearning(cond,dataset,rcond)

pheno=cond;

% Initialize dataset
train_indices=(cond~=4);
test_indices=~train_indices;

while (true)
    [svm_model,svm_features,~,~,~,~]=trainCondition(dataset(:,train_indices),pheno(train_indices),'SVM');
    [nb_model,nb_features,~,~,~,~]=trainCondition(dataset(:,train_indices),pheno(train_indices),'NB');
    [dt_model,dt_features,~,~,~,~]=trainCondition(dataset(:,train_indices),pheno(train_indices),'DT');
    [knn_model,knn_features,~,~,~,~]=trainCondition(dataset(:,train_indices),pheno(train_indices),'KNN');
    
    svm_r=multisvm(dataset(svm_features(:,1),train_indices)',pheno(train_indices),dataset(svm_features(:,1),test_indices)',svm_model);
    [~,nb_r]=nb_model.posterior(dataset(nb_features(:,1),test_indices)');
    dt_r=predict(dt_model,dataset(dt_features(:,1),test_indices)');
    knn_r=predict(knn_model,dataset(knn_features(:,1),test_indices)');
          
    [consensus_r,F,~]=mode([nb_r svm_r dt_r knn_r],2);
    fprintf('%d/%d\n',sum(F>=3),length(F));
    
    predicted=consensus_r;
    histc(predicted,[1,2,3])
    predicted(F<3)=4;
    %fprintf('similarity %f\n',sum(diag(confusionmat(pheno(test_indices),consensus_r))));
    confusionmat(rcond(test_indices),consensus_r);
    pheno(test_indices)=predicted;
   
    train_indices=(pheno~=4);
    test_indices=~train_indices;
    
    if sum(F<3)/length(F) > 0.2
       break; 
    end
end