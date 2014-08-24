function [train_accuracy, test_accuracy, parameter_spaces, post_vector, r_vector, totalcm, mi] = crossValidateCondition(dataset,pheno,nfold,nclasses,algorithm)


indices = crossvalind('Kfold',pheno,nfold); 
ngenes=size(dataset,1);
nparameters=size(100:300:ngenes,2)+2;

fprintf('cross-validation ready\n');
if strcmp(algorithm, 'All')
    train_accuracy=zeros(5,nfold,1);
    test_accuracy=zeros(5,nfold,1);
    mi=zeros(4,nfold,ngenes,2);
else
    train_accuracy=zeros(nfold,1);
    test_accuracy=zeros(nfold,1);
    mi=zeros(nfold,ngenes,2);
end


parameter_spaces=zeros(nfold,nparameters,2);
post_vector=zeros(size(dataset,2),nclasses);
r_vector=zeros(size(dataset,2),1);
totalcm=zeros(nclasses,nclasses);    

for i = 1:nfold 
    test_indices = (indices == i);
    train_indices = ~test_indices;
    
    if strcmp(algorithm, 'All')
        [svm_model,svm_features,~,~,svm_train_r,svm_mi]=trainCondition(dataset(:,train_indices),pheno(train_indices),'SVM');
        [nb_model,nb_features,~,~,nb_train_r,nb_mi]=trainCondition(dataset(:,train_indices),pheno(train_indices),'NB');
        [dt_model,dt_features,~,~,dt_train_r,dt_mi]=trainCondition(dataset(:,train_indices),pheno(train_indices),'DT');
        [knn_model,knn_features,~,~,knn_train_r,knn_mi]=trainCondition(dataset(:,train_indices),pheno(train_indices),'KNN');
        
        mi(1,i,:,1)=svm_mi(:,1); mi(1,i,:,2)=svm_mi(:,2);
        mi(2,i,:,1)=nb_mi(:,1); mi(2,i,:,2)=nb_mi(:,2);
        mi(3,i,:,1)=dt_mi(:,1); mi(3,i,:,2)=dt_mi(:,2);
        mi(4,i,:,1)=knn_mi(:,1); mi(4,i,:,2)=knn_mi(:,2);
        
        nb_test_r=nb_model.predict(dataset(nb_features(:,1),test_indices)');
        svm_test_r=multisvm(dataset(svm_features(:,1),train_indices)',pheno(train_indices),dataset(svm_features(:,1),test_indices)',svm_model);
        dt_test_r=predict(dt_model,dataset(dt_features(:,1),test_indices)');
        knn_test_r=predict(knn_model,dataset(knn_features(:,1),test_indices)');
        
        nb_train_r=nb_model.predict(dataset(nb_features(:,1),train_indices)');
        svm_train_r=multisvm(dataset(svm_features(:,1),train_indices)',pheno(train_indices),dataset(svm_features(:,1),train_indices)',svm_model);
        dt_train_r=predict(dt_model,dataset(dt_features(:,1),train_indices)');
        knn_train_r=predict(knn_model,dataset(knn_features(:,1),train_indices)');
        
        total_train_r=mode([nb_train_r svm_train_r dt_train_r knn_train_r],2);
        cm=confusionmat(pheno(train_indices),total_train_r);
        train_accuracy(5,i)=mean(diag(cm)./sum(cm,2));
        
        nb_train_cm=confusionmat(pheno(train_indices),nb_train_r);
        svm_train_cm=confusionmat(pheno(train_indices),svm_train_r);
        dt_train_cm=confusionmat(pheno(train_indices),dt_train_r);
        knn_train_cm=confusionmat(pheno(train_indices),knn_train_r);
        
        train_accuracy(1,i)=mean(diag(nb_train_cm)./sum(nb_train_cm,2));
        train_accuracy(2,i)=mean(diag(svm_train_cm)./sum(svm_train_cm,2));
        train_accuracy(3,i)=mean(diag(dt_train_cm)./sum(dt_train_cm,2));
        train_accuracy(4,i)=mean(diag(knn_train_cm)./sum(knn_train_cm,2));
        
        fprintf('[Fold %d] NB training accuracy %f\n', i, train_accuracy(1,i));
        fprintf('[Fold %d] SVM training accuracy %f\n', i, train_accuracy(2,i));
        fprintf('[Fold %d] DT training accuracy %f\n', i, train_accuracy(3,i));
        fprintf('[Fold %d] KNN training accuracy %f\n', i, train_accuracy(4,i));
        fprintf('[Fold %d] Consensus training accuracy %f\n', i, train_accuracy(5,i));
        
        nb_test_cm=confusionmat(pheno(test_indices),nb_test_r);
        svm_test_cm=confusionmat(pheno(test_indices),svm_test_r);
        dt_test_cm=confusionmat(pheno(test_indices),dt_test_r);
        knn_test_cm=confusionmat(pheno(test_indices),knn_test_r);
        
        test_accuracy(1,i)=mean(diag(nb_test_cm)./sum(nb_test_cm,2));
        test_accuracy(2,i)=mean(diag(svm_test_cm)./sum(svm_test_cm,2));
        test_accuracy(3,i)=mean(diag(dt_test_cm)./sum(dt_test_cm,2));
        test_accuracy(4,i)=mean(diag(knn_test_cm)./sum(knn_test_cm,2));
        
        [r,f]=mode([nb_test_r svm_test_r dt_test_r knn_test_r],2);        
    else
        [model,features,accuracy,parameter_space,~,all_features]=trainCondition(dataset(:,train_indices),pheno(train_indices),algorithm);
        train_accuracy(i)=accuracy;
        if strcmp(algorithm, 'NB') % NB
            [post,r]=model.posterior(dataset(features(:,1),test_indices)');
            mi(i,:,1)=all_features(:,1);
            mi(i,:,2)=all_features(:,2);
        elseif strcmp(algorithm, 'SVM') % Multiclass SVM
            r=multisvm(dataset(features(:,1),train_indices)',pheno(train_indices),dataset(features(:,1),test_indices)',model);
        elseif strcmp(algorithm, 'DT') % Decision Tree
            r=predict(model,dataset(features(:,1),test_indices)');
        elseif strcmp(algorithm, 'KNN') % KNN
            r=predict(model,dataset(features(:,1),test_indices)');
        else
            fprintf('unknown algorithm\n');
            break;
        end
    end
        
    cm=confusionmat(pheno(test_indices),r);
    totalcm=totalcm+cm;
    r_vector(test_indices)=r;
        
    if strcmp(algorithm,'NB')
        %for j = 1:nclasses
        %    post_vector(test_indices,j)=post(:,j);
        %end
        
        for j = 1:nparameters
            parameter_spaces(i,j,1)=parameter_space(j,1);
            parameter_spaces(i,j,2)=parameter_space(j,2);
        end
        test_accuracy(i)=mean(diag(cm)./sum(cm,2));
    end    
    if strcmp(algorithm,'All')
        post_vector(test_indices)=f/nclasses;
        test_accuracy(5,i)=mean(diag(cm)./sum(cm,2));
    end
        
        
    fprintf('test accuracy : %f\n',test_accuracy(i));
end

end

