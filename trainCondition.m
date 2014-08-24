function [ best_model, features, train_accuracy, parameter_space, best_r, all_features] = trainCondition(dataset, train_answer, algorithm)

ngenes=size(dataset,1);

% measure mutual information for each gene with respect to class label
mi=zeros(ngenes,2);
j=1;
for x = 1:ngenes
    mi(j,1)=x;
    mi(j,2)=mutualInformation(dataset(x,:)',train_answer);
    j=j+1;
end
mi=sortrows(mi,-2);
best_accuracy=0;
best_parameter=0;
best_model=0;
best_r=zeros(size(train_answer,1));

indices = crossvalind('Kfold',train_answer,2); 
test_indices = (indices == 1);
train_indices = ~test_indices;

if strcmp(algorithm,'SVM') % Multiclass SVM (O)
    c_range=[0.5 1 2 4];
    parameter_range=[10 100 200 400];
    parameter_space=zeros(size(parameter_range,2),size(c_range,2),3);
    for i = 1:size(parameter_range,2)
        for j = 1:size(c_range,2)
            
            r=multisvm(dataset(mi(1:parameter_range(i),1),train_indices)',train_answer(train_indices),dataset(mi(1:parameter_range(i),1),test_indices)',c_range(j));  
            model=c_range(j);
            
            cm=confusionmat(train_answer(test_indices),r);
            accuracy=mean(diag(cm)./sum(cm,2));
    
            if best_accuracy < accuracy
                best_accuracy=accuracy;
                best_parameter=parameter_range(i);
                best_model=model;
                best_r=r;
            end
            parameter_space(i,j,1)=parameter_range(i);
            parameter_space(i,j,2)=c_range(j);
            parameter_space(i,j,3)=accuracy;
            fprintf('[Training SVM] NG: %d C: %d Accuracy: %f\n',parameter_range(i),c_range(j),accuracy);
        end
    end
else
    parameter_range=[10 100:300:ngenes ngenes];
    parameter_space=zeros(size(parameter_range,2),2);
    for i = 1:size(parameter_range,2)   
        if strcmp(algorithm,'NB') % NB (O)
            model=NaiveBayes.fit(dataset(mi(1:parameter_range(i),1),train_indices)',train_answer(train_indices),'Dist','mvmn');
            r=model.predict(dataset(mi(1:parameter_range(i),1),test_indices)');   
        elseif strcmp(algorithm,'DT') % Decision Tree (O)
            model=ClassificationTree.fit(dataset(mi(1:parameter_range(i),1),train_indices)',train_answer(train_indices));
            r=predict(model,dataset(mi(1:parameter_range(i),1),test_indices)');
        elseif strcmp(algorithm,'KNN') % KNN (O)
            model=ClassificationKNN.fit(dataset(mi(1:parameter_range(i),1),train_indices)',train_answer(train_indices));
            r=predict(model,dataset(mi(1:parameter_range(i),1),test_indices)');
        else
            fprintf('[error] unknown algorithm\n');
            break;
        end
        cm=confusionmat(train_answer(test_indices),r);
        accuracy=mean(diag(cm)./sum(cm,2));
    
        if best_accuracy < accuracy
            best_accuracy=accuracy;
            best_parameter=parameter_range(i);
            best_model=model;
            best_r=r;
        end
        parameter_space(i,1)=parameter_range(i);
        parameter_space(i,2)=accuracy;
        fprintf('[%s] feature size %d, accuracy %f\n',algorithm,parameter_range(i), accuracy);
    end
    
end     
        
features=mi(1:best_parameter,:);
all_features=mi;
train_accuracy=best_accuracy;
%fprintf('selected feature size: %d (accuracy %f)\n',best_parameter,best_accuracy);

end

