%for binary classification task
function fscore=fscore(label,prediction)

    % INPUT
% group = true class labels
% grouphat = predicted class labels
%
% OR INPUT
% stats = confusionmatStats(group);
% group = confusion matrix from matlab function (confusionmat)
%
% OUTPUT
%               Predicted Classes
%                    p'    n'
%              ___|_____|_____|
%       Actual  p |     |     |
%      Classes  n |     |     |
%
% stats.accuracy = (TP + TN)/(TP + FP + FN + TN) ; the average accuracy is returned
% stats.precision = TP / (TP + FP)                  % for each class label
% stats.sensitivity = TP / (TP + FN)                % for each class label
% stats.specificity = TN / (FP + TN)                % for each class label
% stats.recall = sensitivity                        % for each class label
% stats.Fscore = 2*TP /(2*TP + FP + FN)            % for each class label
%
% TP: true positive, TN: true negative,
% FP: false positive, FN: false negative
%

if nargin < 2
    value1 = label;
else
    [value1,gorder] = confusionmat(label,prediction);
end

numOfClasses = size(value1,1);
totalSamples = sum(sum(value1));

[TP,TN,FP,FN,accuracy,sensitivity,specificity,precision,f_score] = deal(zeros(numOfClasses,1));
for class = 1:numOfClasses
   TP(class) = value1(class,class);
   tempMat = value1;
   tempMat(:,class) = []; % remove column
   tempMat(class,:) = []; % remove row
   TN(class) = sum(sum(tempMat));
   FP(class) = sum(value1(:,class))-TP(class);
   FN(class) = sum(value1(class,:))-TP(class);
end

if numOfClasses > 2
    for class = 1:numOfClasses
        fscores.accuracy(class) = (TP(class) + TN(class)) / totalSamples;
        fscores.sensitivity(class) = TP(class) / (TP(class) + FN(class));
        fscores.specificity(class) = TN(class) / (FP(class) + TN(class));
        fscores.falsepositiverate(class) = FP(class)/ (FP(class)+ TN(class));
        fscores.precision(class) = TP(class) / (TP(class) + FP(class));
        fscores.f_score(class) = 2*TP(class)/(2*TP(class) + FP(class) + FN(class));
    end
else
    fscores.accuracy = (value1(1,1) + value1(2,2)) / totalSamples; %TP+TN/all
    fscores.sensitivity = value1(1,1) / (value1(1,1) + value1(2,1)); % TP/(TP+FN)
    fscores.specificity = value1(2,2) / (value1(1,2) + value1(2,2)); %TN/(FN+TN)
    fscores.falsepositiverate = value1(1,2)/ (value1(1,2)+ value1(2,2)); %FP/(FP+TN)
    fscores.precision= value1(1,1) / (value1(1,1) + value1(1,2)); %TP/(TP+FP)
    fscores.f_score = 2*value1(1,1)/(2*value1(1,1) + value1(1,2) + value1(2,1)); %2*TP/(2*TP+FP+FN)
end
fscore = fscores;
end
