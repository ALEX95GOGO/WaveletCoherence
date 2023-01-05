function PM=PM(labels, predictors, varargin)
  %% Performance Measures (PM)

    % calculates different task performance measures see OUTPUT
    % the relevant values for the calculation can be reduced with the space of interest and cone of influences
    % this enables a more precise estimation

  %% INPUT

    % labels: vector representing the "ground truth" in time
    % predictors: vector or matrix containing the predicted values, must have the same type as the label
    % ts: sampling rate in [s]
    % coi: the cone of influence is a vector where performance measures are not valid due to edge effects of previous Analysis
    % soi: is a two valued vector describing the region of interest in the frequency space, where the performance should be measured
    % space: vector representing the frequency space in scale, periods or HZ, should correspond to coi

  %% OUTPUT

    % pm: a struct containing accuracy, sensitivity, specificity, false-positive-rate, precision, F-score

  %% REQUIREMENTS

    % fscore which requires confusionmat from Matlab

  Args=struct('ts', 1/10,...  % sampling rate
              'coi', [],...  % the cone of influence
              'soi', [],...  % region of interest
              'space',[]...  % frequency domain (can be scale or period)
              );

  Args=parseargs_special(varargin,Args);

  %% consider cone of influence
  % currently the coi is not supported
  % instead the first and last 30s are cutted out to avoid edge effects
  labels = labels(30/Args.ts:end-30/Args.ts);  % cut start and end
  predictors = predictors(:,30/Args.ts:end-30/Args.ts);  % cut start and end

  %% reduce to space / frequency of interest
  if any(Args.soi) & any(Args.space)
   idx = find(Args.space >= Args.soi(1) & Args.space <= Args.soi(2));
   labels = repmat(labels', 1, length(idx));  % repeat vector of ground truth to create matrix of labels
   predictors = predictors(idx,:);
  end

  PM = fscore(labels(:),predictors(:));
end
