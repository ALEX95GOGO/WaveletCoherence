function [Objective]  = OptimizationWCOHelper(ts_x1, ts_x2, labels_vec, x)


  [Objective, Taskperformance, Cardiac_1, Cardiac_2] = WCOExperiment(ts_x1, ts_x2, labels_vec, x);

end
