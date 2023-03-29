function metrics = sumperf(evalmetrics)

    numOfSample = length([evalmetrics.accuracy]);
    mean_perf = mean([evalmetrics.perf]).*100;
    mean_accuracy = mean([evalmetrics.accuracy]).*100;
    std_perf = std([evalmetrics.accuracy]).*100;
    mean_missclass = mean([evalmetrics.missclass]).*100;
    mean_recall_c1 = mean([evalmetrics.recall_c1]).*100;
    mean_recall_c2 = mean([evalmetrics.recall_c2]).*100;
    mean_precision_c2 = mean([evalmetrics.precision_c2]).*100;
    mean_precision_c1 = mean([evalmetrics.precision_c1]).*100;
    mean_f1 = mean([evalmetrics.f1]).*100;
    mean_tpr = mean([evalmetrics.tpr]).*100;
    mean_fpr = mean([evalmetrics.fpr]).*100;
    mean_interval = mean([evalmetrics.interval]).*100;
    
    se = std_perf ./ sqrt(numOfSample);
    ts = tinv( [0.025 0.975] , numOfSample-1);
    CI = mean_perf + ts.*se;

    metrics.perf = mean_perf;
    metrics.accuracy = mean_accuracy;
    metrics.perf_std = std_perf;
    metrics.missclass = mean_missclass;
    metrics.recall_c1 = mean_recall_c1;
    metrics.recall_c2 = mean_recall_c2;
    metrics.precision_c2 = mean_precision_c2;
    metrics.precision_c1 = mean_precision_c1;
    metrics.f1 = mean_f1;
    metrics.tpr = mean_tpr;
    metrics.fpr = mean_fpr;
    metrics.interval = mean_interval;
    metrics.std_error = se;
    metrics.CI = CI;

    
%     interval = 1.96.*sqrt((accuracy.*(1-accuracy))./length(88)); % CI 95%

    
    
    
    
end