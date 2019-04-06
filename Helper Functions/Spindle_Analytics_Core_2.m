function stats = Spindle_Analytics_Core_2(standard_labels, algorithm_labels, totduration, starttime)
temp_cell{1} = standard_labels;
standard_labels = d2s3(temp_cell, totduration, starttime, 1, 1);
standard_clone = standard_labels;
temp_cell{1} = algorithm_labels;
algorithm_labels = d2s3(temp_cell, totduration, starttime, 1, 1);
algorithm_clone = algorithm_labels;

spindles{1} = standard_labels;
spindles{2} = algorithm_labels;
tp_labels = d2s3(spindles, totduration, starttime, 2, 2);
stats.ntp = size(tp_labels, 1);

for i = 1:size(tp_labels, 1)
    tp_start = tp_labels(i, 1);
    tp_stop = tp_labels(i, 2);
    for s = 1:size(standard_clone, 1)
        if((tp_start >= standard_clone(s, 1)  && tp_stop <= standard_clone(s, 2)) || (tp_stop >= standard_clone(s, 1) && tp_stop <= standard_clone(s, 2)))
            standard_clone(s, 1) = -1;
        end
    end
    for a = 1:size(algorithm_clone, 1)
        if((tp_start >= algorithm_clone(a, 1)  && tp_stop <= algorithm_clone(a, 2)) || (tp_stop >= algorithm_clone(a, 1) && tp_stop <= algorithm_clone(a, 2)))
            tp_labels(i, 3) = (tp_stop - tp_start) / (algorithm_clone(a, 2) - algorithm_clone(a, 1)) * 100;
            algorithm_clone(a, 1) = -1;
        end
    end
    standard_clone(standard_clone(:, 1) == -1, :) = [];    
    algorithm_clone(algorithm_clone(:, 1) == -1, :) = [];
end

stats.fp_set = algorithm_clone;
stats.nfp = size(stats.fp_set, 1);
stats.fn_set = standard_clone;
stats.nfn = size(stats.fn_set, 1);

stats.algotot = size(algorithm_labels, 1);
stats.stndtot = size(standard_labels, 1);

stats.fpr = stats.nfp / stats.algotot;
stats.fnr = stats.nfn / stats.stndtot;

stats.tpr = stats.ntp / stats.algotot;
stats.recall = stats.ntp / stats.stndtot;
stats.precision = stats.ntp / stats.algotot;
stats.f1 = harmmean([stats.recall, stats.precision]);

stats.overlap = tp_labels;

if isnan(stats.f1)
    stats.f1 = 0;
end

