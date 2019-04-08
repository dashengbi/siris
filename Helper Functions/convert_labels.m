function spindle_labels = convert_labels(spindle_seg, data_spindle, maxtime)
spindle_labels = [];
for i=1:length(data_spindle)
    offset = size(cat(1, data_spindle{1:i-1}), 1);
    for j=1:size(spindle_seg{i})
        start = spindle_seg{i}(j, 1) + offset;
        stop = spindle_seg{i}(j, 2) + offset;
        if stop > maxtime
            return;
        end
        spindle_labels = cat(1, spindle_labels, [start stop]);
    end
end