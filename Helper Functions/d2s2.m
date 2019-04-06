function spindle_labels = d2s2(spindle_cell, totduration, llevel, ulevel) % Function for marking spindle time indices
spindle_labels = [];
j = 1;
spindle_bool = zeros(1, totduration);
for i=1:length(spindle_cell)
    temp = zeros(1, totduration);
    for k=1:size(spindle_cell{i}, 1)
        if spindle_cell{i}(k, 2) <= totduration
            temp(round(spindle_cell{i}(k, 1)):round(spindle_cell{i}(k, 2))) = 1;
        else
            continue;
        end
    end
    spindle_bool = spindle_bool + temp;
end
while(j <= min(totduration, length(spindle_bool)))
    if (spindle_bool(j) >= llevel) && (spindle_bool(j) <= ulevel)
        start_ind = j;
        while (j <= min(totduration, length(spindle_bool))) && (spindle_bool(j) >= llevel) && (spindle_bool(j) <= ulevel)
            j = j + 1;
        end
        end_ind = j - 1;
        if(end_ind - start_ind <= 3000) %Condition for duration
            spindle_labels = cat(1, spindle_labels, [start_ind, end_ind]);
        end
    end
    j = j + 1;
    
end