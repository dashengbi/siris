function data_all = find_seg(ind)

data_seg = [];
if length(ind) > 1
    diff_ind = diff(ind); diff_ind2 = [diff_ind(2:end);0];
    if diff_ind(1) == 1
        data_seg(1,1) = ind(1);
    end

    count = 1;
    for i=1:length(diff_ind)
        if (diff_ind(i) ~= 1) && (diff_ind2(i) == 1)
            data_seg(count,1) = ind(i+1);
        elseif (diff_ind(i) == 1) && (diff_ind2(i) ~= 1) 
            data_seg(count,2) = ind(i+1);
            count = count + 1;
        end
    end

    if diff_ind(end) == 1
        data_seg(end,2) = ind(end);
    end
    
    data_segs = [];
    for iseg = 1 : size(data_seg, 1)
        data_segs = [data_segs, data_seg(iseg,1): data_seg(iseg,2)];
    end
    data_individuals = setdiff(ind, data_segs);
    data_individual = [];
    count_individual = 1;
    for j = 1: numel(data_individuals)
        data_individual(count_individual, 1) = data_individuals(j);
        data_individual(count_individual, 2) = data_individuals(j);
        count_individual = count_individual + 1;
    end
    data_all = sort([data_seg; data_individual], 1); 
end