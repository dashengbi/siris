function [data_n_cut, ind_NREM_cut] = Spindle_STFT_Cut(data_n, ind_NREM, time)
offset = 0;
for seg=1:size(ind_NREM, 1)
    offset = offset + length(data_n{seg});
    if (offset >= time)
        break;
    end
end
data_n_cut = data_n(1:seg);
ind_NREM_cut = ind_NREM(1:seg, :);
end