%@INPUT
% {PARAMETERS}:
%   Window length = WINLEN
%   Overlap = OVERLAP
%   Freq band = [FBL, FBH]
%   Relative Total Power > LT for period
%   Relative Total Power > UT once
%   Tolerance = TOLERANCE
% DATA_N;
% DATA_F_ABS;
% FREQS;
% IND_NREM;
%
%@OUTPUT
% SPINDLE_POINTS;

function spindle_points = Spindle_STFT_Find(parameters, data_n, data_f_abs, freqs, ind_NREM)
    parameters = strsplit(parameters);
    if length(parameters) ~= 7
        error('# parameters incorrect. Input format: WINLEN OVERLAP FREQBANDLO FREQBANDHI THLO THHI TOLERANCE');
    else
        winlen = str2double(parameters{1});
        overlap = str2double(parameters{2});
        fbl = str2double(parameters{3});
        fbh = str2double(parameters{4});
        lt = str2double(parameters{5});
        ut = str2double(parameters{6});
        tolerance = str2double(parameters{7});
    end
    
    %% Find signal power ratio
    rel_pwr = cell(1,size(ind_NREM,1));
    total_pwr = 0;
    count = 0;
    for i=1:size(ind_NREM, 1)
        total_pwr = total_pwr + sum(sum(cell2mat(data_f_abs{i})));
        count = count + size(cell2mat(data_f_abs{i}), 2);
    end
    total_pwr_mean = total_pwr / count;

    freq_low = find(freqs >= fbl, 1, 'first'); % Lower freq bound
    freq_high = find(freqs <= fbh, 1, 'last'); % Upper freq bound
    
    %% Form logical array
    clear spindles_bool;
    for i=1:size(ind_NREM, 1)
        spindles_bool{i} = [];
        for j=1:size(data_f_abs{i}, 2)
            rel_pwr{i}(j) = sum(data_f_abs{i}{j}(freq_low:freq_high)) / total_pwr_mean;
            if rel_pwr{i}(j) > ut
                spindles_bool{i} = cat(1, spindles_bool{i}, 2); % Exceed upper threshold
            elseif rel_pwr{i}(j) > lt
                spindles_bool{i} = cat(1, spindles_bool{i}, 1); % Exceed lower threshold
            else
                spindles_bool{i} = cat(1, spindles_bool{i}, 0); % Lower than threshold
            end
        end
    end

    %% Count spindles
    spindle_count = []; % Total # of spindles, array for different segments of NREM sleep
    clear spindles_index spindle_count;
    nonoverlap = winlen - overlap;
    for i=1:size(spindles_bool, 2)
        spindles_index{i} = [];
        spindle_count(i) = 0;
        j = 1;
        while (j <= size(spindles_bool{i}, 1))
            if (spindles_bool{i}(j) == 1) || (spindles_bool{i}(j) == 2)
                spindle_count(i) = spindle_count(i) + 1;
                spindles_index{i}(spindle_count(i), 1) = (j - 1) * nonoverlap + 1;
                spindle_start_index = j;
                while true     
                    j = j + 1;
                    %Check BOTH tolerances
                    if (j > size(spindles_bool{i}, 1)) || (isempty(find(spindles_bool{i}(j:min(j + tolerance, size(spindles_bool{i}, 1))) == 1, 1)) && isempty(find(spindles_bool{i}(j:min(j + tolerance, size(spindles_bool{i}, 1))) == 2, 1)))
                        spindles_index{i}(spindle_count(i), 2) = (j - 1) * nonoverlap + winlen; % Window "offset"
                        spindle_end_index = j - 1;

                        %Check for duration + upper threshold                    
                        if (isempty(find(spindles_bool{i}(spindle_start_index:spindle_end_index) == 2, 1)))
                            spindles_index{i}(spindle_count(i), :) = [];
                            spindle_count(i) = spindle_count(i) - 1;
                        else
                            if (spindles_index{i}(spindle_count(i), 2) - spindles_index{i}(spindle_count(i), 1) < 500) || (spindles_index{i}(spindle_count(i), 2) - spindles_index{i}(spindle_count(i), 1) > 3000)
                                spindles_index{i}(spindle_count(i), :) = [];
                                spindle_count(i) = spindle_count(i) - 1;
                            end
                        end

                        break;
                    end
                end
            end
            j = j + 1;
        end
    end

    spindle_points = convert_labels(spindles_index, data_n, length(cell2mat(data_n')));
    spindle_cell{1} = spindle_points;
    spindle_points = d2s2(spindle_cell, length(cell2mat(data_n')), 1, 1);

end