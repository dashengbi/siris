%@INPUT
% {PARAMETERS}:
%   Window length = WINLEN
%   Overlap = OVERLAP
%   Freq band = [FBL, FBH]
% DATA_N;
% DATA_F_ABS;
% FREQS;
% SPINDLE_POINTS;
%
%@OUTPUT
% SPINDLE_DURFREQAMP;

function spindle_durfreqamp = Spindle_STFT_Classification(parameters, data_n, data_f_abs, freqs, spindle_points)
    parameters = strsplit(parameters);
    if length(parameters) ~= 4
        error('# parameters incorrect. Input format: WINLEN OVERLAP FREQBANDLO FREQBANDHI');
    else
        winlen = str2double(parameters{1});
        overlap = str2double(parameters{2});
        fbl = str2double(parameters{3});
        fbh = str2double(parameters{4});
    end

    nonoverlap = winlen - overlap;
    %% Dur, Freq classification
    data_n_stitch = cell2mat(data_n'); % Cell-->Mat
    spindle_durfreqamp = zeros(size(spindle_points, 1), 3); % # of Cluster Properties
    data_offsets = zeros(1, length(data_f_abs));
    
    eff_winlen = winlen - overlap;
    num_amend = overlap / eff_winlen;
    for i=1:length(data_f_abs)
        data_amend = cell(1, num_amend);
        for j=1:num_amend
            data_amend{j} = zeros(length(data_f_abs{i}{1}), 1);
        end
        data_f_abs_amend{i} = cat(2, data_f_abs{i}, data_amend); % Make length good
        data_offsets(i) = sum(length([data_f_abs_amend{1:i-1}])) * nonoverlap; % Calculate offsets of NREM segments
    end
    
    %Change data struct
    data_stretch = [];
    data_mat = [];
    for i=1:length(data_f_abs_amend)
        data_mat = cat(2, data_mat, cell2mat(data_f_abs_amend{i}));
%         for j=1:length(data_f_abs_amend{i})
%             if(j ~= length(data_f_abs_amend{i}))
%                 data_stretch = cat(2, data_stretch, repmat(data_f_abs_amend{i}{j}, [1, eff_winlen]));
%             else
%                 data_stretch = cat(2, data_stretch, repmat(data_f_abs_amend{i}{j}, [1, overlap]));
%             end
%         end
    end
    
    for i=1:size(spindle_points, 1)
        spindle_durfreqamp(i, 1) = round(spindle_points(i, 2)) - round(spindle_points(i, 1)) + 1;
        spindle_durfreqamp(i, 3) = rms(data_n_stitch(round(spindle_points(i, 1)):round(spindle_points(i, 2))));
        
        sum_freqs_mean = 0;
%         for j=spindle_points(i, 1):spindle_points(i, 2)
%             sum_freqs_mean = sum_freqs_mean + (sum(data_stretch(1:end, j) .* freqs) / sum(freqs)); % Calculation of effective mean frequency
%         end
        start_ind = round((round((spindle_points(i, 1)) - 1) / eff_winlen) + 1);
        stop_ind = start_ind + round(((spindle_durfreqamp(i, 1) - winlen) / eff_winlen));
        
        for j=start_ind:stop_ind
            sum_freqs_mean = sum_freqs_mean + (sum(data_mat(1:end, j) .* freqs) / sum(freqs));
        end
        
        spindle_durfreqamp(i, 2) = sum_freqs_mean / (stop_ind - start_ind + 1);
        
%         spindle_durfreqamp(i, 2) =  sum_freqs_mean / (spindle_durfreqamp(i, 1) / eff_winlen);
        
%         temp_index = find(data_offsets <= spindle_points(i, 1), 1, 'last');
%         temp_start = (spindle_points(i, 1) - 1 - data_offsets(temp_index)) / nonoverlap;
%         temp_end = (spindle_points(i, 2) - data_offsets(temp_index)) / nonoverlap - (winlen / nonoverlap);
%         temp_freqs = zeros(round(length(data_f_abs{temp_index}{1})), round(temp_end - temp_start + 1));
%         temp_offset = 0;
%         temp_flag = false;
%         for j=temp_start:temp_end
%             if ((j + 1) > length(data_f_abs_amend{temp_index})) && (~temp_flag)
%                 temp_flag = true;
%                 temp_offset = length(data_f_abs_amend{temp_index});
%                 temp_index = temp_index + 1;
%             end
%             temp_freqs(:, j - temp_start + 1) = data_f_abs_amend{temp_index}{j + 1 - temp_offset};
%         end
%         temp_freqs_sum_sum = zeros(1, fbh);
%         for j=fbl:fbh
%             temp_freq_low = find(freqs >= j, 1, 'first'); % Lower freq bound
%             temp_freq_high = find(freqs <= j + 1, 1, 'last'); % Upper freq bound
%             temp_freqs_sum_sum(j) = sum(temp_freqs(temp_freq_low:temp_freq_high, :));
%         end
%         spindle_durfreqamp(i, 2) = find(temp_freqs_sum_sum == max(temp_freqs_sum_sum));
%         sum_temp_freqs = sum(temp_freqs, 2);
%         spindle_durfreqamp(i, 2) = freqs(find(sum_temp_freqs == max(sum_temp_freqs), 1));
%         spindle_durfreqamp(i, 4) = sum(temp_freqs_sum_sum);
    end
end