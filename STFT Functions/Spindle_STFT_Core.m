%@INPUT
% {PARAMETERS}:
%   Window length = WINLEN
%   Overlap = OVERLAP
% DATA_N;
% IND_NREM;
%
%@OUTPUT
% DATA_F_ABS;
% FREQS;

function [data_f_abs, freqs] = Spindle_STFT_Core(parameters, data_n, ind_NREM)
    parameters = strsplit(parameters);
    if length(parameters) ~= 2
        error('# parameters incorrect. Input format: WINLEN OVERLAP');
    else
        winlen = str2double(parameters{1});
        overlap = str2double(parameters{2});
    end

    %% Define parameters for Short-Time Fourier Transform
    window = hamming(winlen); % Window length

    %% Short-Time Fourier Transform
    clear data_f;
    for i=1:size(ind_NREM, 1)
        [data_f{i}, freqs, ~] = spectrogram(data_n{i}, window, overlap, [], 1000); % Check parameters for this function
    end

    %% Amplitude normalization and conversion to real
    clear data_f_abs;
    for i=1:size(ind_NREM, 1)
        for j=1:size(data_f{i}, 2)
            data_f_abs{i}{j} = abs(data_f{i}(:, j)); % NO NORMALIZATION APPLIED
        end
    end

end