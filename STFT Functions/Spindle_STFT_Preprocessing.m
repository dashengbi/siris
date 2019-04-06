%@INPUT
% CH-->Channel;
% EPOCH-->Length of epoch;
% HDR-->Header of .edf file;
% RECORD--> .edf file data;
% TMP_SCORE--> .txt file scoring data;
%
%@OUTPUT
% DATA_N;
% IND_NREM;

function [data_n, ind_NREM] = Spindle_STFT_Preprocessing(ch, epoch, hdr, record, tmp_score)
    global f fs nfft win;    
    %Properties of raw data
    fs = hdr.frequency(1); % 1000 (Hz) = frequency of data sampling
    win = fs*2;
    nfft = 2^nextpow2(win); % 2048
    f = fs/2*linspace(0,1,nfft/2+1);
    sig_band = find(f>10,1):find(f<16,1,'last');

    %Preprocessing
    ind = find(tmp_score == 1);     % Wake->Future use?
    epoch_WAKE = find_seg(ind');
    ind_WAKE(:,1) = (epoch_WAKE(:,1)-1)*epoch*fs+1;
    ind_WAKE(:,2) = (epoch_WAKE(:,2))*epoch*fs;
    ind = find(tmp_score == 2);     % NREM
    epoch_NREM = find_seg(ind');
    ind_NREM(:,1) = (epoch_NREM(:,1)-1)*epoch*fs+1;
    ind_NREM(:,2) = (epoch_NREM(:,2))*epoch*fs;
    ind = find(tmp_score == 3);
    epoch_REM = find_seg(ind');     % REM->Future use?
    ind_REM(:,1) = (epoch_REM(:,1)-1)*epoch*fs+1;
    ind_REM(:,2) = (epoch_REM(:,2))*epoch*fs;

    for i=1:size(ind_NREM,1)
        ind = ind_NREM(i,1):ind_NREM(i,2);
        data = record(ch,ind);
        data_n{i} = pre_process(data); % Smoothing data and eliminating some noise
    end
end