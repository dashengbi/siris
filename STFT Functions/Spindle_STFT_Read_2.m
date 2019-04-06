%@INPUT
% FPATH-->Path of files;
% FNAME_EDF-->Name of .edf file;
% FNAME_S-->Name of .txt scoring data;
%
%@OUTPUT
% HDR;
% RECORD;
% TMP_SCORE;


function [hdr, record, tmp_score] = Spindle_STFT_Read_2(fpath_edf, fname_edf, fpath_s, fname_s)
    [hdr, record] = edfread([fpath_edf fname_edf]);
    tmp_fname = [fpath_s fname_s];
    fid = fopen(tmp_fname);
    fgets(fid); 
    tline = fgets(fid); 
    count = 0;
    clear tmp_score;
    while tline ~= -1
        tmp = strsplit(tline,',');
        count = count + 1;
        tmp_score(count) = str2num(tmp{end-1});
        tline = fgets(fid);
    end
    fclose(fid);
end