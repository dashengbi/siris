function Spindle_STFT_Revision_Modularized(ch, epoch, fpath_edf, fname_edf, fpath_s, fname_s, b1, b2, wl, wo, ut1, ut2)

global f fs nfft win;

[hdr, record, tmp_score] = Spindle_STFT_Read_2(fpath_edf, fname_edf, fpath_s, fname_s);
[data_n, ind_NREM] = Spindle_STFT_Preprocessing(ch, epoch, hdr, record, tmp_score);

% [data_n, ind_NREM] = Spindle_STFT_Cut(data_n, ind_NREM, time * fs);

        mkdir(strcat('Spindle_Labels_STFT\STFT_Revision_', fname_edf(1:end-4), '\', num2str(wl), '-', num2str(wo), '_', num2str(b1), '-', num2str(b2)));
        [data_f_abs, freqs] = Spindle_STFT_Core(strcat(num2str(wl), " ", num2str(wo)), data_n, ind_NREM);

for ut=ut1:0.01:ut2
    for lt=ut1:0.01:ut %Changeable
        tic;
        spindle_points = Spindle_STFT_Find(strcat(num2str(wl), " ", num2str(wo), " ", num2str(b1), " ", num2str(b2), " ", num2str(lt), " ", num2str(ut), " ", '0'), data_n, data_f_abs, freqs, ind_NREM);
        spindle_durfreqamp = Spindle_STFT_Classification(strcat(num2str(wl), " ", num2str(wo), " ", num2str(b1), " ", num2str(b2)), data_n, data_f_abs, freqs, spindle_points);

        save(strcat('Spindle_Labels_STFT\STFT_Revision_', fname_edf(1:end-4), '\', num2str(wl), '-', num2str(wo), '_', num2str(b1), '-', num2str(b2), '\', fname_edf(1:end-4), 'STFT_', num2str(lt), '-', num2str(ut), '_labels', '.mat'), 'spindle_points', 'spindle_durfreqamp');
        toc;
    end
end
disp('End STFT Processing.');