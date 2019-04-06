function Spindle_Analytics_STFT_Revision_Modularized(fname_sp, fpath_sp, fpath_s, fname_edf, b1, b2, wl, wo, ut1, ut2, maxfn)

temp = load(strcat(fpath_sp, fname_sp));
standard_labels = temp.spindle_points;
fname_sp = fname_edf;

identifier=input('Identifier (text):');
analytics_length = input('Analytics Length (frames):');
starttime = input('Start Time (frames):');

% fpath_sp = strcat(fpath_s, '\', num2str(wl), '-', num2str(wo), '_', num2str(b1), '-', num2str(b2));
disp(strcat('WL:', num2str(wl), '; WO:', num2str(wo), '; B1:', num2str(b1), '; B2:', num2str(b2)))
overlap = cell(100, 100);
tpt = zeros(100, 100);
fpt = zeros(100, 100);
fnt = zeros(100, 100);
froc = cell(100, 100);
f1t = zeros(100, 100);
grph = zeros(100, 100);
fprt = zeros(100, 100);
fnrt = zeros(100, 100);
alsize = zeros(100, 100);

for ut=ut1:0.01:ut2
    for lt=ut1:0.01:ut
        temp = load(strcat(fpath_s, '\', fname_edf(1:end-4), 'STFT_', num2str(lt), '-', num2str(ut), '_labels', '.mat'));
        algorithm_labels = temp.spindle_points;
        
%         analytics_length = ceil(max(max(standard_labels)));

        stats = Spindle_Analytics_Core_2(standard_labels, algorithm_labels, analytics_length, starttime);
        
        rlt = round(lt * 100);
        rut = round(ut * 100);
        
        fpt(rlt, rut) = stats.nfp;
        fnt(rlt, rut) = stats.nfn;
        fprt(rlt, rut) = stats.fpr;
        fnrt(rlt, rut) = stats.fnr;
        alsize(rlt, rut) = stats.algotot;
        
        tpt(rlt, rut) = stats.ntp;
        overlap{rlt, rut} = stats.overlap;
        froc{rlt, rut} = [stats.recall, stats.precision];
        f1t(rlt, rut) = stats.f1;
        
        grph(rlt, rut) = stats.f1 * stats.tpr;
    end
end

u = find(max(f1t) == max(max(f1t)));
[maxval, v] = max(f1t(:, u));
for i=1:length(u)
disp(strcat('LT:', num2str(v), '; UT:', num2str(u), '; F1:', num2str(maxval)))
fname_tp = strcat('\', fname_sp(1:end-4), 'STFT_', num2str(v(i)/100), '-', num2str(u(i)/100), '_labels', '.mat');
temp = load(strcat(fpath_s, fname_tp));
algorithm_labels = temp.spindle_points;
algo_size = size(algorithm_labels, 1);
temp_cell{1} = algorithm_labels;
algorithm_labels = d2s3(temp_cell, analytics_length, starttime, 1, 1); % STFT TESTING
algo_size_comb = size(algorithm_labels, 1);

ntp = tpt(v(i), u(i));
nfp = fpt(v(i), u(i));
nfn = fnt(v(i), u(i));

tpr = ntp / algo_size_comb;
fpr = nfp / algo_size_comb;
fnr = nfn / stats.stndtot;

disp(strcat('algo-size:', num2str(algo_size), '; comb-size:', num2str(algo_size_comb)))

stats_txt = {strcat("TP: ", num2str(ntp), ", rate = ", num2str(tpr));
        strcat("FP: ", num2str(nfp), ", rate = ", num2str(fpr));
        strcat("FN: ", num2str(nfn), ", rate = ", num2str(fnr));
        strcat("Recall: ", num2str(ntp / stats.stndtot));
        strcat("Precision: ", num2str(ntp / algo_size_comb));
        strcat("F1 Score: ", num2str(f1t(v(i), u(i))))};

figure;
ax1 = axes('Position',[0 0 1 1],'Visible','off');
ax2 = axes('Position',[.3 .1 .6 .8]);
histogram(overlap{v(i), u(i)}(:, 3), 0:10:100);
axes(ax1);
text(.025, 0.6, stats_txt);
savefig(strcat(fpath_sp, '\', fname_sp(1:end-4), identifier, 'STFT_', num2str(v(i)/100), '-', num2str(u(i)/100), '_Histogram', '.fig'));
end
close;

figure;
hold on;
for a=1:size(froc, 1)
    for b=1:size(froc, 2)
        if(~isempty(froc{a, b}))
            plot(froc{a, b}(1), froc{a, b}(2), 'k.');
        end
    end
end
xlabel('Recall'); ylabel('Precision');
savefig(strcat(fpath_s, '\', fname_sp(1:end-4), identifier, 'STFT_0.1-0.01-0.8_froc', '.fig'));
close;

save(strcat(fpath_s, '\', fname_sp(1:end-4), identifier, 'STFT_0.1-0.01-0.8_labels_stats', '.mat'), 'f1t', 'overlap', 'froc', 'fnt', 'fpt', 'fprt', 'fnrt', 'alsize');

count = 0;
fnrtag = [];
fnrta = [];

for i=1:100
    for j=1:100
        if(fnrt(i,j) == 0)
            fnrta(i, j) = NaN;
        else
            fnrta(i, j) = fnrt(i,j);
        end
        if(fnrta(i,j) < maxfn)
            count = count + 1;
            fnrtag(count, 1) = fnrta(i, j);
            fnrtag(count, 2) = i;
            fnrtag(count, 3) = j;
            fnrtag(count, 4) = fprt(i, j);
            fnrtag(count, 5) = alsize(i, j);
        end
    end
end
if isempty(fnrtag)
    disp('No Parameters fit Selection Criteria.');
else
[~, index] = min(fnrtag(:, 4));
disp(strcat('Selected Parameters:: LT:', num2str(fnrtag(index, 2)), '; UT:', num2str(fnrtag(index, 3))));
end
disp('===================================================');