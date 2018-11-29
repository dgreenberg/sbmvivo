function oedb_firingrateplot(oedb, ds)
frthresh = 1; %Hz
unitless = ismember(oedb.algnames, {'CFOOPSI', 'FOOPSI'});

frest = arrayfun(@(v) nanforempty(v.fr_est_spiketimes), oedb.stats.byneuron{ds})';
frtrue = arrayfun(@(v) nanforempty(v.fr_true_spiketimes), oedb.stats.byneuron{ds})';
algok = find(any(~isnan(frest) & ~isnan(frtrue), 1) & ~unitless);

maxfr_est = max(max(frest(:, algok)));
maxfr_true = max(max(frtrue(:, algok)));
lowfr = frtrue(:, 1) < frthresh;

if ~any(algok)
   
    warning('no valid data for firing rate plot');
    return;
    
end

[cvals,svals] = deal(nan(1, numel(algok)));
for k = 1:numel(algok)
    
    svals(k) = spearman_rankcorr(frest(:, algok(k)), frtrue(:, algok(k)));
    
    cmat = corrcoef(frest(lowfr, algok(k)), frtrue(lowfr, algok(k)));
    cvals(k) = cmat(2);
    
end

[~, si] = sort(svals,'descend');

figh = figure('name',sprintf('true vs. estimated firing rates for %s', oedb.datasetnames{ds}));

ax(1) = subplot(1, 4, 1, 'parent', figh);
h = plot(frtrue(:, algok), frest(:, algok), 'o', 'parent', ax(1));
set(ax(1), 'xlim', [0 frthresh], 'ylim', [0 maxfr_est * 1.1], 'dataaspectratio', [1 1 1]);
line([0 frthresh], [0 frthresh],'parent',ax(1),'color',[1 1 1] * .6);

ax(2) = subplot(1, 4, 2, 'parent', figh);
h = plot(frtrue(:, algok), frest(:, algok), 'o', 'parent', ax(2));
legend(h, oedb.algnames(algok));
set(ax(2), 'xlim', [frthresh maxfr_true * 1.1], 'ylim', [0 maxfr_est * 1.1]);
line([frthresh maxfr_true * 1.1], [frthresh maxfr_true * 1.1],'parent',ax(2),'color',[1 1 1] * .6);

ax(3) = subplot(1, 4, 3, 'parent', figh);
bar(max(0, svals(si)),'parent',ax(3));
set(ax(3), 'xtick', 1:numel(algok), 'xticklabel', oedb.algnames(algok(si)))
for k = 1:numel(algok)
    
    text(k, max(0, svals(si(k))) * 1.025 + 0.01, num2str(svals(si(k))),'parent',ax(3));
    
end
ylabel('spearman rank corr');

ax(4) = subplot(1, 4, 4, 'parent', figh);
bar(max(0, cvals(si)),'parent',ax(4));
set(ax(4), 'xtick', 1:numel(algok), 'xticklabel', oedb.algnames(algok(si)))
for k = 1:numel(algok)
    
    text(k, max(0, cvals(si(k))) * 1.025 + 0.01, num2str(cvals(si(k))),'parent',ax(4));
    
end
ylabel(sprintf('pearson corr, fr < %f Hz', frthresh));

set(ax(2),'plotboxaspectratio', get(ax(1),'plotboxaspectratio'));
set(ax, 'tickdir', 'out', 'box', 'off');

set(findobj(figh,'marker','o'),'markersize',10,'linewidth',1)

assignin('base', 'frest', frest);
assignin('base', 'frtrue', frtrue);
assignin('base', 'cvals', cvals);
assignin('base', 'svals', svals);