function [nstrue_filtered, nsest_filtered, it_masked, csegment, cneuron] = oedb_masked_filter(nstrue_allsegs, nsest_allsegs, tol, it, clipwin) %calculate filtered correlation for one neuron
if ~exist('clipwin','var')
    clipwin = clipwin_corr;
end
nseg = numel(nstrue_allsegs);
[nstrue_filtered, nsest_filtered, it_masked] = deal(cell(1, nseg));
csegment = nan(1, nseg);
dt = cellfun(@(v) median(diff(v)), it);

for w = 1:nseg
    
    nstrue = reshape(nstrue_allsegs{w},[],1); nsest = reshape(nsest_allsegs{w},[],1);
    [nstrue_filtered{w}, nsest_filtered{w}] = deal(nan(size(nsest)));
    
    if isempty(nsest) || isempty(nstrue), continue; end
    
    filtwidth_bins = tol / dt(w);
    
    nsmask = ~isnan(nstrue) & ~isnan(nsest);
    nsmask(1:min(end, ceil(clipwin(1) / dt(w)))) = false;
    nsmask(max(1, end - ceil(clipwin(2) / dt(w)) + 1):end) = false;
    it_masked{w} = it{w}(nsmask);
    s = find(nsmask & ~[false; nsmask(1:end-1)]);
    e = find(nsmask & ~[nsmask(2:end); false]);
    for kk = 1:numel(s)
        ii = s(kk):e(kk);
        nstrue_filtered{w}(ii) = gauss_filter(nstrue(ii), filtwidth_bins);
        nsest_filtered{w}(ii)  = gauss_filter(nsest(ii), filtwidth_bins);
    end
    nstrue_filtered{w} = nstrue_filtered{w}(nsmask);
    nsest_filtered{w} = nsest_filtered{w}(nsmask);
    cmat = corrcoef(nstrue_filtered{w}, nsest_filtered{w});
    csegment(w) = cmat(2);
    
end
cmat = corrcoef_nonnan(cat(1, nstrue_filtered{:}), cat(1, nsest_filtered{:}));
cneuron = cmat(2);