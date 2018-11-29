gfiltwidth = 0.2; %sec

nseg = numel(data);

neuron_index = cellfun(@(v) v.cell_num, data);
nneurons = max(neuron_index);

gfiltwidth_bin = 0.2 / (1 / data{1}.fps); %FIXME check whether FPS is constant

[nstrue_filtered, nsest_filtered] = deal(cell(1, nneurons));

for k = 1:nseg
    
    j = data{k}.cell_num;
    
    
    nstrue_filtered{j} = [nstrue_filtered{j}; gauss_filter(data{k}.outputs', gfiltwidth_bin)];
    
    n = gauss_filter(data{k}.predictions', gfiltwidth_bin);
    
    extrabins = numel(data{k}.predictions) - numel(data{k}.outputs);
    pad_left = ceil(extrabins / 2);
    pad_right = extrabins - pad_left;
        
    nsest_filtered{j}  = [nsest_filtered{j};  n(1 + pad_left:end - pad_right)];
    
    tt = (1:numel(nstrue_filitered{j}));
    
end

c = nan(1, nneurons);
for j = 1:nneurons
    if isempty(nstrue_filtered{j}), continue; end
    cmat = corrcoef(nstrue_filtered{j}, nsest_filtered{j});
    c(j) = cmat(2, 1);
end