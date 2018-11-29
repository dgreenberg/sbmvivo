function meanfvstdiff(oedb, fileinfo, dsind, neuronind, w)
if ~exist('w','var')
    
    w = 1;
    
end

maxXL_scatter = 0.3;

%data
oerec = fetch_neurondata(oedb, fileinfo, dsind, neuronind);

[~, ~, ~, ~, responsedata, ~, ~, ~, ~, ~, ~, stdiff] = avspikeresponse(oerec);

rsingles = responsedata{1};
asingles = cellfun(@(v) mean(v(v(:,1) >= 0 & v(:,1) <= w, 2)), rsingles, 'uniformoutput', false); asingles = cat(2, asingles{:}); asingles = asingles(:);
%asingles = cellfun(@(v) max(v(v(:,1) >= 0 & v(:,1) <= w, 2)), rsingles, 'uniformoutput', false); asingles = cat(2, asingles{:}); asingles = asingles(:);
fsingles = find(~isnan(asingles) & cellfun(@(v) max(v(:,1)) >= w, rsingles)');
singlemean = mean(asingles(fsingles));

r = responsedata{2};
a = cellfun(@(v) mean(v(v(:,1) >= 0 & v(:,1) <= w, 2)), r, 'uniformoutput', false); a = cat(2, a{:}); a = a(:);
%a = cellfun(@(v) max(v(v(:,1) >= 0 & v(:,1) <= w, 2)), r, 'uniformoutput', false); a = cat(2, a{:}); a = a(:);
d = cat(1,stdiff{2}{:});

f = find(~isnan(a) & cellfun(@(v) max(v(:,1)) >= w, r)');
[~, si] = sort(d(f),'descend');
fs = f(si);


%model
algind = find(strcmpi(oedb.algnames, 'sbm'));
if numel(algind) ~= 1, return; end
results = fetch_results(oedb, fileinfo, dsind, neuronind, 1, algind);
P = results.params;
V = results.opts;
domodel = ~isempty(P) && isa(P, 'struct') && ~isempty(fieldnames(P));
if domodel
    [~, ~, dt] = extract_data_fromoedatasegment(oerec.data(1));
    
    twin = w; %sec
    nsteps = ceil(dt / V.substepsize); %10 ms increments or so
    stepsize = dt / nsteps;
    nsteps_sim = ceil(twin / stepsize) + 1;
    
    maxisi_insteps = ceil(maxXL_scatter / stepsize) + 1;
    
    for k = 1:maxisi_insteps
        
        nsvec = zeros(nsteps_sim, 1);
        nsvec(2) = 1;
        nsvec(2 + k - 1) = nsvec(2 + k - 1) + 1;
        
        [brightness(k, :), ~, ~, tt_brightness] = sbm.model.spiketrainresponse(P,V,dt,nsvec);
        
    end
    
    brightness0 = brightness(1, 1);
    DFF = (brightness - brightness0)' / brightness0;    
    
    ii = tt_brightness >= stepsize & tt_brightness <= stepsize + w;    
    asim = mean(DFF(ii,:), 1)';
    dsim = stepsize * ((1:maxisi_insteps)' - 1);
end

%plots
figh = figure;

spacer = 1 / 3;
subplot(1,4,1);
for k = 1:numel(fs)
    
    j = fs(k);
    nextr = r{j};
    nextd = d(j);
    plot(nextr(:,1), nextr(:,2) + spacer * k); hold on
    plot(nextr([1 end],1), [1;1] * spacer * k, 'k')
    plot([0; 0], spacer * k + [-1; 1] * spacer / 4,'k'); %first spike
    plot([0; 0] + nextd, spacer * k + [-1; 1] * spacer / 4, 'k'); %second spike
    
end
xlim([-0.5 max(w, 1)])

subplot(1,4,2);
plot(d(fs) * 1000, (1:numel(fs))' * spacer);
set(gca,'ytick',  (1:numel(fs))' * spacer, 'yticklabel', num2str((1:numel(fs))'));
xlabel('Interspike interval [ms]')

subplot(1,4,3);
plot(a(fs) * 100, (1:numel(fs))' * spacer);
hold on;
plot([1; 1] * singlemean * 100, ylim', 'm--')
set(gca,'ytick',  (1:numel(fs))' * spacer, 'yticklabel', num2str((1:numel(fs))'));
xlabel('Mean \Delta F/F_0 [%]')

subplot(1,4,4);
plot(d(fs), a(fs), '.', 'markersize', 20)
q = [d(fs) ones(numel(fs), 1)] \ a(fs);
[corr, pval] = corrcoef(a(fs), d(fs));
set(gca,'xtick',[0:.025:maxXL_scatter],'xticklabel',num2str(1000 * [0:.025:maxXL_scatter]'))
axis([0 maxXL_scatter 0 0.1 * ceil(max(a(fs) / 0.1))]);
XL = xlim;
hold on; plot(XL', XL' * q(1) + q(2), 'r')
plot(XL', [1; 1] * singlemean, 'm--')
if domodel
    
    plot(dsim, asim, 'g');
    
end

text(min(xlim), min(ylim) + diff(ylim) / 3, {sprintf('\\Delta F/F_0 = %f%% - %f%% / ms', q(2) * 100, q(1) * 100 / 1000); sprintf('r = %f, p = %f', corr(2), pval(2))});
xlabel('Interspike interval [ms]')
ylabel('Mean \Delta F/F_0 [%]')

set(findobj(figh, 'type', 'axes'), 'tickdir', 'out', 'box', 'off');

if domodel
    
    figh = figure;
    plot([-0.5 tt_brightness], DFF([1 1:end], :));
    axis([-.5 1 -.05 .3])
    
    set(findobj(figh, 'type', 'axes'), 'tickdir', 'out', 'box', 'off');
    
end