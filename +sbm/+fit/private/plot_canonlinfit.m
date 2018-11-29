function figh = plot_canonlinfit(states, states_eq, t_states, it, fBL, f, fhat, figname, blmask, rmse, nsubseg, fdc, R, S, ...
    statenames, perispikewin, basefigsavedir, states_scalefac, mindt)
if ~isempty(basefigsavedir) && basefigsavedir(end) ~= filesep
    
    basefigsavedir(end + 1) = filesep;
    
end    
if ~exist('mindt', 'var')
    
    mindt = 1e-3;
    
end
if ~exist('states_scalefac', 'var')
    
    states_scalefac = 1;
    
end
tdiff = diff(t_states);
tdiff(isnan(tdiff)) = [];
dt_states = median(tdiff);
decimation_factor = 1;

if dt_states < mindt
    
    decimation_factor = ceil(mindt / dt_states);
    
end

nstates = size(states, 1);
nax = nstates + 1;
figh = figure;
if ~isempty(figname)
    
    set(figh, 'name', figname);
    
end

ax(1) = subplot(nax,1,1);

plot(it, fBL, 'k');
hold on
plot(it,f, 'r', 'linewidth',2)
fmasked = f;
fmasked(~blmask) = nan;
plot(it, fmasked, 'm');
plot(it, fhat, '.-')

title(sprintf('rmse = %.8g, %d subsegments, fdc = %.5g, S = %.5g \\muM', rmse, nsubseg, fdc, S * states_scalefac / 1000));
ylabel(sprintf('R = %.5g', R));

XL = [min(it) - perispikewin(2), max(it)]; %FIXME start 1 modeling timestep earlier

for k = 1:nstates
    
    ax(k + 1) = subplot(nax,1, k + 1);
    plot(t_states(1:decimation_factor:end), states_scalefac * states(k, 1:decimation_factor:end)', '.-');
    hold on;
    plot(XL', states_scalefac * states_eq(k) * [1; 1], ':');
    ylabel(statenames{k});
    
end

set(ax, 'xlim', XL)
%linkaxes(ax, 'x');

if ~isempty(figname) && ~isempty(basefigsavedir)
    
    figname = get(figh,'name');
    figname(ismember(figname, '\/:')) = ' ';
    saveas(figh,[basefigsavedir figname]);
    
end
if nargout == 0

    close(figh);
    
end