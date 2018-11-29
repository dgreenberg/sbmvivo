function [states, states_eq] = fspikecountsfit_5stateb_odesolver_multi_kdex( ...
    stepsize, spikecounts, P, states, subseglist, usemex)

odemethod = 'backwardeuler';

if ~exist('usemex', 'var'), usemex = true; end
if ~exist('subseglist', 'var'), subseglist = 1:numel(spikecounts); end

niter_newton = 3;  % for backward euler

nseg = numel(spikecounts);
nbuffers = numel(P(1).Btot);

assert(all(arrayfun(@(v) numel(v.Btot), P) == nbuffers), 'all data must have the same number of buffers');
tau_B = cat(1, P.tau_B);
kd_B = cat(1, P.kd_B);
c0 = cat(1, P.c0);
Btot = cat(1, P.Btot);
A = cat(1, P.A);
S = cat(1, P.S);
g = 1 ./ cat(1, P.tau_ex);
kon = cat(1, P.kon);
koff = cat(1, P.koff);
nbindingsteps = size(kon, 2);
koff_B = 1 ./ tau_B;
kon_B = koff_B ./ kd_B;
kd_ex = cat(1, P.kd_ex);
maxex = g .* kd_ex;

%get equilbirium states
beta = [ones(nseg, 1) cumprod(kon ./ koff, 2)];
c0powers = bsxfun(@power, c0, 0:nbindingsteps);
fracs_eq = beta .* c0powers;
fracs_eq = bsxfun(@rdivide, fracs_eq, sum(fracs_eq, 2));

if nbuffers > 0
    
    B0 = bsxfun(@times, Btot, c0) ./ bsxfun(@plus, c0, kd_B);  % concentration of binding sites bound to calcium for each buffer (this is simply the Hill equation with exponent 1)
    
else
    
    B0 = [];
    [Btot, kon_B, koff_B] = deal(zeros(numel(states), 0));
    
end

states_eq = [c0, bsxfun(@times, S, fracs_eq), B0]';

if strcmpi(odemethod, 'backwardeuler') && usemex
    
    fspikecountsfit_5state_odesolver_be_mainloop_kdex( ...
        spikecounts, A, states_eq, kon', koff', kon_B', koff_B', Btot', maxex', c0, kd_ex', niter_newton, states, subseglist, stepsize);
    return;
    
end
    
for ii = 1:numel(subseglist)
    
    jj = subseglist(ii);
    
    dfunc = @(t, y) ratefunc_5s(y, koff(jj, :)', kon(jj, :)', koff_B(jj, :)', kon_B(jj, :)', maxex(jj, :), c0(jj), Btot(jj, :)', kd_ex(jj, :));
    
    odeopts = odeset('reltol', 1e-3);
    
    Jfunc = @(t, y) ratefuncjac_5s(y, koff(jj, :)', kon(jj, :)', koff_B(jj, :)', kon_B(jj, :)', maxex(jj, :), Btot(jj, :)', kd_ex(jj, :));
    
    if strcmpi(odemethod, '15s')
        
        odeopts = odeset(odeopts, 'NonNegative', (1:size(states_eq, 1))', 'Jacobian', Jfunc);
        odefunc = @(tlist, y0) ode15s(dfunc, tlist, y0, odeopts);
        
    elseif strcmpi(odemethod, '23s')
        
        odeopts = odeset(odeopts, 'Jacobian', Jfunc);  % non-negativity constraint not possible
        odefunc = @(tlist, y0) ode23s(dfunc, tlist, y0, odeopts);
        
    elseif strcmpi(odemethod, '23t')
        
        odeopts = odeset(odeopts, 'NonNegative', (1:size(states_eq, 1))', 'Jacobian', Jfunc);
        odefunc = @(tlist, y0) ode23t(dfunc, tlist, y0, odeopts);
        
    elseif strcmpi(odemethod, '23tb')
        
        odeopts = odeset(odeopts, 'NonNegative', (1:size(states_eq, 1))', 'Jacobian', Jfunc);
        odefunc = @(tlist, y0) ode23tb(dfunc, tlist, y0, odeopts);
        
        
    elseif strcmpi(odemethod, 'backwardeuler')
        
        odefunc = @(tlist, y0) odesolve_backward_euler(tlist, y0, dfunc, Jfunc, niter_newton);
        
    else
        
        error('unrecognized odemethod')
        
    end
    
    states{jj} = usesolver_5s(odefunc, spikecounts{jj}, A(jj), stepsize(jj), states_eq(:, jj));    
    
end


function states = usesolver_5s(odefunc, spikecounts, A, stepsize, yprev)

ii_spikes = find(spikecounts > 0);
if isempty(ii_spikes)
    
    states = yprev * ones(1, numel(spikecounts));
    return;
    
end

tlist = (1:numel(spikecounts)) * stepsize;
states = nan(numel(yprev), numel(spikecounts));
%copy equilbirium state until the first spike:
states(:, 1:ii_spikes(1) - 1) = yprev * ones(1, ii_spikes(1) - 1);

if ii_spikes(1) == 1
    
    tprev = 0;
    
else
    
    tprev = tlist(ii_spikes(1) - 1);
    
end

for j = 1:numel(ii_spikes)
    
    yprev(1) = yprev(1) + A * spikecounts(ii_spikes(j));
    
    if j == numel(ii_spikes)
        
        ii_next = numel(spikecounts);
        
    else
        
        ii_next = ii_spikes(j + 1) - 1;
        
    end
    
    ii = ii_spikes(j):ii_next;
    [~, nexty] = odefunc([tprev, tlist(ii)], yprev);
    
    if numel(ii) == 1  % solver will return all points in this case
        
        states(:, ii) = nexty(end, :)';
        
    else
        
        states(:, ii) = nexty(2:end, :)';
        
    end
    
    tprev = tlist(ii_next);
    yprev = states(:, ii(end));
    
end


function [t, y] = odesolve_backward_euler(t, y0, ratefunc, jacfunc, niter)
n = numel(y0);
y = nan(n, numel(t));
y(:, 1) = y0;

for j = 2:numel(t)
    
    h = t(j) - t(j - 1);
    ynext = y(:, j - 1);
    yprev = y(:, j - 1);
    
    for k = 1:niter % Newton iterations
        
        f = ratefunc(t + h, ynext);
        J = jacfunc(t + h, ynext);
        z = yprev - ynext + h * f;
        dy = (eye(n) - h * J) \ z;
        ynext = ynext + dy;
        
    end
    
    y(:, j) = ynext;
    
end

y = y';


function dy_dt = ratefunc_5s(y, koff, kon, koff_B, kon_B, maxex, c0, Btot, kdex)
nbindingsteps = numel(kon);
c = y(1); s = y(2:2 + nbindingsteps); b = y(2 + nbindingsteps + 1:end);

r_ex = sum(maxex .* (c ./ (c + kdex) - c0 ./ (c0 + kdex)));

r = kon .* s(1:end - 1) * c - koff .* s(2:end);
ds = [0; r] - [r; 0];

db = (Btot - b) .* kon_B * c - b .* koff_B;  % derivative of bound state only
dc = -sum(r) - sum(db) - r_ex;  % free calcium
dy_dt = [dc; ds; db];


function J = ratefuncjac_5s(y, koff, kon, koff_B, kon_B, maxex, Btot, kdex)
nbindingsteps = numel(kon);
c = y(1); s = y(2:2 + nbindingsteps); b = y(2 + nbindingsteps + 1:end);
ns = numel(s);
nb = numel(b);

d2c_dtdc = -sum(maxex .* kdex ./ (c + kdex) .^ 2) - kon' * s(1:ns - 1) - kon_B' * (Btot - b);
d2c_dtds = [0,   koff'] - [kon' * c,   0];
d2c_dtdb = koff_B' + kon_B' * c;  % only bound states are included in y

d2s_dtdc = [0; s(1:ns - 1) .* kon] - [s(1:ns - 1) .* kon; 0];
d2s_dtds = zeros(ns);
d2s_dtds(2:ns + 1:end) = c * kon;  % d2s(i)_dtds(i-1)
d2s_dtds(1:ns +1:end) = [0; -koff] - [c * kon; 0];  % d2s(i)_dtds(i)
d2s_dtds(ns + 1:ns + 1:end) = koff;  % d2s(i)_dtds(i + 1)
d2s_dtdb = zeros(ns, nb);

d2b_dtdc = kon_B .* (Btot - b);
d2b_dtds = zeros(nb, ns);
d2b_dtdb = diag(-d2c_dtdb);  % only bound states are included in y

J = [d2c_dtdc d2c_dtds d2c_dtdb; ...
    d2s_dtdc d2s_dtds d2s_dtdb; ...
    d2b_dtdc d2b_dtds d2b_dtdb];