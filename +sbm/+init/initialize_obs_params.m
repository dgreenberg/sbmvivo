function P = initialize_obs_params(P, V, F, dt, nA2D)
%P = initialize_obs_params(P, V, F, dt, nA2D)
%initialize observation parameters. requires parameters and options for calcium and dye-binding dynamics to be set already
init_window = 120; %sec
init_percentile = 0.1;
if ~exist('nA2D','var')
    
    nA2D = ones(1, numel(P));
    
end
nseg = numel(F);
halfnormal_median = 2 * erfinv(0.5); %median of the absolute difference of two independent normal variables with zero mean and unit variance, see https://en.wikipedia.org/wiki/Half-normal_distribution

if numel(intersect(fieldnames(P), {'zeta', 'gain'})) < 2 || numel(intersect(fieldnames(P), {'mu_r_init', 'sigma_r_init'})) < 2
    pg = nan(numel(P), 1); %photon flux times gain
    [Fvar, init_F] = deal(nan(nseg, 1));
    
    %estimate fluorescence noise using a median absolute differences, which are robust to both large outliers due to spiking and slow baseline drifts:
    for u = 1:numel(P)
        T = numel(F{u});
        init_Fvals = F{u}(1:min(T, ceil(init_window / dt(u))));
        init_Fvals = sort(init_Fvals);
        percentile_index = ceil(numel(init_Fvals) * init_percentile);
        init_F(u) = init_Fvals(percentile_index);
        
        abs_diff_F_median = median(abs(diff(F{u}))); %permute randomly so we're dealign with a property of the distrubtion. FIXME lookat autocorrelation etc.???
        %abs_diff_F_median = median(abs(diff(F{u}(randperm(numel(F{u})))))); %permute randomly so we're dealign with a property of the distrubtion. FIXME lookat autocorrelation etc.???
        total_F_sd = abs_diff_F_median / halfnormal_median; %a median-based estimate of fluorescence noise s.d.
        
        total_F_sd = max(total_F_sd, std(F{u}) / 20);
        
        Fvar(u) = total_F_sd ^ 2;
        assert(Fvar(u) > 0, 'constant fluroescence signals cannot be analyzed');
        
        if init_F(u) < 0
            
            init_F(u) = total_F_sd;
            
        end
        
    end
    
    brightness = sbm.model.spikeresponse(P(1), V, min(dt) / 10); %saturation transient for 1 AP. use a finer timestep to make sure we get the max
    b01 = [brightness(1); max(brightness)];  % fixme: below we just use the eq value, so why bother with the simulation?
    
    for u = 1:nseg
        
        P(u).mu_r_init = log(init_F(u) / b01(1));
        
    end
    
    for u = 1:nseg
        P(u).sigma_r_init = log(1.05); %P(u).sigma_r * sqrt(5);
        pg(u) = exp(P(u).mu_r_init) * b01(1); %photon flux times gain at initial baseline calcium and equilibrium binding
        P(u).fdc = 0;
    end
    
    %FIXME include measured darknoise appropriately if available
    %FIXME should take into account segment length when we have multiple segments
    %FIXME fluorescence noise should not account for all variance in F. consider Gaussian mixture on C-based estimate based on tau_c, firing rate, baseline drift, etc.???
    
    [~, ~, settingsindex_renumbered] = unique(V.settingsindex);
    n_settings = max(settingsindex_renumbered);
    
    if true && V.darknoise && ~V.shotnoise
        for u = 1:numel(P)
            
            P(u).gain = 0;
            P(u).zeta = nA2D(u) * Fvar(u) * V.zetafac;
            
        end
    else
    for si = 1:n_settings
        ii = find(settingsindex_renumbered(:) == si);
        nA2D_ii = reshape(nA2D(ii), [], 1);
        if V.shotnoise && V.darknoise
            
            minzeta = min(Fvar(ii) .* nA2D_ii) / 10;
            if numel(ii) > 1
                gz = ([pg(ii) ones(size(ii))] ./ (nA2D_ii * [1 1])) \ Fvar(ii);
                if any(gz < 0) %no negative noise
                    gz = [(pg(ii) ./ nA2D_ii) \ (Fvar(ii) - minzeta ./ nA2D_ii) minzeta];
                end
            else
                gz = [(Fvar(ii) - minzeta ./ nA2D_ii) / (pg(ii) / nA2D_ii) minzeta];
            end
        elseif V.shotnoise
            gz = [(pg(ii) ./ nA2D_ii) \ Fvar(ii) 0];
        elseif V.darknoise
            gz = [0 (1 ./ nA2D_ii) \ Fvar(ii)];
        end
        for u = reshape(ii, 1, [])
            [P(u).gain, P(u).zeta] = deal(gz(1), gz(2));
        end
    end
    end
    
end
P = orderfields(P);