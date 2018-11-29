function P = mstep_obs(f, moments, nA2D, P, V)

si_unique = unique(reshape(V.settingsindex, 1, []));

for sind = 1:numel(si_unique)
    ii = reshape(find(V.settingsindex == si_unique(sind)), 1, []);
    [num, denom] = deal(0);
    for u = ii
        EFhat    = moments(u).gp_mean + P(u).fdc; %expectation of particles of predicted fluorescence, no noise model taken into account yet
        EFhat_sq = EFhat .^ 2 + moments(u).gp_sd .^ 2;
        num = num + sum(f{u} .^ 2 + EFhat_sq - 2 * f{u} .* EFhat) * nA2D(u);
        denom = denom + numel(f{u});
    end
    zeta = num / denom;
    for u = ii
        P(u).zeta = zeta;
    end
end
