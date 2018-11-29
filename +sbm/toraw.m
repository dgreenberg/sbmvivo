function toraw(outputfile,F,dt,V,P,true_st,nA2D,q_spike)
if ~exist('nA2D','var'), nA2D = 1; end
if ~exist('true_st','var'), true_st = []; end
if ~exist('q_spike','var') || isempty(q_spike)
    q_spike = [];
end
lik = empty_likstruct;

%fixme use a standard function to assign params and settings etc.
if ~isfield(V, 'ResampleThreshold')
    V.ResampleThreshold = V.Nparticles / 2;
end

[Fadj, vF, T, randtype, spikesknown, nsteps, stepsize, true_n_sub, p_spike, ...
    nrand_r, quickresample, randseed, guardfac, baseresamp, parents_noresamp, equal_log_w, ...
    q_spike, log_pq_spike, log_pq_nospike, log_pq_spiking, ntimepoints_pre] ...
    = pfilter_precalc(P, V, dt, F, true_st, nA2D, q_spike);

K = ceil(V.filtersmoother_window / dt);
min_prewin = 2.5;
ntimepoints_pre = ceil(max(V.filtersmoother_window, min_prewin) / dt);

[c_pre, x_pre, s_pre, n_pre, r_pre, gp, expr, log_w, log_w_corrected, vFtotal] = sample_initial_states(P, V, F, dt, true_st, nA2D, ...
    q_spike(1:nsteps * (ntimepoints_pre + 1)));

[w_pre, log_w_pre,           ~,     ~] = normalize_weights(log_w);
[~,     log_w_corrected_pre, neff0, log_sum_raw_w0] = normalize_weights(log_w_corrected);

fid = fopen(outputfile,'w');
fwrite(fid, T, 'int32');
fwrite(fid, V.Nparticles, 'int32');
fwrite(fid, nsteps, 'int32');
fwrite(fid, K, 'int32');
fwrite(fid, ntimepoints_pre, 'int32');
fwrite(fid, dt, 'single');
fwrite(fid, V.ResampleThreshold, 'single');
fwrite(fid, neff0, 'single');
fwrite(fid, log_sum_raw_w0, 'single');
fwrite(fid, F, 'single');
fwrite(fid, rand(1, T, 'single'), 'single'); %u
fwrite(fid, q_spike, 'single');
fwrite(fid, log_pq_spike, 'single');
fwrite(fid, log_pq_nospike, 'single');

pvec = zeros(0, 1);

if strmcpi(V.model, '5s2b')
    
    pvec(end + 1:end + 4, 1) = P.dbrightness;
    pvec(end + 1, 1) = P.zeta / nA2D;
    pvec(end + 1, 1) = P.sigma_r;    
    pvec(end + 1, 1) = P.fr * dt; %lambda, the expected spike count per time bin. not to be confused with p_spike = lambda / nsteps
    pvec(end + 1, 1) = P.S;
    pvec(end + 1, 1) = P.R + 1; % FBGp1
    pvec(end + 1, 1) = P.gain / nA2D;    
    pvec(end + 1, 1) = P.kd_ex;    
    pvec(end + 1, 1) = P.kd_ex / P.tau_ex;
    pvec(end + 1, 1) = P.c0;
    pvec(end + 1, 1) = P.fdc;
    pvec(end + 1:end + 2, 1) = P.Btot;
    pvec(end + 1:end + 4, 1) = P.kon;
    pvec(end + 1:end + 4, 1) = P.koff;
    pvec(end + 1:end + 4, 1) = kon_B(1);
    pvec(end + 1:end + 4, 1) = koff_B(1);
    pvec(end + 1:end + 4, 1) = kon_B(2);
    pvec(end + 1:end + 4, 1) = koff_B(2);
    
else
    
    error('unrecognized model: %s', V.model);
    
end

pvec = single(pvec);
fwrite(fid, pvec, 'single');

fwrite(fid,c_pre(:,end),'single');
fwrite(fid,s_pre(:,end),'single');
fwrite(fid,x_pre(:,end),'single');
fwrite(fid,r_pre,'single');
fwrite(fid,w_pre,'single');
fwrite(fid,log_w_pre,'single');
fwrite(fid,log_w_corrected_pre,'single');
fwrite(fid,uint8(n_pre),'uint8');
fwrite(fid,gp,'single');

fclose(fid);
