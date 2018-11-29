function [r, expr, log_w, log_w_corrected, gp, vFtotal] = sample_initial_r(P, V, nA2D, F, Fadj, brightness, brightness_eq, nrand_r, randtype, log_pq_spiking, vF, randseed)
%sample r conditioned on s and Fadj(1)
sigma_r_init = log(4);

F_init = max(abs(F(1) * 0.05), Fadj(1) * 0.95); % always nonnegative

if V.shotnoise

    vFtotal_init = vF + (P.gain / nA2D) * F_init;  % always positive
    
else
    
    vFtotal_init = vF;
    
end

assert(F_init > 0, 'invalid fluorescence data');
mu_r_init = log(F_init / brightness_eq);

%fixme a lot of these factors cancel?! also see if cuda code can be
%simplified

rmaxexp = F_init ./ brightness;  % always positive
rmax = log(rmaxexp);  % max likelihood r given laplace approximation at fcorr_approx
dfac = rmaxexp .* brightness; %denominator factor for calculating a Gaussian approximation to P[F | r]
likvarr = vFtotal_init ./ dfac .^ 2; %variance of a Gaussian approximation of likelihood P[F | r]
likmeanr = F_init ./ dfac + rmax - 1.0; %mean of a Gaussian approximation of likelihood P[F | r]

%for each particle combine P[F | r, brightness] and P[r] to get P[r | F, brightness] using the formula for a product of Gaussians
v_r = 1 ./ (1 ./ likvarr + 1 / sigma_r_init ^ 2);
m_r = v_r .* (likmeanr ./ likvarr + mu_r_init / sigma_r_init ^ 2);
if V.singleprecision
    
    randseed = quickrandn(nrand_r, randseed); %#ok<NASGU>
    
else
    
    nrand_r = randn(V.Nparticles, 1, randtype); %normal variates for state r
    
end
r = m_r + sqrt(v_r) .* nrand_r;

%calculate probability weights
log_priorp_r    = -0.5 * ((r - mu_r_init) .^ 2  / sigma_r_init ^ 2 + log(2 * pi * sigma_r_init ^ 2));
log_samplingp_r = -0.5 * ((r -       m_r) .^ 2 ./ v_r              + log(2 * pi * v_r));
expr = exp(r);
gp = expr .* brightness;
if V.shotnoise
    
    vFtotal = vF + (P.gain / nA2D) * gp;
    log_w = -0.5 * ((gp - Fadj(1)) .^2 ./ vFtotal + log(vFtotal * 2 * pi)) + log_priorp_r - log_samplingp_r;
    
else
    
    vFtotal = nan(V.Nparticles, 1, randtype);
    log_w = -0.5 * ((gp - Fadj(1)) .^2  / vF      + log(vF      * 2 * pi)) + log_priorp_r - log_samplingp_r;
    
end
log_w_corrected = log_w + log_pq_spiking;