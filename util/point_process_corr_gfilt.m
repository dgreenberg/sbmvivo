function [corr, eab, ea, eb, easq, ebsq, corr_eachseg] = point_process_corr_gfilt(ta, tb, w, twin, ha, hb)
% [corr, eab, ea, eb, easq, ebsq, corr_eachseg] = point_process_corr_gfilt(ta, tb, w, twin, ha, hb)
%
% calculate the correlation between two point processes after Gaussian filtering with standard deviation w
% correlation is calculated over a fixed time window using definite integrals of gaussian functions
maxspikes = 10000; %fixme?
if iscell(ta)
    
    n = numel(ta);
    if ~exist('ha', 'var')
        
        [ha, hb] = deal(cell(1, n));
        
    end
    
    corr_eachseg = nan(1, n);
    
    assert(all(size(twin) == [n, 2]) && numel(tb) == n, 'invalid inputs');
    
    if numel(w) == 1
        
        w = w * ones(1, n);
        
    else
        
        assert(numel(w) == n, 'invalid inputs');
        
    end
    
    [sab, sa, sb, sasq, sbsq, Ttot] = deal(0); %sums, not expectations
    
    for j = 1:n
        
        if ~(twin(j, 2) > twin(j, 1)) %this skips nans as well
            
            continue;
            
        end
        
        [corr_eachseg(j), eab, ea, eb, easq, ebsq] = point_process_corr_gfilt(ta{j}, tb{j}, w(j), twin(j, :), ha{j}, hb{j});
        
        T    = max(0, diff(twin(j, :)));
        Ttot = Ttot + T;
        sab  = sab  + eab  * T;
        sa   = sa   + ea   * T;
        sb   = sb   + eb   * T;
        sasq = sasq + easq * T;
        sbsq = sbsq + ebsq * T;
        
    end
    
    eab  = sab  / Ttot;
    ea   = sa   / Ttot;
    eb   = sb   / Ttot;
    easq = sasq / Ttot;
    ebsq = sbsq / Ttot;
    
else
    
    if ~exist('ha', 'var') || isempty(ha)
        
        [ha, hb] = deal(ones(numel(ta), 1), ones(numel(tb), 1));
        
    else
        
        ha = ha(:);
        hb = hb(:);
        
    end
    
    ta = ta(:);    tb = tb(:);
    T = diff(twin);
    
    if T <= 0 || numel(ta) > maxspikes || numel(tb) > maxspikes
        
        %warning('invalid data');
        [ea, easq, eb, ebsq, eab] = deal(nan);
        
    else
        
        g_integrals_a = (normcdf_forgfilt(twin(2), ta, w) - normcdf_forgfilt(twin(1), ta, w)) .* ha;   % integral of each gaussian function on the interval
        gprod_integrals_a = gaussfuncproddefint(ta, ta, w ^ 2, w ^ 2, twin(1), twin(2)) .* (ha * ha'); % integral of products of pairs of gaussian functions the interval (including self-products)
        
        ea = sum(g_integrals_a) / T;          % expectation of (sum of weighted gaussian functions) on the interval        
        easq = sum(gprod_integrals_a(:)) / T; % expectation of (sum of weighted gaussian functions) squared on the interval        
        
        g_integrals_b = (normcdf_forgfilt(twin(2), tb, w) - normcdf_forgfilt(twin(1), tb, w)) .* hb;   % integral of each gaussian function on the interval
        gprod_integrals_b = gaussfuncproddefint(tb, tb, w ^ 2, w ^ 2, twin(1), twin(2)) .* (hb * hb'); % integral of products of pairs of gaussian functions the interval (including self-products)
        
        eb = sum(g_integrals_b) / T;          % expectation of (sum of weighted gaussian functions) on the interval        
        ebsq = sum(gprod_integrals_b(:)) / T; % expectation of (sum of weighted gaussian functions) squared on the interval        
        
        gprod_integrals_ab = gaussfuncproddefint(ta, tb, w ^ 2, w ^ 2, twin(1), twin(2)) .* (ha * hb'); % integral of products of pairs of gaussian functions the interval
        eab = sum(gprod_integrals_ab(:)) / T; % expectation of [(sum of weighted gaussian functions)_a * (sum of weighted gaussian functions)_b] on the interval                
        
    end
    
end

va = easq - ea ^ 2;
sda = sqrt(max(0, va));

vb = ebsq - eb ^ 2;
sdb = sqrt(max(0, vb));

cov = eab - ea * eb;
corr = cov / (sda * sdb);


function p = normcdf_forgfilt(x, mu, sigma)
y = (x - mu) / sigma;
p = quickncdf(y);