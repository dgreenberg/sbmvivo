function [c, s, b] = simulate(n_newton_iterations, cprev, sprev, bprev, stepsize, c0, g, kd_ex, kon, koff, kon_B, koff_B, Btot)
nbindingsteps = numel(koff);
nbuffers = numel(kon_B);
nparticles = numel(cprev);

npz = zeros(nparticles, 1);
npo = ones(nparticles, 1);

[c, s, b] = deal(cprev, sprev, bprev);

npz_npo_koffmat = [npz npo * koff];

if ~isempty(kd_ex)

    maxex = g * kd_ex;
    c_leak = maxex * c0 / (c0 + kd_ex);
    
end

if n_newton_iterations == 0  % Forward Euler
    
    kon_times_s = bsxfun(@times, kon, s(:, 1:end - 1));
    
    rate = bsxfun(@times, kon_times_s, c) - bsxfun(@times, koff, s(:, 2:end));
    
    B_unbound = bsxfun(@minus, Btot, b);
    
    rate_B = bsxfun(@times, kon_B, c) .* B_unbound - bsxfun(@times, koff_B, b);
    
    ds_dt = [npz, rate] - [rate, npz];
    db_dt = rate_B;
    
    if isempty(kd_ex)
        
        dc_dt = -g * (c - c0) - sum(rate, 2) - sum(rate_B, 2);
        
    else
    
        dc_dt = -maxex * c ./ (c + kd_ex) + c_leak - sum(rate, 2) - sum(rate_B, 2);
        
    end
    
    c = c + dc_dt * stepsize;
    s = s + ds_dt * stepsize;
    b = b + db_dt * stepsize;
    
    return;
    
end

for i = 1:n_newton_iterations  % Backward Euler
    
    kon_times_s = bsxfun(@times, kon, s(:, 1:end - 1));
    
    rate = bsxfun(@times, kon_times_s, c) - bsxfun(@times, koff, s(:, 2:end));
    
    B_unbound = bsxfun(@minus, Btot, b);
    
    rate_B = bsxfun(@times, kon_B, c) .* B_unbound - bsxfun(@times, koff_B, b);
    
    ds_dt = [npz, rate] - [rate, npz];
    db_dt = rate_B;
    
    if isempty(kd_ex)
        
        dc_dt = -g * (c - c0) - sum(rate, 2) - sum(rate_B, 2);
        M_c_c = 1.0 + stepsize * (g + sum(kon_times_s, 2) + sum(bsxfun(@times, kon_B, B_unbound), 2));
                
    else
        
        dc_dt = -maxex * c ./ (c + kd_ex) + c_leak - sum(rate, 2) - sum(rate_B, 2);
        M_c_c = 1.0 + stepsize * (maxex * kd_ex ./ (c + kd_ex) .^ 2 + sum(kon_times_s, 2) + sum(bsxfun(@times, kon_B, B_unbound), 2));
        
    end
    
    M_c_s = stepsize * ([bsxfun(@times, c, kon) npz] - npz_npo_koffmat);
    M_s_c = stepsize * ([kon_times_s npz] - [npz kon_times_s]);
    
    M_s_s_diag = 1.0 + stepsize * ([bsxfun(@times, c, kon) npz] + npz_npo_koffmat);
    M_s_s_diagp1 = -stepsize * npo * koff;
    M_s_s_diagm1 = -stepsize * bsxfun(@times, c, kon);
    
    M_c_b = -stepsize * bsxfun(@plus, bsxfun(@times, kon_B, c), koff_B);
    M_b_c = -stepsize * bsxfun(@times, kon_B, B_unbound);
    M_b_b = 1 - M_c_b;
    
    z_c = cprev - c + stepsize * dc_dt;
    z_s = sprev - s + stepsize * ds_dt;
    z_b = bprev - b + stepsize * db_dt;
    
    % forward pass sets below-diagonal to zero and diagonal to one. first row and column and above-diagonal all change but remain dense
    z_s(:, 1) = z_s(:, 1) ./ M_s_s_diag(:, 1);
    M_s_c(:, 1) = M_s_c(:, 1) ./ M_s_s_diag(:, 1);
    M_s_s_diagp1(:, 1) = M_s_s_diagp1(:, 1) ./ M_s_s_diag(:, 1);
    for j = 2:nbindingsteps + 1
        
        a = M_s_s_diagm1(:, j - 1);
        denom = M_s_s_diag(:, j) - a .* M_s_s_diagp1(:, j - 1);
        z_s(:, j) = (z_s(:, j) - a .* z_s(:, j - 1)) ./ denom;
        M_s_c(:, j) = (M_s_c(:, j) - a .* M_s_c(:, j - 1)) ./ denom;
        if j < nbindingsteps + 1
            
            M_s_s_diagp1(:, j) = M_s_s_diagp1(:, j) ./ denom;
            
        end
        
    end
    
    % backward pass sets above-diagonal and M_c_s to zero
    for j = nbindingsteps + 1:-1:1
        
        if j > 1
            
            z_s(:, j - 1)   = z_s(:, j - 1)   - M_s_s_diagp1(:, j - 1) .*   z_s(:, j);
            M_s_c(:, j - 1) = M_s_c(:, j - 1) - M_s_s_diagp1(:, j - 1) .* M_s_c(:, j);
            
        end
        z_c   = z_c   - M_c_s(:, j) .*   z_s(:, j);
        M_c_c = M_c_c - M_c_s(:, j) .* M_s_c(:, j);
        
    end
    
    % forward pass only for buffers, sets M_b_b to one M_c_b to zero
    for j = 1:nbuffers
        
        z_b(:, j)   = z_b(:, j)   ./ M_b_b(:, j);
        M_b_c(:, j) = M_b_c(:, j) ./ M_b_b(:, j);
        
        z_c   =   z_c - M_c_b(:, j) .*   z_b(:, j);
        M_c_c = M_c_c - M_c_b(:, j) .* M_b_c(:, j);
        
    end
    
    z_c = z_c ./ M_c_c; % set M_c_c to one
    
    % forward pass eliminates M_s_c
    for j = 1:nbindingsteps + 1
        
        z_s(:, j) = z_s(:, j) - M_s_c(:, j) .* z_c;
        
    end
    
    % forward pass eliminates M_b_c
    for j = 1:nbuffers
        
        z_b(:, j) = z_b(:, j) - M_b_c(:, j) .* z_c;
        
    end
    
    c = c + z_c;
    s = s + z_s;
    b = b + z_b;
    
end