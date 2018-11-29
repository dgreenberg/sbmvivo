function [seq, beq] = equilibriumstates(params, opts)
if params.c0 == 0
    
    seq = zeros(numel(params.kon), 1);  % calcium-free state omitted
    beq = zeros(numel(params.Btot), 1);
    return;
    
end

Ka = params.kon ./ params.koff;
AKterms = cumprod([1, Ka * params.c0]);  % Adair-Klotz
seq = params.S * AKterms(2:end) / sum(AKterms);  % calcium-free state omitted

Ka_B = 1 ./ params.kd_B;
beq = params.Btot .* Ka_B .* params.c0 ./ (1 + Ka_B * params.c0);

seq = reshape(seq, [], 1);
beq = reshape(beq, [], 1);