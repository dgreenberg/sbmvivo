function print_likelihood(lik, prestring)
if nargin < 2
    prestring = '';
end

fprintf('%slog P(F) %8.6g\n', prestring, sum(cat(2,lik.logmarglik)));
    