function y = gaussfuncproddefint(mu1, mu2, V1, V2, a, b)
%definite integral of product of Gaussian functions
%y = gaussfuncproddefint(mu1, mu2, V1, V2, a, b)
%
%y = \integral_{a}^b \mathbb{N}(x, \mu_1, V_1) \mathbb{N}(x, \mu_2, V_2) dx
%  = \mathbb{N}(\mu_1, \mu_2, V_1 + V2) \integral_{a}^b \mathbb{N}(x, \mu, V) dx
%where
%V = 1 / (1 / V_1 + 1 / V_2)
%and
%\mu = V * (\mu_1 / V_1 + \mu_2 / V_2)
%so we have
%y = \mathbb{N}(\mu_1, \mu_2, V_1 + V2) * (\Phi((b - \mu) / sqrt(V)) - \Phi((a - \mu) / sqrt(V)))

%equivalent but slower: z = normpdf(0, bsxfun(@minus, mu1(:), mu2(:)'), sigma); %normalization factor
sigma = sqrt(bsxfun(@plus, V1(:), V2(:)'));
z = quicknpdf(bsxfun(@minus, mu1(:), mu2(:)') ./ sigma) ./ sigma;

V = 1 ./ bsxfun(@plus, 1 ./ V1(:), 1 ./ V2(:)');
mu = V .* bsxfun(@plus, (mu1(:) ./ V1(:)), (mu2(:) ./ V2(:))');
sd = sqrt(V);

%equivalent but slower: y = z .* (normcdf(b, mu, sd) - normcdf(a, mu, sd));
y = z .* (quickncdf((b - mu) ./ sd) - quickncdf((a - mu) ./ sd));