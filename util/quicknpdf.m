function p = quicknpdf(x)
p = exp(-0.5 * x .^ 2) / sqrt(2 * pi);