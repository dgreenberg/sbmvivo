function n = strnmatchingchars(s1, s2)
m = min(numel(s1), numel(s2));
n = sum(cumprod(double(s1(1:m) == s2(1:m))));