function gdata = gauss_filter(data,s,L)
%gdata = gauss_filter(data,s,L)
%gaussian filters each column of data
%s and L are in bins, #bins==size(data,1)
%s is standard deviation
%L is width of one side fo the filter
%if L is not given, then
%L = min(ceil(10*s),size(data,1))
%for real data only
if size(data,1) < 2
    gdata = data;
    return;
end
if nargin < 3
    L = min(ceil(10*s),size(data,1)-1);
elseif L >= size(data,1)
    L = size(data,1)-1; %for speed, at no cost to accuracy
end
if ~s
    gdata = data;
    return;
end
L = round(L);
if L<=0, L=1; end;

x = (-L:L)';
f = exp(-(x.^2)/(2*s^2))/(s*sqrt(2*pi)); %guassian filter
f = f/sum(f); %normalize

r = zeros(size(data,1) + 2*L + 1,1);
r(1:L+1) = f(L+1:end);
r(end-L+1:end) = f(1:L);
fftr = fft(r);

gdata = zeros(size(data));
for j = 1:size(data,2) %for loop sacrifices a small amount of speed in return for memory.  with many short columns, this may not be desirable
    fft_padded = fft([data(1,j) * ones(L,1); data(:,j); data(end,j) * ones(L+1,1)]);
    %fft_padded = fft([flipud(data(1:L,j)); data(:,j); flipud(data(end-L:end,j))]);
    conv = real(ifft(fft_padded .* fftr));
    gdata(:,j) = conv(L + 1:end - L - 1,1);
end