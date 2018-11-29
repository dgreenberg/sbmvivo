function [dff, baseline] = kinetic2dff_windowed(F,p,w,blgsd)
%converts a sequence of raw fluorescence values F to df/f0 values dff
%baseline is calculated as the percentile value of raw fluorescence over a moving window
%
%simplest version:
%dff = kinetic2dff_windowed(F)
%
%full version:
%[dff, baseline] = kinetic2dff_windowed(F,p,w,blgsd)
%
%INPUTS:
%F - a column vector of fluorescence values, or a matrix with 1 column for each neuron
%p - percentile value (default: 50)
%w - window size, in number of time steps (default: 200). w must be odd, or else it will be increased by 1 automatically. A window size of 20 sec is recommended.
%blgsd - smoothing s.d. for baseline, in time steps (default: 1). An s.d. of 100 ms is recommended.
%
%OUTPUTS:
%dff - dff values, same size as F. defined as 100 * (F - baseline) ./ baseline
%baseline - baseline fluorescence (F0) values
%
%EXAMPLE:
%Suppose we have fluorescence values recorded at 20 Hz, stored in a variable F.
%then we can use the following commands to calculate and plot baseline and df/f0
%
% dt = 1/ 20;
% p = 50;
% w = ceil(20 / dt);
% blgsd = ceil(0.1 / dt);
% [dff, baseline] = kinetic2dff_windowed(F,p,w,blgsd);
% figure;
% t = (1:numel(F))' * dt;
% ax1 = subplot(2,1,1);
% hf = plot(t, F, 'k');
% hold on;
% hb = plot(t, baseline, 'g');
% legend([hf hb], 'Raw fluorescence', 'Baseline fluorescence','location','best');
% xlabel('Time (s)');
% ylabel('Fluorescence (gsv)');
% ax2 = subplot(2,1,2);
% hd = plot(t, dff);
% xlabel('Time (s)');
% ylabel('Fluorescence (\Delta F / F_0)')
% linkaxes([ax1 ax2], 'x');
%
%Department of Brain and Behavior Organization
%Research Institut caesar
%https://www.caesar.de/
%
%last updated 22.08.18 by David Greenberg
if ~exist('p', 'var')
    p = 50;
end
if ~exist('w', 'var')
    w = 200;    
end
if ~exist('blgsd', 'var')
    blgsd = 1;
end
if size(F,1) < 2
    warning('not enough data');
    dff = nan + zeros(size(F));
    baseline = dff;
    return;
end
w = max(2,min(w, floor((size(F,1) - 1) / 2)));
if ~mod(w,2)
    w = w+1;
end
ws = (w-1)/2; %number of bins on one side of baselining window
if ws < 2
    dff = nan + F;
    baseline = nan + F;
    return;
end
n = size(F,1);
m = size(F,2);
dff = zeros(size(F));
if nargout > 1
    bl = zeros(size(F));
end
if nargin < 4
    blgsd = 0;
end
if blgsd > 0
    gr = gauss_filter(F,blgsd);
    [init_sort,init_perm] = sort(gr(1:ws,:));
    %keyboard
    baseline = kinetic2dff_windowed_mainloop(gr,p,w,init_perm,init_sort);
else
    [init_sort,init_perm] = sort(F(1:ws,:));
    baseline = kinetic2dff_windowed_mainloop(F,p,w,init_perm,init_sort);
end

dff = 100 * (F - baseline) ./ baseline; %convert to df/f0