function data = oerecdata_extract_framerange(data, ii)
data.f = data.f(ii);
if ~isempty(data.f_removed)
    
    data.f_removed = data.f_removed(ii);
    
end
data.f_mask = data.f_mask(ii);
data.t = data.t(ii);
data.frameindex = data.frameindex(ii);
if data.spiketimes_present
    %keep all the spikes...why not?
    %if k > 1
    %   data.spiketimes_window(1) = max(data.spiketimes_window(1), data.t(ii(1)));
    %end
    %if k < numel(s)
    %	data.spiketimes_window(2) = min(data.spiketimes_window(2), data.t(ii(end)));
    %end
    %data.spiketimes = data.spiketimes(data.spiketimes >= data.spiketimes_window(1) & data.spiketimes <= data.spiketimes_window(2));
end
if data.apcounts_present
    
    data.apcounts = data.apcounts(ii);
    
end