function oerec2c2sfile(filename, oerec, preprocess)
%oerec2c2sfile(filename, oerec, preprocess)
if ~exist('preprocess', 'var') || isempty(preprocess)
    
    preprocess = true;

end

data = {};
for u = 1:numel(oerec)
    
    [f, it, ~, ns, ~, ~, st, ~] = extract_data_fromoerec(oerec(u));
    
    for s = 1:numel(oerec(u).data)
        
        tdiff = diff(it{s});
        dt = median(tdiff); %sec FIXME check units
        assert(all(abs(tdiff - dt) < dt / 100), 'data must be evenly spaced in time');
        
        t0 = it{s}(1);
        
        data{1, end + 1} = orderfields(struct( ...
            'calcium', f{s}, ...
            'fps', 1 / dt, ...
            'cell_num', u ...
            )); %#ok<AGROW>
        
        if ~isempty(st{s})
           
            data{end}.spike_times = (st{s} - t0) * 1000; %spike times seem to be relative to the first fluorescence measurement, based on c2s.py's preprocess function. in ms
            
        elseif ~isempty(ns{s})
            
            data{end}.spikes = ns{s};
            
        end
        
        data{end} = orderfields(data{end});
        
    end
end

if preprocess
    
    save([filename '_raw.mat'], 'data', '-mat');    
    commandstr = sprintf('c2s preprocess %s_raw.mat %s -v 0', filename, filename);    
    assert(system(commandstr) == 0, 'command failed');
    delete([filename '_raw.mat']);
    
else
    
    save(filename,'data','-mat');
    
end