function oedb = update_oedb_version(oedb, fileinfo, oedb_version, min_oedb_version)
assert(oedb.version <= oedb_version, 'odb file version is unrecognized, software update required');
if oedb.version < min_oedb_version
    error(['file versions earlier than ' num2str(oedb_version) ' are not supported']); %this should change with release or whenever we have significant results in odb format
end

nalgs = numel(oedb.algnames);

% if oedb.version < XXX
%     
% end
oedb.version = oedb_version;
oedb = orderfields(oedb);
%fixme we should have updaters for opts and params for each alg etc.


function s = newstatfield(s, fn)
for j = 1:numel(s)
    
    for k = 1:numel(fn)
        
        if isfield(s, fn{k}), continue; end
        
        s(j).(fn{k}) = nan;
        s(j).([fn{k} '_sd']) = nan;
        
    end
    
end

s = orderfields(s);