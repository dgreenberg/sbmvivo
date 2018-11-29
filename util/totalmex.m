function totalmex(includedirs)
if ~exist('includedirs','var')
    includedirs = {};
end
olddir = cd;
p = path;
f = find(p == ';' | p == ':');
z = 0;
pp = cell(1, numel(f));
for k=1:length(f)
    pp{k} = p(z + 1:(f(k) - 1));
    if pp{k}(end) ~= filesep
        pp{k} = [pp{k} filesep];
    end
    z = f(k);
end
c99flag = false;
if isunix
    mcc = mex.getCompilerConfigurations;
    for k = 1:numel(mcc)
        if strcmpi(mcc(k).Name, 'GNU C')
            c99flag = true;
        end
    end
end
for k = 1:length(pp)
    if ~isempty(strfind(pp{k},'\toolbox\')) %don't mess with matlab toolboxes
        continue;
    end
    cd(pp{k});
    d = [dir('*.c'); dir('*.cpp')];
    if ~isempty(d)
        disp(['Directory ' cd]);
    end    
    for j = 1:length(d)
        [~, dnoext, ~] = fileparts(d(j).name);
        if ismember(exist([pp{k} dnoext '.' mexext],'file'), [2 3])
            continue; %file has already been compiled
        end
        estring = ['mex ' d(j).name];
        if c99flag
            estring = [estring ' CFLAGS="\$CFLAGS -std=c99"']; %#ok<AGROW> %prevent compilation errors from lines starting with // coments                    
        end
        for u = 1:numel(includedirs)
            estring = [estring ' -I' includedirs{u}]; %#ok<AGROW>
        end
        disp(['Compiling ' d(j).name]);
        eval(estring);
    end    
end
cd(olddir);