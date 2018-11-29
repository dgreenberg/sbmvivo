function [f,p] = putfile_fromp(startp, ext, dialogtitle, filename)
olddir = cd;
try cd(startp); end; %#ok<TRYNC>
if nargin  > 3 && ~isempty(filename)
    [f,p] = uiputfile(ext, dialogtitle, filename);
else
    [f,p] = uiputfile(ext, dialogtitle);
end
cd(olddir);
if ~isnumeric(f) && p(end) ~= filesep
    p = [p filesep];
end