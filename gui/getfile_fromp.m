function [f,p] = getfile_fromp(startp, ext, dialogtitle)
olddir = cd;
try cd(startp); end %#ok<TRYNC>
[f,p] = uigetfile(ext, dialogtitle);
cd(olddir);
if ~isnumeric(f) && p(end) ~= filesep
    p = [p filesep];
end