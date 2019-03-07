function [config_dir] = config_dir()
config_dir = fullfile(char(java.lang.System.getProperty('user.home')),'.sbmvivo');
end

