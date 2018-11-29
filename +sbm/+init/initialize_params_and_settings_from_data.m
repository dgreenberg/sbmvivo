function [P,V] = initialize_params_and_settings_from_data( ...
    F,dt,P,V,indicatorstring,nA2D)

[P, V] = sbm.init.assign_default_settings_and_params(P,V,indicatorstring); %fill in any nonexistent fields for which we have a default value

if ~isfield(V, 'settingsindex') || isempty(V.settingsindex)
    V.settingsindex = ones(1, numel(F));
end
V.settingsindex = checkvectorinput(V.settingsindex, numel(F), 'V.settingsindex');
if ~isfield(V, 'sessionindex') || isempty(V.sessionindex)
    V.sessionindex = 1:numel(F);
end
V.sessionindex = checkvectorinput(V.sessionindex, numel(F), 'V.sessionindex');
if ~isfield(V, 'focalplaneindex') || isempty(V.focalplaneindex)
    V.focalplaneindex = ones(1, numel(F));
end
V.focalplaneindex = checkvectorinput(V.focalplaneindex, numel(F), 'V.focalplaneindex');

P = sbm.init.initialize_obs_params(P, V, F, dt, nA2D); %choose alpha, gain, zeta to match some simple properties of the fluorescence traces if they're not already available

assert(~isfield(P, 'kd_ex') || numel(P.kd_ex) < 2, 'multiple extrusion mechanisms are not yet implemented for particle filtering')