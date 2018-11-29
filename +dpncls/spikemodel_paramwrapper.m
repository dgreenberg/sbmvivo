function [amp, model, bl, blockrej, initamp] = spikemodel_paramwrapper(rawdata, dt, params, maxblocksize)
assign_defaults('maxblocksize',inf);
desiredversion = 1.0;
if params.version ~= desiredversion
    error(['this version of spikemodel_paramwrapper expects a parameter structure of version ' num2str(desiredversion) '; input version was ' num2str(params.version)]);
end
[amp, model, bl, blockrej, initamp] = spikemodel(rawdata, dt, params.regwin, params.fitwin, params.minth, params.minindamp, params.minpeakdff, ...
    params.minpeakderiv, params.minpeakderivhigh, params.highthresh, params.tau, params.blwin, params.blgsd, ...
    params.minregpeak, params.pval, params.nearnoregthresh, params.minstradle, params.mindepbestamp, params.reblsd, params.reblsd2, ...
    params.reblwin, params.mindepmean, params.minpeakdoubderiv, params.areawin, params.mindepdoubderiv, params.afterspikewin, maxblocksize);