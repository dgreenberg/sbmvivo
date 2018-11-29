function stop = canonlinfit_opt_outputfcn(x, optimValues, state, figh)
stop = false;
if ~strcmpi(state, 'iter'), return; end
%record history of params + error
allp = getappdata(figh, 'allp');
allp = cat(2, allp, x(:));
setappdata(figh, 'allp', allp);
alle = getappdata(figh, 'alle');
stop = numel(alle) > 0 && any(alle <= optimValues.fval);
alle = cat(2, alle, optimValues.fval);
setappdata(figh, 'alle', alle);