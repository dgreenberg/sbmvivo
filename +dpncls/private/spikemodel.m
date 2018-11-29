function [amp, model, bl, blockrej, initamp] = spikemodel(rawdata, dt, regwin, fitwin, minth, minindamp, minpeakdff, minpeakderiv, minpeakderivhigh, highthresh, tau, blwin, blgsd, ...
    minregpeak, pval, nearnoregthresh, minstradle,mindepbestamp,reblsd,reblsd2,reblwin,mindepmean,minpeakdoubderiv,areawin,mindepdoubderiv,afterspikewin,maxblocksize)
assign_defaults('nearnoregthresh',0.27,'\afterspikewin','fitwin','maxblocksize',inf);

model = nan + zeros(size(rawdata));
amp = model;
bl = model;
initamp = amp;

if any(isnan(rawdata))
    gf = ~isnan(rawdata);
    s = unique([min(find(gf)); find(diff(gf) == 1) + 1]);
    e = unique([max(find(gf)); find(diff(gf) == -1)]);
    for k = 1:length(s)
        [amp(s(k):e(k)), model(s(k):e(k)), bl(s(k):e(k)), blockrej(s(k):e(k)), initamp(s(k):e(k))] = spikemodel(rawdata(s(k):e(k)), dt, regwin, fitwin, minth, minindamp, minpeakdff, minpeakderiv, ...
            minpeakderivhigh, highthresh, tau, blwin, blgsd, minregpeak, pval, nearnoregthresh, minstradle,mindepbestamp,reblsd,reblsd2,reblwin,...
            mindepmean,minpeakdoubderiv,areawin,mindepdoubderiv,afterspikewin,maxblocksize);
    end
    return;
end

ndata = length(rawdata);
regabins = round(regwin(2) / dt);
regbbins = round(-regwin(1) / dt);
fitwinbins = round(fitwin / dt);  
areawinbins = floor(areawin / dt);
taubins = tau / dt;
wbins = ceil(blwin / dt);
blgsdbins = blgsd / dt;
nearnoregbins = ceil(nearnoregthresh / dt);

fitcomp = exp(-(0:fitwinbins - 1) / taubins)';

dff = kinetic2dff_windowed(rawdata, pval, wbins, blgsdbins);
passed = calc_filter(dff, dt, fitwinbins, fitcomp, minindamp, minpeakdff, minpeakderiv, regbbins, regabins, minpeakderivhigh, highthresh, minstradle,minpeakdoubderiv);

tested = ~isnan(passed);
tested(1:regbbins) = 0;
tested(max(1,end-regabins + 1):end) = 0;
reblsdbins = reblsd / dt;
reblsd2bins = reblsd2 / dt;
reblwinbins = ceil(reblwin / dt);
aspikebins = round(afterspikewin / dt);
blockrej = zeros(size(rawdata));

n = partition_vector(passed, fitwinbins - 1);
if ~isempty(n) && any(n > 0)
    disp(['points per cluster: ' num2str(histc(n(find(passed == 1)),1:max(n))')]);
    pause(0.05)
else
    amp(tested & ~(amp > 0)) = 0;
    model = makemodel(fitcomp, amp(amp > 0), find(amp > 0), ndata);
    [dff, bl, mask] = rebaseline(amp, rawdata, wbins, fitwinbins, passed, reblsdbins, reblsd2bins, aspikebins);
    bl(isnan(amp)) = nan;
    return;
end

for i = 1:max(n)
    f = find(n == i);    
    data = dff(f);
    mask = passed(f) == 1;
    if sum(mask) > maxblocksize
        blockrej(f) = 1;
        passed(f) = nan;
        continue; %leave amp and model as nan
    end
    [amp(f), model(f)] = deconv_mincoef_leastsquares_tree(data, fitcomp, minth, mask);
    %f(amp(f) > 0)
end
%amp(1:regbbins) = nan;
%amp(max(1,end-regabins + 1):end) = nan;
tested = ~isnan(passed);
initamp = amp;

poss = amp > 0;
while(1)
    prev_poss = poss;
    [dff, bl, blpoints] = rebaseline(amp, rawdata, reblwinbins, fitwinbins, passed, reblsdbins, reblsd2bins, aspikebins);  %adjust baseline and df/f
    passed = calc_filter(dff, dt, fitwinbins, fitcomp, minindamp, minpeakdff, minpeakderiv, regbbins, regabins, minpeakderivhigh, highthresh, minstradle,minpeakdoubderiv);
    
    n = partition_vector(poss, fitwinbins - 1);
    for i = 1:max(n)
        f = find(n == i);
        data = dff(f);
        mask = poss(f) == 1;
        amp(f) = deconv_mincoef_leastsquares_tree(data, fitcomp, minth, mask);
    end                
    poss(amp==0 | isnan(amp)) = 0;
    
    if ~any(poss)
        break;
    elseif all(poss == prev_poss)
        f = find(amp > 0);
        model = makemodel(fitcomp, amp(f), f, ndata);
        [exreglist,exbestcomp,exdoubderiv,exarea] = calc_postfilt(poss, dff, model, taubins, minregpeak, regabins, regbbins, fitcomp, amp, areawinbins);

        f = find(amp > 0);
        g = find(diff(f) <= nearnoregbins);
        g = unique([g; g + 1]);
        f(g) = []; %nearby transients make each other immune to this test
        exreglist(g) = [];
        if ~isempty(f) && any(exreglist < minregpeak)
            [garb, minind] = min(exreglist);
            poss(f(minind)) = 0;
        elseif any(exbestcomp < mindepbestamp)
            f = find(poss);
            [garb, minind] = min(exbestcomp);
            poss(f(minind)) = 0;
        elseif any(exdoubderiv < mindepdoubderiv)
            f = find(poss);
            [garb, minind] = min(exdoubderiv);
            poss(f(minind)) = 0;
        elseif any(exarea < mindepmean)
            f = find(poss);
            [garb, minind] = min(exarea);
            poss(f(minind)) = 0;
        else
            break; %all tests passed
        end
        poss = poss & (passed == 1);
        amp(~poss & ~isnan(amp)) = 0;
    end
end
amp(tested & ~(amp > 0)) = 0;
model = makemodel(fitcomp, amp(amp > 0), find(amp > 0), ndata);
bl(isnan(amp)) = nan;

function m = makemodel(template, amp, fm, n)
m = zeros(n,1);
m(fm) = amp;
m = filter(template,1,m);

function passed = calc_filter(dff, dt, fitwinbins, fitcomp, minindamp, minpeakdff, minpeakderiv, regbbins, regabins, ...
    minpeakderivhigh, highthresh, minstradle,minpeakdoubderiv)
deriv = [nan; diff(dff)];
doubderiv = [nan; nan; dff(3:end) - dff(1:end-2)];
stradle = [nan; dff(3:end) - dff(1:end-2); nan];
bestamp = nan + zeros(size(dff));
passed = zeros(size(dff));
noreg = logical(ones(size(dff)));
fitcompdenom = 1 / (fitcomp' * fitcomp);
gf = ~isnan(dff);
s = unique([min(find(gf)); find(diff(gf) == 1) + 1]);
e = unique([max(find(gf)); find(diff(gf) == -1)]);
for j = 1:length(s)
    fitpoints = s(j):e(j) - fitwinbins + 1; %time bins where we will test for the presence of transients
    %fitpoints = fitpoints(find(fitpoints > fitwinbins));
    if ~isempty(fitpoints)
        fitnpoints = length(fitpoints);
        datablock = dff(repmat((fitpoints(1):fitpoints(1) + fitwinbins - 1)', 1, fitnpoints) + repmat(0:fitnpoints - 1,fitwinbins,1));
        bestamp(fitpoints) = (fitcomp' * datablock) * fitcompdenom;
    end    
    regresspoints = s(j) + regbbins:e(j) - regabins; %time bins where we will test for the presence of transients
    noreg(regresspoints) = 0;
end    
passed(stradle > minstradle & doubderiv > minpeakdoubderiv & bestamp >= minindamp & dff >= minpeakdff & ((deriv >= minpeakderivhigh & dff >= highthresh) | deriv >= minpeakderiv)) = 1;
passed(isnan(stradle) | isnan(bestamp) | isnan(dff) | isnan(deriv) | noreg | isnan(doubderiv)) = nan;

function [dff, bl, mask] = rebaseline(amp, rawdata, wbins, fitwinbins, passed, reblsdbins, reblsd2bins, aspikebins)
ws = ceil((wbins - 1) / 2);
mask = ~isnan(rawdata) & ~isnan(passed);
for jj = find(amp > 0)'
    mask(max(1,jj - 1):min(end,jj + aspikebins - 1 + 1)) = 0;
end
bl = masked_baseline(rawdata, mask, ws, reblsdbins);
gf = ~isnan(bl);
s = unique([min(find(gf)); find(diff(gf) == 1) + 1]);
e = unique([max(find(gf)); find(diff(gf) == -1)]);
for j = 1:length(s)
    bl(s(j):e(j)) = gauss_filter(bl(s(j):e(j)), reblsd2bins);
end
dff = 100 * (rawdata - bl) ./ bl;

function [exregcomp,exbestcomp,exdoubderiv,exarea] = calc_postfilt(poss, dff, model, taubins, minregpeak, regabins, regbbins, fitcomp, amp, areawinbins)
fitwinbins = length(fitcomp);
regnbins = regbbins + 1 + regabins;
regcomps = [[zeros(regbbins,1); exp(-(0:regabins) / (taubins))'] ones(regnbins,1) [zeros(regbbins,1); ones(regabins + 1,1)]];
fitcomp = exp(-(0:fitwinbins - 1) / taubins)';
fitcompdenom = 1 / (fitcomp' * fitcomp);

fulldiff = dff - model;

exregcomp = zeros(sum(poss),1); 
fm = find(poss);
for jj = 1:length(fm)
    if fm(jj) > regbbins && fm(jj) + regabins <= length(dff) && ~any(isnan(fulldiff(fm(jj) - regbbins:fm(jj) + regabins)))
        exdiff = fulldiff(fm(jj) - regbbins:fm(jj) + regabins);
        L = min(fitwinbins, regabins + 1);
        exdiff(regbbins + 1:regbbins + L) = exdiff(regbbins + 1:regbbins + L) + fitcomp(1:L) * amp(fm(jj));
        r = regcomps \ exdiff;
        if r(2) < 0 && r(3) < 0
            r = regcomps(:,1) \ exdiff;
        elseif r(2) < 0 %don't let the offset help, only hurt
            r = regcomps(:,[1 3]) \ exdiff;
            if r(2) < 0
                r = regcomps(:,1) \ exdiff;
            end
        elseif r(3) < 0
            r = regcomps(:,[1 2]) \ exdiff;
            if r(2) < 0
                r = regcomps(:,1) \ exdiff;
            end
        end
        exregcomp(jj) = r(1);
        exbestcomp(jj) = fitcomp(1:L)' * exdiff(regbbins + 1:regbbins + L) * fitcompdenom;        
        exarea(jj) = sum(exdiff(regbbins + 1:regbbins + 1 + min(areawinbins,regabins))) / min(areawinbins,regabins);
        exderiv(jj) = exdiff(regbbins + 1) - exdiff(regbbins);
        exdoubderiv(jj) = exdiff(regbbins + 1) - exdiff(regbbins - 1);
    else
        exregcomp(jj) = nan;        
        exbestcomp(jj) = nan;
        exarea(jj) = nan;
        exderiv(jj) = nan;
        exdoubderiv(jj) = nan;
    end
end
return;

ii = 2;
[amp, model] = spikemodel(mudata.kin.raw{ii}, mudata.kin.dt(ii), mudata.params.regresswin, mudata.params.fitwin, mudata.params.minth, mudata.params.minindh, ...
    mudata.params.minpeakval, mudata.params.minderiv, mudata.params.tau, mudata.params.blwin, mudata.params.blgsd, mudata.params.minregpeak, mudata.params.pval, 0.37);

minpeakderivhigh = 2;
minpeakderiv = 4.5;
highthresh = 9;
kk = 2;
q = bobdata.roi.kdata(:,kk);
q(~bobdata.image.goodframes | bobdata.roi.cellrej(:,kk) == 1) = nan;
[amp, model] = spikemodel(q, 0.096, mudata.params.regresswin, mudata.params.fitwin, mudata.params.minth, mudata.params.minindh, ...
    mudata.params.minpeakval, minpeakderiv, minpeakderivhigh, highthresh, mudata.params.tau, mudata.params.blwin, mudata.params.blgsd, mudata.params.minregpeak, mudata.params.pval, 0.37);