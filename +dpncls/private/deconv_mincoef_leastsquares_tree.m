function [amp, model, ntests] = deconv_mincoef_leastsquares_tree(data, template, mincoef, mask, minmax_thresh, max_updates, reg_thresh)
%template is required to be an exponential decay starting at amlitude one, for now

%we might improve the search by testing for contradictions in tree
%branches, e.g. when a point is assumed to have a trans and then later the
%error could be improved by removing that transient from minmodel and
%similarly for silence etc. this should return infinite error

%might be desirable to aim for and use segmentation/recursion

%regressions could be sped up by precalculating and using the template's
%autocorrelation function somehow

%global ww
%ww = 0;
%global rc
%rc = 0;
assign_defaults('minmax_thresh',1,'max_updates',30,'reg_thresh',0.01);

data = reshape(data,[],1);
amp = zeros(size(data));
model = zeros(size(data));

n = length(data);
r = length(template);
amp(n - r + 1:n) = nan;

if ~exist('mask','var')
    mask = true(size(data));
    mask(n - r + 2:n) = 0;
    maskinput = 0;
else
    mask = reshape(mask,[],1);
    maskinput = 1;
end
fm = find(mask);

if any(isnan(data))
    ind = repmat(fm', r, 1) + repmat((0:r - 1)',1,length(fm));
    datablock  = data(ind);
    naninblock = any(isnan(datablock),1);
    mask(naninblock) = 0;
    fm = find(mask);
    if any(naninblock) && maskinput
        warning('Input mask contains time points with NaN data in the template window, ignoring these points');
    end
end

if isempty(fm)
    return;
end

templateinvnorm = 1 / (template' * template);
maxamp = max(mincoef, regressamp_ind(data, template, fm));
minamp = zeros(size(fm));
confirmed = false(size(fm));
unsolved = true(size(fm));
interactmat = abs(repmat(fm,1,length(fm)) - repmat(fm',length(fm),1)) < r;
interactmat(1:size(interactmat,1) + 1:end) = 0; %don't need to know that a point interacts with itself
indexmat = repmat(fm', r, 1) + repmat((0:r-1)',1,length(fm));
next2check = (1:length(fm))';
minnext = 1;
maxmodel = makemodel(template, maxamp, fm, length(data));
minmodel = makemodel(template, minamp, fm, length(data));
[minamp, maxamp, confirmed, unsolved, minmodel, maxmodel, active] = update_minmax(data, minamp, maxamp, confirmed, mincoef, template, ...
    max_updates, minmax_thresh, minnext, indexmat, next2check, minmodel, maxmodel, interactmat, unsolved, templateinvnorm);
maskout = false(size(mask));
part_ind = partition_vector(mask, r - 1);
ntests = 0;
for k = 1:max(part_ind)
    f = find(part_ind == k); %indices into data and mask
    g = find(part_ind(fm) == k); %indices into fm
    subactive = active(g);
    g = g(subactive);
    submask = find(mask(f));
    submask = submask(subactive);
    [nextamp, err, nextfm, nextntests] = best_fit(data(f), template, mincoef, submask, minmax_thresh, max_updates, minamp(g), maxamp(g), confirmed(g), reg_thresh, interactmat(g,g), indexmat(:,g) - f(1) + 1, unsolved(g), templateinvnorm, minmodel(f), maxmodel(f));
    amp(f(nextfm)) = nextamp;
    maskout(f(nextfm)) = 1;
    ntests = ntests + 1;
end
fmout = find(maskout);
model = makemodel(template, amp(fmout), fmout, length(data));

function m = makemodel(template, amp, fm, n)
m = zeros(n,1);
m(fm) = amp;
m = filter(template,1,m); %USE THIS for speed

function amp = regressamp_ind(data, template, fm)
datablock = data(repmat((0:length(template) - 1)', 1, length(fm)) + repmat(fm', length(template), 1));
amp = (template' * datablock)' / (template' * template);

function [minamp, maxamp, confirmed, unsolved, minmodel, maxmodel, active] = update_minmax(data, minamp, maxamp, confirmed, mincoef, template, max_updates, minmax_thresh, minnext, indexmat, next2check, minmodel, maxmodel, interactmat, unsolved, templateinvnorm)
%global ww vv
%ww = ww + 1;
iter = 1;
%keyboard
while(~isempty(next2check) && iter < max_updates)
    if minnext
        d = data - maxmodel;
        newminamp = maxamp(next2check) + (template' * d(indexmat(:,next2check)) * templateinvnorm)';
        
        changed = newminamp - minamp(next2check) >= minmax_thresh & newminamp >= mincoef / 2; %also throws out cases where minamp would be rounded down
        if ~any(changed)
            break;
        end
        next2check = next2check(changed);
        newminamp = newminamp(changed);
        
        newminamp(newminamp < mincoef) = mincoef; %round up where needed
        %update model
        for k = 1:length(next2check)
            i = next2check(k);
            minmodel(indexmat(:,i)) = minmodel(indexmat(:,i)) + (newminamp(k) - minamp(i)) * template;            
        end
        minamp(next2check) = newminamp;
        
        confirmed(next2check) = minamp(next2check) >= mincoef;
        unsolved(next2check) = maxamp(next2check) - minamp(next2check) >= minmax_thresh;

        minnext = 0;
        next2check = find(any(interactmat(:,next2check),2) & unsolved); %vector of indices
    else
        d = data - minmodel;
        newmaxamp = minamp(next2check) + (template' * d(indexmat(:,next2check)) * templateinvnorm)';
        %round up and down where needed
        rounddown = newmaxamp < mincoef / 2;
        roundup = (newmaxamp < mincoef) & ~rounddown;
        newmaxamp(rounddown & ~confirmed(next2check)) = 0; %note that confirmed points can't be rounded down, though I'm not sure if that could actually happen
        newmaxamp(roundup) = mincoef;
        
        changed = maxamp(next2check) - newmaxamp >= minmax_thresh;
        if ~any(changed)
            break;
        end
        next2check = next2check(changed);
        newmaxamp = newmaxamp(changed);
        %update model
        for k = 1:length(next2check)
            i = next2check(k);
            maxmodel(indexmat(:,i)) = maxmodel(indexmat(:,i)) + (newmaxamp(k) - maxamp(i)) * template;
        end
        maxamp(next2check) = newmaxamp;
        
        maxmodel(maxmodel < 0) = 0;
        unsolved(next2check) = maxamp(next2check) - minamp(next2check) >= minmax_thresh;

        minnext = 1;
        next2check = find(any(interactmat(:,next2check),2) & unsolved);
    end
    iter = iter + 1;
end
active = maxamp > 0;
%vv(ww) = iter;
%iter

function [amp, err] = regress_constrained(data, template, fm, mincoef, reg_thresh)
%points with unrestrained regression coefficients between mincoef / 2 and
%mincoef can be fixed at the start to mincoef

%global rc z
%rc = rc + 1;
%could use minamp and  maxamp
r = length(template);
npoints = length(data);
ncomps = length(fm);
regcomp = zeros(npoints, ncomps);
for k = 1:ncomps
    regcomp(fm(k):fm(k) + r - 1, k) = template;
end
covmat = regcomp' * regcomp;
amp = covmat \ (regcomp' * data); %unconstrained regression
if ~any(amp < mincoef)
    err = sum((regcomp * amp - data).^2);
    return;
end
denom = 1 / (template' * template);
%keyboard
regiter = 0;
while(1) %begin gradient descent. the set is convex.
    regiter = regiter + 1;
    prevamp = amp;
    amp = max(mincoef, amp);
    f = find(amp > mincoef);
    g = find(amp == mincoef);

    tmpdata = data - sum(regcomp(:,g) * mincoef,2);
    tmpamp = covmat(f,f) \ (regcomp(:,f)' * tmpdata);
    h = find(tmpamp < mincoef); %indices into f
    if ~isempty(h)
        [adjustratio, ind] = min((amp(f(h)) - mincoef) ./ (amp(f(h)) - tmpamp(h))); %the ratio of the distance to the minimum along which we can advance before hitting the limit imposed by mincoef
        amp(f) = amp(f) + (tmpamp - amp(f)) * adjustratio;
        amp(f(h)) = max(mincoef, amp(f(h))); %to avoid rounding errors
    else
        amp(f) = tmpamp;
    end

    tmpdata = data - sum(regcomp * amp,2);
    
    for k = g'
        coefadj = regcomp(:,k)' * tmpdata * denom; %least-squares optimal change in coefficient k, which is currently valued at mincoef
        if coefadj > 0
            amp(k) = amp(k) + coefadj;
            tmpdata = tmpdata - regcomp(:,k) * coefadj;
        end
    end
    if all(abs(amp - prevamp) < reg_thresh)
        break; %convergence
    end
end
err = sum((regcomp * amp - data).^2);
%z(rc) = regiter;

function [amp, err, fm, ntests] = best_fit(data, template, mincoef, fm, minmax_thresh, max_updates, minamp, maxamp, confirmed, reg_thresh, interactmat, indexmat, unsolved, templateinvnorm, minmodel, maxmodel)
ndata = length(data);
npoints = length(fm);
nambig = sum(~confirmed & maxamp > 0);
if ~nambig
    ntests = 1;
    [amp, err] = regress_constrained(data, template, fm, mincoef, reg_thresh);
    return;
end

silentbranched = false(nambig,1); %we will follow down and resolve silent branches before following down transient branches
trans_const = logical([0 1]);
silence_const = logical([1 0]);

nposs_hyp = zeros(nambig, 2); %first column for silence, second column for transient
nposs_tree = zeros(nambig + 1, 2); %first column for silence, second column for transient
nposs_max = 2^nambig;
nposs_tree(1, silence_const) = nposs_max;

tree2hyp_indmat = repmat((1:npoints)',[1, nambig + 1, 2]);
tree2hyp_indmat_tmp = tree2hyp_indmat;

tree2hyp_indmat_model = repmat((1:ndata)',[1,nambig + 1, 2]);
tree2hyp_indmat_model_tmp = tree2hyp_indmat_model;

powersof2 = 2.^(0:nambig);

%could keep all of these in a single giant array to increase speed
err_tree = zeros(nambig + 1,2);
finalamp_tree = zeros(npoints,nambig + 1,2);
minamp_tree = repmat(minamp,[1,nambig + 1,2]);
maxamp_tree = repmat(maxamp,[1,nambig + 1,2]);
confirmed_tree = repmat(confirmed,[1,nambig + 1,2]);
unsolved_tree = repmat(unsolved,[1,nambig + 1,2]);
active_tree = true(npoints,nambig + 1,2);
minmodel_tree = repmat(minmodel,[1,nambig + 1,2]);
maxmodel_tree = repmat(maxmodel,[1,nambig + 1,2]);

minamp_hyp = repmat(minamp,[1,nambig,2]);
maxamp_hyp = repmat(maxamp,[1,nambig,2]);
confirmed_hyp = repmat(confirmed,[1,nambig,2]);
unsolved_hyp = repmat(unsolved,[1,nambig,2]);
active_hyp = true(npoints,nambig,2);
minmodel_hyp = repmat(minmodel,[1,nambig,2]);
maxmodel_hyp = repmat(maxmodel,[1,nambig,2]);

ntests = 0;
focuslevel = 1;
focuspos = silence_const; %above the top of the tree, there abides only eternal silence
while(1) %branch the node with focus, or solve it and propogate the result up the tree
    if all(confirmed_tree(unsolved_tree(:, focuslevel, focuspos), focuslevel, focuspos)) %all transient times are known
        %solve for transient heights
        [amp, err_tree(focuslevel,focuspos)] = regress_constrained(data, template, indexmat(1,confirmed_tree(:,focuslevel,focuspos))', mincoef, reg_thresh);
        ntests = ntests + 1;
        finalamp_tree(:, focuslevel, focuspos) = 0;
        finalamp_tree(confirmed_tree(:,focuslevel,focuspos), focuslevel, focuspos) = amp;
        %keyboard
        if focuspos(1) %just solved a silent case
            nposs_tree(focuslevel, silence_const) = 1;
            focuspos = trans_const; %now move the focus to the transient case on the same tree level
            continue;
        else %just solved a transient case
            node2resolve = find(silentbranched, 1, 'last' ); %deepest silent-hypothesis node which is branched. this can be the top node
            superbranch = node2resolve+1:focuslevel; %focus levels for which silent-hypothesis nodes shall be compared to the present transient-hypothesis node
            [best_silent_err, best_silent_choice] = min(err_tree(superbranch, silence_const));

            %could probably do this without a flow-control statement:
            if best_silent_err < err_tree(focuslevel, focuspos) %the newly solved transient case is inferior to all silent cases above it and below the previous branched silent node
                bestfocuslevel = superbranch(best_silent_choice);
                besterr = best_silent_err;
                bestfocuspos = silence_const;
            else %the newly solved transient case is superior to all silent cases above it and below the previous branched silent node
                bestfocuslevel = focuslevel;
                besterr = err_tree(focuslevel, focuspos);
                bestfocuspos = trans_const;
            end
            %now resolve the preceding branched silent-hypothesis node
            err_tree(node2resolve, silence_const) = besterr;
            finalamp_tree(:, node2resolve, silence_const) = finalamp_tree(:, bestfocuslevel, bestfocuspos);
            nposs_tree(node2resolve, silence_const) = 1; %we set this value to zero when we branched out, now we set it to one as we collapse the branch

            if node2resolve == 1
                ampnonzero = finalamp_tree(:,1,silence_const) > 0;
                fm = indexmat(1,ampnonzero); %indicies into data corresponding to transient times
                amp = finalamp_tree(ampnonzero,1,silence_const);
                err = besterr;
                return;
            else
                focuspos = trans_const;
            end
            silentbranched(node2resolve) = 0;
            focuslevel = node2resolve;
            %keyboard
            nposs_tree(superbranch, :) = 0;
        end
    else %some transient times are not known, so branch
        branch_cand = find(active_tree(:,focuslevel,focuspos) & ~confirmed_tree(:,focuslevel,focuspos)); %places where the occurence of a transient is unknown
        n_branch_cand = length(branch_cand);

        tree2hyp_indmat_tmp(:,1:n_branch_cand,:) = tree2hyp_indmat(:,1:n_branch_cand,:) + (focuslevel - 1) * npoints + focuspos(2) * npoints * (nambig + 1);
        tree2hyp_indmat_model_tmp(:,1:n_branch_cand,:) = tree2hyp_indmat_model(:,1:n_branch_cand,:) + (focuslevel - 1) * ndata + focuspos(2) * ndata * (nambig + 1);
                
        maxamp_hyp(:,branch_cand,:) = maxamp_tree(tree2hyp_indmat_tmp(:,1:n_branch_cand,:));
        minamp_hyp(:,branch_cand,:) = minamp_tree(tree2hyp_indmat_tmp(:,1:n_branch_cand,:));
        confirmed_hyp(:,branch_cand,:) = confirmed_tree(tree2hyp_indmat_tmp(:,1:n_branch_cand,:));
        unsolved_hyp(:,branch_cand,:) = unsolved_tree(tree2hyp_indmat_tmp(:,1:n_branch_cand,:));
        active_hyp(:,branch_cand,:) = active_tree(tree2hyp_indmat_tmp(:,1:n_branch_cand,:));
        
        maxmodel_hyp(:,branch_cand,:) = maxmodel_tree(tree2hyp_indmat_model_tmp(:,1:n_branch_cand,:));
        minmodel_hyp(:,branch_cand,:) = minmodel_tree(tree2hyp_indmat_model_tmp(:,1:n_branch_cand,:));
        
        for k = branch_cand' %index into minamp, maxamp, etc.            
            next2check = find(interactmat(:,k) & unsolved_tree(:, focuslevel, focuspos)); %note that inactive points are considered solved
            %some initializations might be vectorized outside the for-loop

            %consider the case where no transient occurs
            maxamp_hyp(k, k, silence_const) = 0;
            confirmed_hyp(k, k, silence_const) = 0;
            active_hyp(k, k, silence_const) = 0;
            unsolved_hyp(k, k, silence_const) = 0;
            maxmodel_hyp(indexmat(:,k), k, silence_const) = max(0, maxmodel_hyp(indexmat(:,k), k, silence_const) - template * maxamp_tree(k, focuslevel, focuspos));            
            
            minnext = 1;
            [minamp_hyp(:,k,silence_const), maxamp_hyp(:,k,silence_const), confirmed_hyp(:,k,silence_const), unsolved_hyp(:,k,silence_const), ...
                minmodel_hyp(:,k,silence_const), maxmodel_hyp(:,k,silence_const), active_hyp(:,k,silence_const)] = update_minmax(...
                data, minamp_hyp(:,k,silence_const), maxamp_hyp(:,k,silence_const), confirmed_hyp(:,k,silence_const), mincoef, template, max_updates, minmax_thresh, ...
                minnext, indexmat, next2check, minmodel_hyp(:,k,silence_const), maxmodel_hyp(:,k,silence_const), interactmat, unsolved_hyp(:,k,silence_const), templateinvnorm);

            nposs_hyp(k,silence_const) = powersof2(1 + sum(active_hyp(:, k, silence_const) & ~confirmed_hyp(:, k, silence_const)));

            %consider the case where a transient occurs
            minamp_hyp(k, k, trans_const) = mincoef;
            confirmed_hyp(k, k, trans_const) = 1;
            unsolved_hyp(k, k, trans_const) = maxamp_hyp(k, k, trans_const) - mincoef >= minmax_thresh;
            minmodel_hyp(indexmat(:,k), k, trans_const) = minmodel_hyp(indexmat(:,k), k, trans_const) + template * mincoef;
            minnext = 0;
            
            [minamp_hyp(:,k,trans_const), maxamp_hyp(:,k,trans_const), confirmed_hyp(:,k,trans_const), unsolved_hyp(:,k,trans_const), ...
                minmodel_hyp(:,k,trans_const), maxmodel_hyp(:,k,trans_const), active_hyp(:,k,trans_const)] = update_minmax(...
                data, minamp_hyp(:,k,trans_const), maxamp_hyp(:,k,trans_const), confirmed_hyp(:,k,trans_const), mincoef, template, max_updates, minmax_thresh, ...
                minnext, indexmat, next2check, minmodel_hyp(:,k,trans_const), maxmodel_hyp(:,k,trans_const), interactmat, unsolved_hyp(:,k,trans_const), templateinvnorm);

            nposs_hyp(k,trans_const) = powersof2(1 + sum(active_hyp(:, k, trans_const) & ~confirmed_hyp(:, k, trans_const)));
        end
        besttotalposs = min(sum(nposs_hyp(branch_cand,:),2)); %find the branch point candidate whose hypotheses have the least number of possibilities
        bestcandind = floor(median(find(sum(nposs_hyp(branch_cand,:),2) == besttotalposs)));
        bestcand = branch_cand(bestcandind);

        %update nposs_tree here
        silentbranched(focuslevel) = focuspos(1);
        silentbranched(focuslevel + 1) = 0;

        %these operations assign values to both branches simultaneously:
        maxamp_tree(:,focuslevel + 1,:) = maxamp_hyp(:,bestcand,:);
        minamp_tree(:,focuslevel + 1,:) = minamp_hyp(:,bestcand,:);
        confirmed_tree(:,focuslevel + 1,:) = confirmed_hyp(:,bestcand,:);
        unsolved_tree(:,focuslevel + 1,:) = unsolved_hyp(:,bestcand,:);
        active_tree(:,focuslevel + 1,:) = active_hyp(:,bestcand,:);
        maxmodel_tree(:,focuslevel + 1,:) = maxmodel_hyp(:,bestcand,:);
        minmodel_tree(:,focuslevel + 1,:) = minmodel_hyp(:,bestcand,:);

        nposs_tree(focuslevel, focuspos) = 0;
        nposs_tree(focuslevel + 1,:) = nposs_hyp(bestcand,:);
           
        focuslevel = focuslevel + 1; %now solve the just-branched-out nodes
        focuspos = silence_const; %deal with the silent case first
    end
end