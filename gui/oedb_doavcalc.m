function oedb_doavcalc(oedb, fileinfo)
mintbefore_lowfr = 1.5; %sec
mintbefore_highfr = 0.3; %sec
frthresh = 1; %Hz
meanresponsewin = [0 0.5]; %sec
peakmaxwin = [0 0.3]; %1ap peak
peakmaxwin2 = [0 0.5]; %2ap burst peak
snrnoise_tbefore = -0.1; %minimum time of bin center before spike
[datasetav, datasetsd, meanresponse, peakheight, peakheight2, snrval, noisesdval, peaksdval, av, av_sd, av_sem, responsedata, naps_isolated, stdiff, peaklat] = deal(cell(1, 0));
%dslist = 1:oedb.ndatasets;
dslist = find(ismember(oedb.datasetnames, {'Simon_rawaverage', 'Wallace_galvo_Tuebingen', 'Wallace_resonance_Bonn', 'Czubayko', 'Packer'}));
[dslist,ok] = listdlg('ListString',oedb.datasetnames,'selectionmode','multiple','name','Calculate averages','InitialValue',dslist);
if ~ok, return; end;
for jj = 1:numel(dslist)
    
    dsind = dslist(jj);
    oerec = fetch_dataset_oerecarray(oedb, fileinfo, dsind);
    
    totalT = arrayfun(@(x) sum(arrayfun(@(d) max(d.t) - min(d.t), x.data)), oerec);
    totalaps = arrayfun(@(x) sum(arrayfun(@(d) sum(d.apcounts), x.data)), oerec);
    fr = totalaps ./ totalT;
    
    mintbefore = mintbefore_lowfr * ones(1, numel(oerec));
    mintbefore(fr > frthresh) = mintbefore_highfr;
    
    [av{jj}, tav, av_sd{jj}, av_sem{jj}, responsedata{jj}, naps_isolated{jj}, naps_total{jj}, fr_total{jj}, T_total{jj}, tevent{jj}, segind{jj}, stdiff{jj}] = avspikeresponse(oerec, mintbefore);
    tav = reshape(tav{1}, [], 1);
    maxspikes = max(cellfun(@numel, av{jj}));
    [datasetav{jj}, datasetsqav, datasetn] = deal(zeros(numel(tav), maxspikes));
    
    for ns = 1:maxspikes
        
        for neuronind = 1:oedb.nneurons(dsind)
            
            if numel(av{jj}{neuronind}) < ns, continue; end
            nextav = av{jj}{neuronind}{ns};
            avmask = ~isnan(nextav);
            datasetn(:,ns) = datasetn(:,ns) + double(avmask);
            datasetav{jj}(avmask, ns) = datasetav{jj}(avmask, ns) + nextav(avmask);
            datasetsqav(avmask,ns) = datasetsqav(avmask,ns) + nextav(avmask) .^ 2;
            
        end
        
        datasetav{jj}(:, ns) = datasetav{jj}(:, ns) ./ datasetn(:,ns);
        datasetsqav(:,ns) = datasetsqav(:,ns) ./ datasetn(:,ns);
        datasetsd{jj}(:,ns) = sqrt(datasetsqav(:,ns) - datasetav{jj}(:, ns) .^ 2);
        
    end
    
    %calculate transient height and area and SNR for single APs
    [meanresponse{jj}, peakheight{jj}, peakheight2{jj}, snrval{jj}, noisesdval{jj}, peaksdval{jj}, peaklat{jj}] = deal(nan(1, oedb.nneurons(dsind)));
    for neuronind = 1:oedb.nneurons(dsind)
        if numel(av{jj}{neuronind}) < 1, continue; end
        meanresponserange = tav > meanresponsewin(1) & tav < meanresponsewin(2);
        meanresponse{jj}(neuronind) = mean(av{jj}{neuronind}{1}(meanresponserange)); %intentionally returns NaN if any elements of the av in this range are NaN.
        peakrange = find(tav > peakmaxwin(1) & tav < peakmaxwin(2));
        peakrange2 = find(tav > peakmaxwin2(1) & tav < peakmaxwin2(2));
        if ~all(isnan(av{jj}{neuronind}{1}(peakrange)))
            
            [peakheight{jj}(neuronind), peakind] = max(av{jj}{neuronind}{1}(peakrange));
            peakt = tav(peakrange(peakind));
            peaklat{jj}(neuronind) = peakt;
            
        else
            
            peakt = nan;
            
        end
        if numel(av{jj}{neuronind}) > 1 && ~all(isnan(av{jj}{neuronind}{2}(peakrange2)))
            
            [peakheight2{jj}(neuronind), ~] = max(av{jj}{neuronind}{2}(peakrange2));
            
        end
        nevents = numel(responsedata{jj}{neuronind}{1});
        [noisef, peakf] = deal(nan(nevents, 1));
        for eventind = 1:nevents
            tobs = responsedata{jj}{neuronind}{1}{eventind}(:,1);
            [~, snr_noiseind] = find(tobs < min(-median(diff(tobs)) * 0.99, snrnoise_tbefore),1,'last');
            if ~isempty(snr_noiseind)
                noisef(eventind) = responsedata{jj}{neuronind}{1}{eventind}(snr_noiseind, 2);
            end
            if isnan(peakt), continue; end
            [tdiff, peakind_event] = min(abs(tobs - peakt));
            if tdiff > median(diff(tobs)) * 1.01, continue; end
            peakf(eventind) = responsedata{jj}{neuronind}{1}{eventind}(peakind_event, 2);
        end
        if nevents > 1
            noisesdval{jj}(neuronind) = std_nonnan(noisef);
            peaksdval{jj}(neuronind) = std_nonnan(peakf);
        end
        snrval{jj}(neuronind) = peakheight{jj}(neuronind) / noisesdval{jj}(neuronind);
    end
    
end
totalmaxspikes = max(cellfun(@(x) size(x, 2), datasetav));
totalav = nan(numel(tav), totalmaxspikes); %across datasets, for each spike count
for ns = 1:totalmaxspikes
    totalavdata = zeros(numel(tav), 0);
    for jj = 1:numel(dslist)
        if size(datasetav{jj},2) >= ns
            totalavdata = [totalavdata datasetav{jj}(:,ns)];
        end
    end
    totalav(:, ns) = mean_nonnan(totalavdata,2);
end

co = [     0         0    1.0000
    0    0.5000         0
    1.0000         0         0
    0    0.7500    0.7500
    0.7500         0    0.7500
    0.7500    0.7500         0
    0.2500    0.2500    0.2500];
col = nan(numel(dslist), 3);
for jj = 1:numel(dslist)
    col(jj,:) = co(mod(jj - 1, size(co, 1)) + 1,:);
end

%1AP vs 2AP peak heights
figure('name', '1ap vs. 2ap peak heights');
for jj = 1:numel(dslist)
    
    plot(peakheight{jj}, peakheight2{jj}, '.', 'color', col(jj,:)); hold on;
    xlabel('1AP peak height');
    ylabel('2AP peak height');
    title(oedb.datasetnames{dslist(jj)})
    
end
%make a plot of basic summary statistics
figure('name','basic statistics');
prevneurons = 0;
ph = [];
for jj = 1:numel(dslist)
    dsind = dslist(jj);
    subplot(5,1,1,'tickdir','out','box','off'); hold on;
    ph(1, end + 1) = plot((1:oedb.nneurons(dsind)) + prevneurons, snrval{jj},'.','color',col(jj,:)); ylabel('SNR');
    subplot(5,1,2,'tickdir','out','box','off'); hold on;
    plot((1:oedb.nneurons(dsind)) + prevneurons, peakheight{jj},'.','color',col(jj,:)); ylabel('Peak h.');
    subplot(5,1,3,'tickdir','out','box','off'); hold on;
    plot((1:oedb.nneurons(dsind)) + prevneurons, noisesdval{jj},'.','color',col(jj,:)); ylabel('noise s.d.');
    subplot(5,1,4,'tickdir','out','box','off'); hold on;
    plot((1:oedb.nneurons(dsind)) + prevneurons, peaksdval{jj},'.','color',col(jj,:)); ylabel('peak s.d.');
    subplot(5,1,5,'tickdir','out','box','off'); hold on;
    plot((1:oedb.nneurons(dsind)) + prevneurons, peaksdval{jj},'.','color',col(jj,:)); ylabel(['mean, ' num2str(meanresponsewin(1))  '-' num2str(meanresponsewin(2))]);
    prevneurons = prevneurons + oedb.nneurons(dsind);
end
subplot(5,1,1);
L = legend(ph, oedb.datasetnames(dslist));
set(findobj(L,'type','text'),'interpreter','none');

maxspikes_plot = 5;
figlist = [];
%first, plot spike responses for each spike count. show each dataset's average on the plot
maxspikes = min(maxspikes_plot, totalmaxspikes);
if numel(dslist) > 1
    figlist(1, end + 1) = figure('name','Dataset comparison');
    for ns = 1:maxspikes
        subplot(maxspikes, 1, ns,'nextplot','add');
        for jj = 1:numel(dslist)
            if size(datasetav{jj},2) >= ns
                nextph = plot(tav, datasetav{jj}(:,ns),'color',col(jj,:));
                if ns == 1
                    ph(jj) = nextph;
                end
            end
        end
        plot(tav, totalav(:,ns),'k','linewidth',2);
        title([num2str(ns) ' APs']);
        if ns == 1
            L = legend(ph, oedb.datasetnames(dslist));
            set(findobj(L,'type','text'),'interpreter','none');
        end
    end
end
%now again, but with individual neurons' averages visible too
maxspikes = min(maxspikes_plot, totalmaxspikes);
if numel(dslist) > 1
    figlist(1, end + 1) = figure('name','Dataset comparison w/ individual neurons');
    for ns = 1:maxspikes
        subplot(maxspikes, 1, ns,'nextplot','add');
        for jj = 1:numel(dslist)
            dsind = dslist(jj);
            for neuronind = 1:oedb.nneurons(dsind)
                if numel(av{jj}{neuronind}) >= ns
                    plot(tav, av{jj}{neuronind}{ns},'color',col(jj,:));
                end
            end
        end
        for jj = 1:numel(dslist)
            if size(datasetav{jj},2) >= ns
                nextph = plot(tav, datasetav{jj}(:,ns),'color',col(jj,:),'linewidth',2);
                if ns == 1
                    ph(jj) = nextph;
                end
            end
        end
        title([num2str(ns) ' APs']);
        if ns == 1
            L = legend(ph, oedb.datasetnames(dslist));
            set(findobj(L,'type','text'),'interpreter','none');
        end
    end
end
%and a figure showing just the single AP responses
if numel(dslist) > 1
    figlist(1, end + 1) = figure('name','Dataset comparison w/ individual neurons (1AP only)');
    ns = 1;
    axes('nextplot','add');
    for jj = 1:numel(dslist)
        dsind = dslist(jj);
        for neuronind = 1:oedb.nneurons(dsind)
            if numel(av{jj}{neuronind}) >= ns
                plot(tav, av{jj}{neuronind}{ns},'color',col(jj,:));
            end
        end
    end
    for jj = 1:numel(dslist)
        if size(datasetav{jj},2) >= ns
            nextph = plot(tav, datasetav{jj}(:,ns),'color',col(jj,:),'linewidth',2);
            if ns == 1
                ph(jj) = nextph;
            end
        end
    end
    title([num2str(ns) ' APs']);
    if ns == 1
        L = legend(ph, oedb.datasetnames(dslist));
        set(findobj(L,'type','text'),'interpreter','none');
    end
end
%next, do plots for each dataset showing the total average and the individual neurons' averages
for jj = 1:numel(dslist)
    dsind = dslist(jj);
    figlist(1, end + 1) = figure('name',oedb.datasetnames{dsind});
    maxspikes = min(maxspikes_plot, size(datasetav{jj},2));
    ax = []; lineobj = []; origlinewidth = []; index = [];
    for ns = 1:maxspikes
        neuronsused = 0;
        ax(ns) = subplot(maxspikes, 1, ns,'nextplot','add');
        for neuronind = 1:oedb.nneurons(dsind)
            if numel(av{jj}{neuronind}) >= ns
                lineobj(1, end+1) = plot(tav, av{jj}{neuronind}{ns},'b');
                index(1, end+1) = neuronind;
                neuronsused = neuronsused + 1;
            end
        end
        plot(tav, datasetav{jj}(:,ns),'color','k','linewidth',2);
        title([num2str(ns) ' APs (n = ' num2str(neuronsused) ' neurons)']);
    end
    names = oedb.neuronnames{dsind}';
    makeavsliders(figlist(1, end), lineobj, index, ax, names)
end

%next, do plots for each neuron showing the average and individual trials
for jj = 1:numel(dslist)
    dsind = dslist(jj);
    for neuronind = 1:oedb.nneurons(dsind)
        figlist(1, end + 1) = figure('name',[oedb.datasetnames{dsind} ' ' oedb.neuronnames{dsind}{neuronind}]);
        maxspikes = min(maxspikes_plot, numel(av{jj}{neuronind}));
        ax = []; lineobj = []; origlinewidth = []; index = []; names = {};
        for ns = 1:maxspikes
            ax(ns) = subplot(maxspikes, 1, ns,'nextplot','add');
            for eventind = 1:numel(responsedata{jj}{neuronind}{ns})
                lineobj(1, end + 1) = plot(responsedata{jj}{neuronind}{ns}{eventind}(:,1), responsedata{jj}{neuronind}{ns}{eventind}(:,2),'b');
                index(1, end + 1) = numel(index) + 1;
                segname = oedb.segmentnames{dsind}{neuronind}{segind{jj}{neuronind}{ns}{eventind}};
                names{end + 1, 1} = [num2str(ns) 'AP, t = ' num2str(tevent{jj}{neuronind}{ns}{eventind}) ', ' segname];
            end
            plot(tav, av{jj}{neuronind}{ns},'k','linewidth',2);
            title([num2str(ns) ' APs (n = ' num2str(numel(responsedata{jj}{neuronind}{ns})) ' events)']);
        end
        if maxspikes > 0
            makeavsliders(figlist(1, end), lineobj, index, ax, names);
        end
    end
end
vn = {'datasetav', 'datasetsd', 'meanresponse', 'peakheight', 'peakheight2', 'snrval', 'noisesdval', 'peaksdval', 'av', 'av_sd', 'av_sem', 'responsedata', 'tav', 'stdiff', 'peaklat', 'tevent'};
for u = 1:numel(vn)
    assignin('base',vn{u}, eval(vn{u}));
end


function makeavsliders(figh, lineobj, index, ax, names)
for j = 1:numel(ax)
    set(ax(j),'units','normalized','tickdir','out','box','off');
    YL = get(ax(j),'ylim'); XL = get(ax(j), 'xlim');
    set(ax(j),'ytick',[0 floor(YL(2) / 0.01) * 0.01]);
    set(ax(j),'xtick',[0 floor(XL(2) / 0.1) * 0.1]);
    p = get(ax(j),'position');
    p(1) = 0.2; p(3) = 0.8;
    set(ax(j),'position',p);
end
origlinewidth = [];
origlinecolor = zeros(0, 3);
for k = 1:numel(lineobj)
    origlinewidth(k) = get(lineobj(k),'linewidth');
    origlinecolor(k, :) = get(lineobj(k),'color');
end
u = uicontrol('Style','listbox','units','normalized','position',[0 0 .15 1],'parent',figh,'value',1,'string',names);
set(u,'callback', @(hObject, eventdata) avlistboxcb(hObject, eventdata, lineobj, index, ax, names, origlinewidth, origlinecolor));
avlistboxcb(u, [], lineobj, index, ax, names, origlinewidth, origlinecolor);


function avlistboxcb(hObject, eventdata, lineobj, index, ax, names, origlinewidth, origlinecolor)
v = get(hObject,'value');
for j = 1:numel(ax)
    ch = get(ax(j),'children');
    selectedobj = intersect(ch, lineobj(index == v));
    avobj = setdiff(findobj(ax(j), 'linewidth', 2), lineobj);
    backobj = setdiff(ch, union(selectedobj, avobj));
    set(ax(j), 'children', [reshape(selectedobj, [], 1); reshape(avobj, [], 1); reshape(backobj, [], 1)]);
end
if isempty(ax), return; end
title(ax(1), names{v},'interpreter','none');
for k = reshape(find(index ~= v),1,[])
    set(lineobj(k), 'linewidth', origlinewidth(k),'color',origlinecolor(k, :));
end
for u = reshape(find(index == v),1,[])
    set(lineobj(u), 'linewidth', 5, 'color', [1 .5 0]);
end