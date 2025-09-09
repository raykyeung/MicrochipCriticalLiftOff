function h = B2K_daviolinplot_withOverlay02(Y,varargin)
% B2K_daviolinplot_withOverlay02  daviolinplot + full-violin overlay support
%
%   Exactly the same functionality as daviolinplot.m, plus:
%     'overlay'  logical (false by default). If true, draws all groups
%                at each condition as full violins, centered at the same
%                x-position, with FaceAlpha = violinalpha.
%
%   Usage:
%     h = B2K_daviolinplot_withOverlay02(Y, 'groups', G, ..., 'overlay', true);

h = struct;
p = inputParser;

%% 1) Original Name/Value options
addOptional(p, 'groups', []);          
addOptional(p, 'violin', 'half');       
addOptional(p, 'bins', 10);
addOptional(p, 'colors', get(gca,'colororder'));
addOptional(p, 'violinalpha', 1);
addOptional(p, 'violinwidth', 1);
addOptional(p, 'smoothing','default');
addOptional(p, 'box', 2);
addOptional(p, 'boxcolors', 'k');
addOptional(p, 'whiskers', 1);
addOptional(p, 'scatter', 0);
addOptional(p, 'scattersize', 20);
addOptional(p, 'scattercolors', 'k');
addOptional(p, 'scatteralpha', 1);
addOptional(p, 'jitter', 0);
addOptional(p, 'jitterspacing', 1);
addOptional(p, 'outliers', 1);
addOptional(p, 'outfactor', 1.5);
addOptional(p, 'outsymbol', 'k*');
addOptional(p, 'boxalpha', 1);
addOptional(p, 'boxspacing', 1);
addOptional(p, 'boxwidth', 1);
addOptional(p, 'linkline', 0);
addOptional(p, 'withinlines', 0);
addOptional(p, 'xtlabels', []);
addOptional(p, 'legend', []);

%% 2) New overlay flag
addOptional(p, 'overlay', false, @(x) islogical(x) || ismember(x,[0,1]));

parse(p, varargin{:});
confs = p.Results;

%% A) Build grouping index Gi and num_groups
if ~isempty(confs.groups)
    [Gi,~,Gv] = grp2idx(confs.groups);
    num_groups = numel(Gv);

elseif iscell(Y)
    num_groups = numel(Y);
    Gi = [];
    y  = [];
    for g = 1:num_groups
        y  = [y;    Y{g}];
        Gi = [Gi; g*ones(size(Y{g},1),1)];
    end
    Y = y;

else
    Gi = ones(size(Y,1),1);
    num_groups = 1;
end

%% B) If overlay, force full violins
if confs.overlay
    confs.violin = 'full';
end

%% C) Condition positions
if any(size(Y)==1)
    Y = Y(:);
    cpos = 1;
else
    cpos = 1:size(Y,2);
end
num_locs = numel(cpos);

%% D) Compute group positions gpos & box_width
if confs.overlay
    % all groups centered at each condition
    gpos     = repmat(cpos, [num_groups,1]);
    box_width = 0.1 * confs.boxwidth;
else
    % original spacing logic
    if num_locs == 1
        gpos      = (1:num_groups)';
        box_width = 0.1*confs.boxwidth;
        cpos      = gpos;
    else
        if num_groups == 1
            gpos      = cpos;
            box_width = 0.1*confs.boxwidth;
        else
            box_width = 0.2/(num_groups+1) * confs.boxwidth;
            loc_sp    = 4*box_width * confs.boxspacing;
            gpos = [];
            for g = 1:num_groups
                gpos = [gpos;
                        cpos + (g-(num_groups+1)/2)*(box_width + loc_sp)]; %#ok<AGROW>
            end
        end
    end
end

h.gpos = gpos;
h.cpos = cpos;

%% E) Pre-compute colors for boxes & scatters
if ~strcmp(confs.boxcolors,'same')
    bcol = confs.boxcolors;
end
if ~strcmp(confs.scattercolors,'same')
    scol = confs.scattercolors;
end

if exist('scol','var') && strcmp(scol,'k')
    ecol = 'w';
else
    ecol = 'k';
end

if exist('bcol','var') && strcmp(bcol,'k')
    mcol = 'w';
else
    mcol = 'k';
end

%% F) Plotting loop
for g = 1:num_groups
    % percentiles & IQR
    pt  = prctile(Y(Gi==g,:), [2 9 25 50 75 91 98]);
    if size(pt,1)==1, pt=pt'; end
    IQR = pt(5,:) - pt(3,:);
    
    % box coordinates
    y25       = reshape([pt(3,:);pt(3,:)], 1, []);
    y75       = reshape([pt(5,:);pt(5,:)], 1, []);
    x1        = [gpos(g,:)-box_width/2; gpos(g,:)-box_width/2];
    x2        = [gpos(g,:)+box_width/2; gpos(g,:)+box_width/2];
    box_ycor  = [y75; y25];
    box_medcor= reshape([pt(4,:);pt(4,:)],1,[]);
    
    % boxplot & whisker x coords
    if strcmp(confs.violin,'full')
        box_xcor = reshape([x1;x2],2,[]);
        whi_xcor = [gpos(g,:); gpos(g,:)];
    else
        if     confs.box==3, bpos = -15;
        elseif confs.box==2, bpos = 2;
        else                 bpos = 0.5;
        end
        box_xcor = reshape([x1;x2],2,[]) - box_width/bpos;
        whi_xcor = [gpos(g,:); gpos(g,:)] - box_width/bpos;
    end
    
    % per-condition plotting
    for k = 1:num_locs
        data_vals = Y(Gi==g,k);
        
        % density estimate
        if strcmp(confs.smoothing,'default')
            [f,xi] = ksdensity(data_vals);
        else
            [f,xi] = ksdensity(data_vals,'Bandwidth',confs.smoothing);
        end
        f = confs.violinwidth * (f./max(f)) * (21*box_width/(num_groups+7));
        
        % draw violin
        if strcmp(confs.violin,'full')
            h.ds(k,g) = fill([ f, -fliplr(f)] + gpos(g,k), ...
                              [ xi, fliplr(xi)], ...
                              confs.colors(g,:));
        else
            if confs.box==3
                h.ds(k,g) = fill(1.3*f + gpos(g,k), xi, confs.colors(g,:));
            else
                h.ds(k,g) = fill(   f + gpos(g,k), xi, confs.colors(g,:));
            end
        end
        set(h.ds(k,g),'FaceAlpha',confs.violinalpha);
        hold on
        
        % draw boxplot
        if confs.box~=0
            win_k = (1:2) + 2*(k-1);
            Yy    = box_ycor(1:2,win_k);
            Xx    = box_xcor(1:2,win_k);
            
            if strcmp(confs.boxcolors,'same')
                facecol = confs.colors(g,:);
            else
                facecol = bcol;
            end
            
            h.bx(k,g) = fill([Xx(:,1)' Xx(1,:) Xx(:,2)' Xx(2,[2,1])], ...
                             [Yy(:,1)' Yy(1,:) Yy(:,2)' Yy(2,:)], ...
                             facecol);
            set(h.bx(k,g),'FaceAlpha',confs.boxalpha);
            
            % median
            h.md(k,g) = line(Xx(1,:), box_medcor(win_k), ...
                             'Color', mcol, 'LineWidth', 2);
            
            % whiskers
            if confs.whiskers==1
                ol = data_vals < (pt(3,k)-confs.outfactor*IQR(k));
                ou = data_vals > (pt(5,k)+confs.outfactor*IQR(k));
                whi_ycor(:,1,k) = [min(data_vals(~ol)), pt(3,k)];
                whi_ycor(:,2,k) = [max(data_vals(~ou)), pt(5,k)];
                h.wh(k,g,:) = plot( ...
                  whi_xcor(:,k), whi_ycor(:,1,k), 'k-', ...
                  whi_xcor(:,k), whi_ycor(:,2,k), 'k-', ...
                  'LineWidth',1.5);
            end
        end
        
        % jitter / histogram for scatter
        if     confs.jitter==0
            xdata = gpos(g,k)*ones(numel(data_vals),1);
        elseif confs.jitter==1
            xdata = gpos(g,k)*ones(numel(data_vals),1) + ...
                    (box_width/1.5)*(0.5-rand(numel(data_vals),1));
        else  % jitter==2
            [N,E] = histcounts(data_vals, confs.bins);
            bin_w = E(2)-E(1);
            bm    = E(1:end-1) + bin_w/2;
            hx=[]; hy=[];
            for i=1:numel(N)
                for j=1:N(i)
                    hy(end+1) = bm(i); %#ok<AGROW>
                    hx(end+1) = (-j+0.8)*confs.jitterspacing/(8*numel(gpos)) + gpos(g,k); %#ok<AGROW>
                end
            end
            data_vals = hy;
            xdata     = hx;
        end
        
        % shift half-violin
        if strcmp(confs.violin,'half') && confs.jitter~=2
            xdata = xdata - box_width/2;
        end
        
        % scatter offset
        if confs.scatter==1
            xdata = xdata - 1.3*box_width;
        elseif confs.scatter==2
            xdata = xdata - box_width/4;
        elseif confs.scatter==3
            xdata = xdata + box_width;
        end
        
        scdata{g}(:,:,k) = [xdata(:), data_vals(:)]; %#ok<AGROW>
        
        % outliers
        if confs.outliers==1
            ox = data_vals < (pt(3,k)-confs.outfactor*IQR(k)) | ...
                 data_vals > (pt(5,k)+confs.outfactor*IQR(k));
            h.ot(k,g) = scatter(xdata(ox), data_vals(ox), ...
                                confs.scattersize, confs.outsymbol);
        else
            ox = false(size(data_vals));
        end
        
        % scatter
        if confs.scatter~=0
            if strcmp(confs.scattercolors,'same')
                mfc = confs.colors(g,:);
                mec = 'k';
            else
                mfc = scol;
                mec = ecol;
            end
            h.sc(k,g) = scatter(xdata(~ox), data_vals(~ox), ...
                                confs.scattersize, ...
                                'MarkerFaceColor', mfc, ...
                                'MarkerEdgeColor', mec, ...
                                'MarkerFaceAlpha', confs.scatteralpha);
        end
    end
    
    % link-lines
    if confs.linkline==1
        h.ln(g) = line(gpos(g,:), pt(4,:), ...
                       'Color', confs.colors(g,:), ...
                       'LineStyle','-.', ...
                       'LineWidth',1.5);
    end
    
    % within-lines
    if confs.withinlines==1
        for s = 1:size(scdata{g},1)
            h.wl(g) = plot( scdata{g}(s,1,:), scdata{g}(s,2,:), ...
                            'Color',[0.8 0.8 0.8]);
            uistack(h.wl(g),'bottom');
        end
    end
end

% move link-lines behind
if confs.linkline==1
    uistack([h.ln],'bottom');
end

% legend
if ~isempty(confs.legend)
    h.lg = legend(h.ds(1,:), confs.legend);
end

% ticks & labels
set(gca,'XTick',cpos,'XTickLabels',cpos,'Box','off');
if ~isempty(confs.xtlabels)
    set(gca,'XTickLabels',confs.xtlabels,'XTick',cpos);
end

% x-limits
xlim([gpos(1)-3*box_width, gpos(end)+3*box_width]);

end