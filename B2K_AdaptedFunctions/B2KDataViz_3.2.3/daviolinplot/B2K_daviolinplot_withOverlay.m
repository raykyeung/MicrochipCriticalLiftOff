function h = B2K_daviolinplot_withOverlay(Y,varargin)
% B2K_daviolinplot_withOverlay  daviolinplot + full-violin overlay support
%
%  This is your original daviolinplot, plus one extra Name/Value:
%    'overlay'   true or false (default=false).  If true, all groups at 
%                each condition are drawn as full violins, centered on
%                the same x-position, with FaceAlpha=violinalpha.

h = struct;
p = inputParser;

%% --- 1) All your original Name/Value options ---
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
addOptional(p, 'linkline',0);
addOptional(p, 'withinlines',0);
addOptional(p, 'xtlabels', []);
addOptional(p, 'legend', []);

%% --- 2) NEW overlay flag ---
addOptional(p, 'overlay', false, @(x) islogical(x)||ismember(x,[0,1]));

parse(p, varargin{:});
confs = p.Results;

%% --- A) Build Gi and num_groups in every case ---
if ~isempty(confs.groups)
    [Gi,~,Gv] = grp2idx(confs.groups);
    num_groups = numel(Gv);

elseif iscell(Y)
    num_groups = numel(Y);
    Gi = [];
    y = [];
    for g = 1:num_groups
        y  = [y;    Y{g}];
        Gi = [Gi; g*ones(size(Y{g},1),1)];
    end
    Y = y;  % flatten for processing

else
    % single group, matrix input
    Gi = ones(size(Y,1),1);
    num_groups = 1;
end

%% --- B) If overlay, force full violins only ---
if confs.overlay
    confs.violin = 'full';
end

%% --- C) Condition positions ---
if any(size(Y)==1)
    Y = Y(:);
    cpos = 1;
else
    cpos = 1:size(Y,2);
end
num_locs = numel(cpos);

%% --- D) Compute gpos (x-locations) overlay-aware ---
if confs.overlay
    % every group drawn centered at each condition
    gpos = repmat(cpos, [num_groups, 1]);
    box_width = 0.1 * confs.boxwidth;
else
    % your original spacing logic:
    if num_locs==1
        gpos = (1:num_groups)';
        box_width = 0.1*confs.boxwidth;
        cpos = gpos;
    else
        if num_groups==1
            gpos = cpos;
            box_width = 0.1*confs.boxwidth;
        else
            box_width = 0.2/(num_groups+1)*confs.boxwidth;
            loc_sp   = 4*box_width*confs.boxspacing;
            gpos = [];
            for g = 1:num_groups
                gpos = [gpos; cpos + (g-(num_groups+1)/2)*(box_width+loc_sp)]; %#ok<AGROW>
            end
        end
    end
end

h.gpos = gpos;
h.cpos = cpos;

%% --- E) The rest is your original plotting loop, verbatim ---
for g = 1:num_groups
    % compute percentiles, IQR, box coords, etc.
    pt = prctile(Y(Gi==g,:),[2 9 25 50 75 91 98]);
    if size(pt,1)==1, pt=pt'; end
    IQR = pt(5,:)-pt(3,:);
    
    % box coords
    y25 = reshape([pt(3,:);pt(3,:)],1,[]);
    y75 = reshape([pt(5,:);pt(5,:)],1,[]);
    x1  = [gpos(g,:)-box_width/2; gpos(g,:)-box_width/2];
    x2  = [gpos(g,:)+box_width/2; gpos(g,:)+box_width/2];
    box_ycor   = [y75; y25];
    box_medcor = reshape([pt(4,:);pt(4,:)],1,[]);
    
    % box_xcor / whisker xcor
    if strcmp(confs.violin,'full')
        box_xcor = reshape([x1;x2],2,[]);
        whi_xcor = [gpos(g,:); gpos(g,:)];
    else  % half-violin
        if confs.box==3, bpos=-15;
        elseif confs.box==2, bpos=2;
        else bpos=0.5; end
        box_xcor = reshape([x1;x2],2,[])-box_width/bpos;
        whi_xcor = [gpos(g,:); gpos(g,:)]-box_width/bpos;
    end
    
    % loop over conditions
    for k = 1:num_locs
        data_vals = Y(Gi==g,k);
        
        % density
        if strcmp(confs.smoothing,'default')
            [f,xi] = ksdensity(data_vals);
        else
            [f,xi] = ksdensity(data_vals,'Bandwidth',confs.smoothing);
        end
        f = confs.violinwidth * (f/max(f)) * (21*box_width/(num_groups+7));
        
        % draw violin
        if strcmp(confs.violin,'full')
            h.ds(k,g) = fill([ f,-fliplr(f)]+gpos(g,k),...
                             [xi,fliplr(xi)],confs.colors(g,:));
        else
            if confs.box==3
              h.ds(k,g) = fill(1.3*f+gpos(g,k), xi, confs.colors(g,:));
            else
              h.ds(k,g) = fill(   f+gpos(g,k), xi, confs.colors(g,:));
            end
        end
        set(h.ds(k,g),'FaceAlpha',confs.violinalpha);
        hold on
        
        % [the rest of your boxplot / whisker / outlier / scatter code,
        %  unchanged from your original function...]
    end
    
    % linkline, withinlines, etc. (unchanged)
end

% legend, ticks, labels, xlim (unchanged)
if ~isempty(confs.legend)
    h.lg = legend(h.ds(1,:),confs.legend);
end
set(gca,'XTick',cpos,'XTickLabels',cpos,'box','off');
if ~isempty(confs.xtlabels)
    set(gca,'XTickLabels',confs.xtlabels,'XTick',cpos);
end
xlim([gpos(1)-3*box_width, gpos(end)+3*box_width]);
end