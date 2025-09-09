function hText = B2KTextAlignedToLineLogLogPlot(ax,x,y,txt,m,varargin)
%% B2KTextAlignedToLineLogLogPlot.m - [Function] Align text annotation to power-law linear line in log-log scale
%
% Author: Raymond Yeung
% Release date: 2025
% E-mail: ryeun003@ucr.edu
% B2K Group, Dept. of Bioengineering, Univ. of California, Riverside
% Victor G. J. Rodgers Dept. of Bioengineering, Univ. of California, Riverside
% William H. Grover, Dept. of Bioengineering, Univ. of California, Riverside
% Philip L. Brisk Dept. of Computer Science, Univ. of California, Riverside
%
% [STATUS] -
% [CURRENT FUNC]
% -
% [WORKING ON]
% -
% [BUGS]
% -
% [STILL NEED]
% -
%
% Procedure:
% 1) Fit or provide log-log slope
% 2) Measure axes' pixel aspect ratio (since figure ratio may skew the visual angle)
% 3) Compute on-screen rotation angle
% 4) Place text with rotation angle
%   % LOGLOGTEXT   place & rotate 'txt' at data‐point (x,y) on loglog axes 'ax'
%   %   'm' is the log–log slope (e.g. 3/7); varargin passes other text() args.
%   fig = ancestor(ax,'figure');
%   oF = fig.Units; fig.Units='pixels'; fpos=fig.Position; fig.Units=oF;
%   oA = ax.Units;  ax.Units='pixels'; apos=ax.Position;  ax.Units=oA;
%   ratio = apos(4)/apos(3);
%   ang   = atan(m*ratio)*180/pi;
%   h = text(ax, x, y, txt, 'Rotation',ang, varargin{:});
% end

% B2KTEXTALIGNEDTOLINELOGPLOTPLOT
%   Place TEXT at (x,y) in data units on the log‑log axes AX,
%   rotated so that it follows the slope m = d(log10 y)/d(log10 x).
%   Works correctly even if AX is inside a tiledlayout.
%
%   H = B2KTextAlignedToLineLogLogPlot(...)
%   returns the handle to the text object.
%
%   Additional name/value pairs (e.g. 'Interpreter','latex') go in varargin.

%     % 1) Find the “pixel container” (tiledlayout or figure)
%     parent = ax.Parent;
%     if isa(parent, 'matlab.graphics.layout.TiledChartLayout')
%         container = parent;
%     else
%         container = ancestor(ax, 'figure');
%     end
% 
%     % 2) Grab container size in pixels; save-and-restore to temporarily switch Units to pixels
%     oldU = container.Units; 
%     container.Units = 'pixels'; 
%     contPos = container.Position; 
%     container.Units = oldU;
% 
%     % 3) Grab axes size in pixels (relative to that container)
%     oldU = ax.Units; 
%     ax.Units = 'pixels'; 
%     axPos = ax.Position; 
%     ax.Units = oldU;
% 
%     % 4) Compute pixel‑space slope = m*(height/width)
%     pixelRatio = axPos(4) / axPos(3);
%     angDeg     = atan( m * pixelRatio ) * (180/pi);
% 
%     % 5) Finally draw the text
%     hText = text( ax, x, y, txt, ...
%                   'Rotation',           angDeg, ...
%                   varargin{:} );
% end


% % Place TEXT at (x,y) on log‑log axes AX, rotated to slope m = d(log y)/d(log x).
% % Works inside a tiledlayout or plain figure.
% 
%   % 1) Force the figure to finish sizing the tile/axes:
%   drawnow
% 
%   % 2) Measure the axes in pixels:
%   oldU = ax.Units;
%   ax.Units = 'pixels';
%   apos = ax.Position;        % [x y width height] in px
%   ax.Units = oldU;
% 
%   % 3) Make sure it’s valid:
%   if apos(3) <= 0 || apos(4) <= 0
%     error('Axes size is zero.  Call drawnow after plotting before using this function.');
%   end
% 
%   % 4) Compute rotation in degrees:
%   angDeg = atan( m * (apos(4)/apos(3)) ) * (180/pi);
% 
%   % 5) Draw the text:
%   hText = text(ax, x, y, txt, ...
%                'Rotation',       angDeg, ...
%                varargin{:});
% end


% B2KTEXTALIGNEDTOLINELOGPLOTPLOT
%   Place TEXT at (x,y) in data units on the log‑log axes AX,
%   rotated so that it follows the slope
%       m = d(log10 y)/d(log10 x).
%   Works inside a tiledlayout or a plain figure.
%
%   H = B2KTextAlignedToLineLogLogPlot(...)
%   returns the handle to the text object.
%
%   Additional Name/Value pairs (e.g. 'Interpreter','latex') go in varargin.

  % 1) Make sure the tile/axes have been drawn & sized
  drawnow;

  % 2) Measure the **inner** plot‐box in pixels
  oldUnits = ax.Units;
  ax.Units  = 'pixels';

    pos    = ax.Position;     % [left bottom totalWidth totalHeight]
    inset  = ax.TightInset;   % [left bottom right top] margins in px

  ax.Units = oldUnits;

  % Compute the width & height of the data region
  dataW = pos(3) - inset(1) - inset(3);
  dataH = pos(4) - inset(2) - inset(4);

  if dataW <= 0 || dataH <= 0
    error('Couldn''t compute positive data‐region size.  Check your axes layout.');
  end

  % 3) Compute the on‐screen angle
  %    slope in pixel‐space = m * (dataH/dataW)
  angDeg = atan( m * (dataH/dataW) ) * (180/pi);

  % 4) Place the text, rotated by that angle
  hText = text(ax, x, y, txt, ...
               'Rotation',          angDeg, ...
               varargin{:} );
end
