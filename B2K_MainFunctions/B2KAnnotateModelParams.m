function hHandles = B2KAnnotateModelParams(ax, syms, beta, SE, x0, y0, deltaSpacing, varargin)
%% B2KAnnotateModelParams.m - [Function] Outputs text annotation for model parameters with uncertainties to 1 significant figure
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
% B2KANNOTATEMODELPARAMS  Stacked, aligned annotations for model parameters.
%
% H = B2KAnnotateModelParams(AX,SYMS,BETA,SE,X0,Y0,DELTASPACING,...)
%   Plots “symbol = value ± se” at (X0,Y0) on axes AX, spaced by DELTASPACING.
%
% Name–value options:
%   'Interpreter'        (default 'latex')
%   'FontSize'           (default 12)
%   'VerticalAlignment'  (default 'bottom')
%   'AlignEqual'         align all “=” columns         (default false)
%   'AlignPM'            align all “±” columns         (default false)
%   'YSpacing'           'log'|'linear'|numeric        (default auto)
%   'AlignDigit'         align integer‐digit columns   (default false)
%   'AlignDecimalPoint'  align decimal/fraction cols   (default false)
%   'AlignPosNeg'        phantom-pad “+” for positives (default false)

  %% 1) Parse inputs
  p = inputParser;
  p.KeepUnmatched   = true;
  p.PartialMatching = false;
  addRequired(p,'ax',           @(h) ishghandle(h,'axes'));
  addRequired(p,'syms',         @(c) iscellstr(c)||isstring(c));
  addRequired(p,'beta',         @isnumeric);
  addRequired(p,'SE',           @isnumeric);
  addRequired(p,'x0',           @isnumeric);
  addRequired(p,'y0',           @isnumeric);
  addRequired(p,'deltaSpacing', @isnumeric);
  addParameter(p,'Interpreter','latex',  @ischar);
  addParameter(p,'FontSize',   12,       @isnumeric);
  addParameter(p,'VerticalAlignment','bottom',@ischar);
  addParameter(p,'AlignEqual', false,    @islogical);
  addParameter(p,'AlignPM',    false,    @islogical);
  addParameter(p,'YSpacing',   '',       @(s) ischar(s)||isnumeric(s));
  addParameter(p,'AlignDigit',        false,@islogical);
  addParameter(p,'AlignDecimalPoint', false,@islogical);
  addParameter(p,'AlignPosNeg',       false,@islogical);
  parse(p,ax,syms,beta,SE,x0,y0,deltaSpacing,varargin{:});
  opts = p.Results;

  fn = fieldnames(p.Unmatched);
  fv = struct2cell(p.Unmatched);
  extraArgs = reshape([fn.'; fv.'],1,[]);

  %% 2) Compute vertical positions
  n = numel(opts.syms);
  assert(numel(opts.beta)>=n && numel(opts.SE)>=n, 'Lengths must be ≥ %d', n);

  if isempty(opts.YSpacing)
    if strcmpi(get(ax,'YScale'),'log')
      spacingMode = 'log';
    else
      spacingMode = 'linear';
    end
  else
    spacingMode = opts.YSpacing;
    if isstring(spacingMode), spacingMode = char(spacingMode); end
  end

  if ischar(spacingMode)
    m = lower(spacingMode);
    if strcmp(m,'log')
      yVals = opts.y0 * 10 .^ (-(0:n-1)*opts.deltaSpacing);
    elseif strcmp(m,'linear')
      yVals = opts.y0 - (0:n-1)*opts.deltaSpacing;
    else
      error('YSpacing must be ''log'', ''linear'', or numeric vector.');
    end
  else
    assert(numel(spacingMode)>=n,'YSpacing vector must have ≥ %d elements',n);
    yVals = spacingMode(1:n);
  end

  wasHold = ishold(ax);
  hold(ax,'on');

  %% 3) SIMPLE MODE
  if ~opts.AlignEqual
    % choose inline ± token
    if ~opts.AlignPM && ~opts.AlignDigit && ...
       ~opts.AlignDecimalPoint && ~opts.AlignPosNeg
      pmToken = '\,\pm\,';
    else
      pmToken = '\,\pm\!\!\!\!';
    end

    hHandles = gobjects(n,1);
    for i = 1:n
      [seR, dp] = oneSigFigAndDP(opts.SE(i));
      vR = round(opts.beta(i), dp);

      % format value, force "0" if zero
      if vR == 0
        rv = '0';
      else
        if dp >= 0
          fmt = ['%.' num2str(dp) 'f'];
        else
          fmt = '%.0f';
        end
        rv = sprintf(fmt, vR);
        if dp == 0
          rv = stripZeros(rv);
        end
      end

      % format SE, force "0" if zero
      if seR == 0
        re = '0';
      else
        if dp >= 0
          fmt = ['%.' num2str(dp) 'f'];
        else
          fmt = '%.0f';
        end
        re = sprintf(fmt, seR);
        if dp == 0
          re = stripZeros(re);
        end
      end

      % phantom-pad sign if requested
      if opts.AlignPosNeg
        if ~startsWith(rv,'-')
          rv = ['\phantom{-}', rv];
        end
        if ~startsWith(re,'-')
          re = ['\phantom{-}', re];
        end
      end

      txt = ['$', opts.syms{i}, '=\,', rv, pmToken, re, '$'];
      hHandles(i) = text(ax, x0, yVals(i), txt, ...
                         'Interpreter',opts.Interpreter, ...
                         'FontSize',   opts.FontSize, ...
                         'HorizontalAlignment','center', ...
                         'VerticalAlignment',  opts.VerticalAlignment, ...
                         extraArgs{:});
    end

    if ~wasHold, hold(ax,'off'); end
    return
  end

  %% 4) ALIGNED MODE: SYMBOL | = VALUE | ± ERROR

  % a) SYMBOL column
  hSym = gobjects(n,1);
  for i = 1:n
    hSym(i) = text(ax,x0,yVals(i),sprintf('$%s$',opts.syms{i}), ...
                   'Interpreter',opts.Interpreter, ...
                   'FontSize',   opts.FontSize, ...
                   'HorizontalAlignment','right', ...
                   'VerticalAlignment',opts.VerticalAlignment, ...
                   extraArgs{:});
  end
  drawnow;
  origU   = get(ax,'Units'); set(ax,'Units','pixels');
  ext     = cell2mat(get(hSym,'Extent'));
  maxSymW = max(ext(:,3));
  set(ax,'Units',origU);

  % b) Precompute raw strings & absolute lengths
  rawV  = cell(n,1); rawE = cell(n,1);
  intV  = zeros(n,1); fracV = zeros(n,1);
  intE  = zeros(n,1); fracE = zeros(n,1);
  for i = 1:n
    [seR, dp] = oneSigFigAndDP(opts.SE(i));
    vR        = round(opts.beta(i), dp);

    % value
    if vR == 0
      sv = '0';
    else
      if dp >= 0
        fmt = ['%.' num2str(dp) 'f'];
      else
        fmt = '%.0f';
      end
      sv = sprintf(fmt, vR);
      if dp == 0
        sv = stripZeros(sv);
      end
    end

    % SE
    if seR == 0
      seStr = '0';
    else
      if dp >= 0
        fmt = ['%.' num2str(dp) 'f'];
      else
        fmt = '%.0f';
      end
      seStr = sprintf(fmt, seR);
      if dp == 0
        seStr = stripZeros(seStr);
      end
    end

    rawV{i} = sv;
    rawE{i} = seStr;

    absV = regexprep(sv,'^-','');
    absE = regexprep(seStr,'^-','');
    % integer part
    tokensV = regexp(absV,'^(\d+)','tokens','once');
    if isempty(tokensV)
      intV(i) = 0;
    else
      intV(i) = strlength(tokensV{1});
    end
    tokensE = regexp(absE,'^(\d+)','tokens','once');
    if isempty(tokensE)
      intE(i) = 0;
    else
      intE(i) = strlength(tokensE{1});
    end
    % fraction part
    fracTokensV = regexp(absV,'\.(\d+)$','tokens','once');
    if isempty(fracTokensV)
      fracV(i) = 0;
    else
      fracV(i) = strlength(fracTokensV{1});
    end
    fracTokensE = regexp(absE,'\.(\d+)$','tokens','once');
    if isempty(fracTokensE)
      fracE(i) = 0;
    else
      fracE(i) = strlength(fracTokensE{1});
    end
  end
  maxIntV  = max(intV);   maxFracV = max(fracV);
  maxIntE  = max(intE);   maxFracE = max(fracE);

  % c) "= VALUE" column
  xVal0 = x0 + maxSymW;
  hVal  = gobjects(n,1);
  for i = 1:n
    sv = rawV{i};

    % sign phantom
    if opts.AlignPosNeg
      if startsWith(sv,'-')
        sp = '-'; mag = sv(2:end);
      else
        sp = '\phantom{-}'; mag = sv;
      end
    else
      if startsWith(sv,'-')
        sp = '-'; mag = sv(2:end);
      else
        sp = '';  mag = sv;
      end
    end

    % integer pad
    if opts.AlignDigit
      pd = maxIntV - intV(i);
      if pd>0
        padI = repmat('\phantom{0}',1,pd);
      else
        padI = '';
      end
    else
      padI = '';
    end

    % fraction pad
    if opts.AlignDecimalPoint
      if fracV(i)==0
        padF = ['\phantom{.' repmat('0',1,maxFracV) '}'];
      elseif fracV(i)<maxFracV
        padF = repmat('\phantom{0}',1,maxFracV-fracV(i));
      else
        padF = '';
      end
    else
      padF = '';
    end

    alignedV = [sp, padI, mag, padF];
    txtV     = ['$=\,', alignedV, '$'];
    hVal(i)  = text(ax, xVal0, yVals(i), txtV, ...
                   'Interpreter',opts.Interpreter, ...
                   'FontSize',   opts.FontSize, ...
                   'HorizontalAlignment','left', ...
                   'VerticalAlignment',opts.VerticalAlignment, ...
                   extraArgs{:});
  end

  % d) "± ERROR"
  if ~opts.AlignPM
    for i = 1:n
      seStr = rawE{i};

      % sign phantom
      if opts.AlignPosNeg
        if startsWith(seStr,'-')
          spE = '-'; mgE = seStr(2:end);
        else
          spE = '\phantom{-}'; mgE = seStr;
        end
      else
        if startsWith(seStr,'-')
          spE = '-'; mgE = seStr(2:end);
        else
          spE = '';  mgE = seStr;
        end
      end

      % integer pad
      if opts.AlignDigit
        pdE = maxIntE - intE(i);
        if pdE>0, padE = repmat('\phantom{0}',1,pdE);
        else      padE = '';
        end
      else
        padE = '';
      end

      % fraction pad
      if opts.AlignDecimalPoint
        if fracE(i)==0
          padFE = ['\phantom{.' repmat('0',1,maxFracE) '}'];
        elseif fracE(i)<maxFracE
          padFE = repmat('\phantom{0}',1,maxFracE-fracE(i));
        else
          padFE = '';
        end
      else
        padFE = '';
      end

      alignedE = [spE, padE, mgE, padFE];
      old = get(hVal(i),'String');
      if old(end)=='$', old=old(1:end-1); end
      newStr = [old, '\,\pm\,', alignedE, '$'];
      set(hVal(i),'String',newStr);
    end
    hHandles = hVal;

  else
    % separate ± column
    drawnow;
    origU = get(ax,'Units'); set(ax,'Units','pixels');
    ext2  = cell2mat(get(hVal,'Extent'));
    maxValW = max(ext2(:,3));
    set(ax,'Units',origU);

    xErr0 = xVal0 + maxValW;
    hErr  = gobjects(n,1);
    for i = 1:n
      seStr = rawE{i};

      % sign phantom
      if opts.AlignPosNeg
        if startsWith(seStr,'-')
          spE = '-'; mgE = seStr(2:end);
        else
          spE = '\phantom{-}'; mgE = seStr;
        end
      else
        if startsWith(seStr,'-')
          spE = '-'; mgE = seStr(2:end);
        else
          spE = '';  mgE = seStr;
        end
      end

      % integer pad
      if opts.AlignDigit
        pdE = maxIntE - intE(i);
        if pdE>0, padE = repmat('\phantom{0}',1,pdE);
        else      padE = '';
        end
      else
        padE = '';
      end

      % fraction pad
      if opts.AlignDecimalPoint
        if fracE(i)==0
          padFE = ['\phantom{.' repmat('0',1,maxFracE) '}'];
        elseif fracE(i)<maxFracE
          padFE = repmat('\phantom{0}',1,maxFracE-fracE(i));
        else
          padFE = '';
        end
      else
        padFE = '';
      end

      alignedE = [spE, padE, mgE, padFE];
      txtE     = ['$', '\,\pm\,', alignedE, '$'];
      hErr(i)  = text(ax, xErr0, yVals(i), txtE, ...
                      'Interpreter',opts.Interpreter, ...
                      'FontSize',   opts.FontSize, ...
                      'HorizontalAlignment','left', ...
                      'VerticalAlignment',opts.VerticalAlignment, ...
                      extraArgs{:});
    end
    hHandles = hErr;
  end

  %% 5) Restore hold
  if ~wasHold, hold(ax,'off'); end
end

%% Helpers
function [rSE, dp] = oneSigFigAndDP(x)
  if x == 0
    rSE = 0; dp = 0;
  else
    e   = floor(log10(abs(x)));
    rSE = round(x / 10^e) * 10^e;
    dp  = -e;
  end
end

function s = stripZeros(s)
  s = regexprep(s,'(\.\d*?[1-9])0+$','$1');
  s = regexprep(s,'\.0+$','');
  s = regexprep(s,'^-0$','0');
end