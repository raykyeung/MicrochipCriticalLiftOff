function stringCell = B2KFormatValuesToStrOneSigFigMaxUnc(values, uncertainties, varargin)
%% B2KFormatValuesToStrOneSigFigMaxUnc.m - [Function] Rounds values based on 1 significant figure of the maximum uncertainty for all values and formats as strings
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
% Rounds entries in VALUES according to one significant figure of the
% maximum entry in UNCERTAINTIES, unless overridden by 'Precision'.
% Uses “round half up” (ties → away from zero).
%
%   stringCell = B2KFormatValuesToStrOneSigFigMaxUnc(values, uncertainties)
%   stringCell = B2KFormatValuesToStrOneSigFigMaxUnc(..., 'Precision', P)
%
% Inputs:
%   values        numeric array
%   uncertainties numeric array (same size as values)
%
% Name-Value:
%   'Precision'   either [nBefore,nAfter] or 'nBefore.nAfter'
%
% Behavior (no Precision given):
%   • Let maxUnc = max(abs(uncertainties(:)))
%   • exponent      = floor(log10(maxUnc))
%   • decimalPlaces = –exponent
%     – if exponent<0 → decimalPlaces>0 → round to that many decimals
%     – if exponent=0 → decimalPlaces=0 → round to integer
%     – if exponent>0 → decimalPlaces<0 → round to tens/hundreds/etc.
%
% Output formatting always displays max(decimalPlaces,0) digits
% after the decimal point.

    %— parse optional
    p = inputParser;
    addParameter(p, 'Precision', [], @(x) ...
        (isnumeric(x) && numel(x)==2) || ischar(x) || isstring(x));
    parse(p, varargin{:});
    prec = p.Results.Precision;

    %— validate required inputs
    assert(isnumeric(values) && isnumeric(uncertainties) && ...
        numel(values)==numel(uncertainties), ...
        'values and uncertainties must be numeric arrays of the same size.');

    if ~isempty(prec)
        %— user override
        if ischar(prec) || isstring(prec)
            parts = split(string(prec),'.');
            assert(numel(parts)==2, 'Precision must be "nBefore.nAfter".');
            nBefore = str2double(parts(1));
            nAfter  = str2double(parts(2));
        else
            nBefore = prec(1);
            nAfter  = prec(2);
        end

        if nAfter > 0
            decimalPlaces = nAfter;
        else
            % e.g. [2 0] → 2 digits before decimal → nearest 10^(2-1)
            decimalPlaces = -max(0, nBefore-1);
        end

    else
        %— one-sig-fig logic on the max uncertainty
        maxUnc = max(abs(uncertainties(:)));
        assert(maxUnc > 0, 'Maximum uncertainty must be positive.');
        exponent      = floor(log10(maxUnc));
        decimalPlaces = -exponent;
    end

    %— build format (never negative in the fmt)
    dp  = max(decimalPlaces, 0);
    fmt = sprintf('%%.%df', dp);

    %— round half-up & format
    stringCell = arrayfun(@(x) ...
        num2str( roundHalfUp(x, decimalPlaces), fmt ), ...
        values, 'UniformOutput', false);
end

function y = roundHalfUp(x, N)
% Round x to N digits after (N>=0) or before (N<0) the decimal, half up.
    if N >= 0
        f = 10^N;
        y = sign(x) .* floor(abs(x).*f + 0.5) ./ f;
    else
        f = 10^(-N);
        y = sign(x) .* floor(abs(x)./f + 0.5) .* f;
    end
end
