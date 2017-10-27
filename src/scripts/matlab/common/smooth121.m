function [c,ww] = smooth121(varargin)
%SMOOTH121  Smooth data using 1-2-1 type filter
%   Z = SMOOTH121(Y) smooths data Y using a 1-2-1 filter.
%
%   Z = SMOOTH121(Y,SPAN) smooths data Y using SPAN as the number of points used
%   to compute each element of Z.
%
%   Z = SMOOTH121(X,Y,SPAN) additionally specifies the X coordinates.  If X is
%   not provided, methods that require X coordinates assume X = 1:N, where
%   N is the length of Y.
%
%   Notes:
%   1. When X is given and X is not uniformly distributed, a warning will
%   be printed; the 1-2-1 type filters are not recommended.
%
%   2. SPAN must be odd so the filter is symmetric. If an even SPAN is
%   specified, it is reduced by 1.
%
%   3. If SPAN is greater than the length of Y, it is reduced to the
%   length of Y.
%
%   For example:
%
%   Z = SMOOTH121(Y) uses the 1-2-1 filter and X=1:length(Y).
%
%   Z = SMOOTH121(Y,5) uses the 1-2-3-2-1 filter and X=1:length(Y).
%
%   Z = SMOOTH121(X,Y) where X is unevenly distributed uses the 1-2-1 filter,
%   but a warning is printed first.
%
%
%   Leah Grant (2016) copied and modified this code from the smooth function


%   Copyright 2001-2012 The MathWorks, Inc.

if nargin < 1
    error('SMOOTH121 needs at least one argument');
end

if nargout > 1 % Called from the GUI cftool
    ws = warning('off', 'all'); % turn warning off and record the previous warning state.
    [lw,lwid] = lastwarn;
    lastwarn('');
else
    ws = warning('query','all'); % Leave warning state alone but save it so resets are no-ops.
end

% is x given as the first argument?
if nargin==1 || ( nargin > 1 && (length(varargin{2})==1 || ischar(varargin{2})) )
    % smooth121(Y) | smooth121(Y,span)
    is_x = 0; % x is not given
    y = varargin{1};
    y = y(:);
    x = (1:length(y))';
else % smooth121(X,Y,...)
    is_x = 1;
    y = varargin{2};
    x = varargin{1};
    y = y(:);
    x = x(:);
end

% is span given?
span = [];
if nargin == 1+is_x
    % smooth121(Y), smooth121(X,Y)
    is_span = 0;
else
    % smooth121(...,SPAN)
    is_span = 1;
    span = varargin{2+is_x};
end

t = length(y);
if t == 0
    c = y;
    ww = '';
    if nargout > 1
        ww = lastwarn;
        lastwarn(lw,lwid);
        warning(ws);  % turn warning back to the previous state.
    end
    return
elseif length(x) ~= t
    warning(ws); % reset warn state before erroring
    error('X and Y must be the same length');
end

% check if x is uniformly distributed
diffx = diff(x);
if ~uniformx(diffx,x,y)
    disp('Warning: X is not uniform; this function is not recommended! See smooth')
end

% realize span
if span < 1
    warning(ws); % reset warn state before erroring
    error('SPAN must be >= 1');
end
if isempty(span), span = 3; end % smooth121(Y) | smooth121(X,Y)

idx = 1:t;

sortx = any(diff(isnan(x))<0);   % if NaNs not all at end
if sortx || any(diff(x)<0) % sort x
    [x,idx] = sort(x);
    y = y(idx);
end

c = NaN(size(y),'like',y);
ok = ~isnan(x);
c(ok) = moving121(x(ok),y(ok),span);

c(idx) = c;

if nargout > 1
    ww = lastwarn;
    lastwarn(lw,lwid);
    warning(ws);  % turn warning back to the previous state.
end
%--------------------------------------------------------------------
function c = moving121(x,y, span)
% moving average of the data.

ynan = isnan(y);
span = floor(span);
n = length(y);
span = min(span,n);
width = span-1+mod(span,2); % force it to be odd
xreps = any(diff(x)==0);
if width==1 && ~xreps && ~any(ynan), c = y; return; end
if ~xreps && ~any(ynan)
    % simplest method for most common case
    % number of endpoints
    nendpts = floor(width/2);
    % initialize variables so they have the same non-singleton dimension
    weights=zeros(width,1); [cbegin,cend]=deal(zeros(nendpts,1));
    % set up symmetric weights
    weights(1:ceil(width/2)) = 1:ceil(width/2);
    weights(ceil(width/2)+1:width) = fliplr(weights(1:floor(width/2)));
    % get filtered y using weights (where the sum of weights adds to 1)
    c = filter(weights/sum(weights),1,y);
    % filter outputs data averaged from previous values, not centered. e.g.
    % if filter is 1-2-3-2-1, c(n) will be averaged from values at
    % y(n), y(n-1), ... y(n-4). So for this example, c needs to be shifted 
    % left by 2 so c(n) is averaged from y(n-2)...y(n+2).
    % Also, endpoints n=1, n=2, n=end-1, n=end need to be re-calculated 
    % from weights (where weights always sum to 1)
    for i=1:nendpts
        cbegin(i) = sum( y(1:nendpts+i) .* ...
                         weights(width-nendpts-i+1:width)/sum(weights(width-nendpts-i+1:width)) );
        cend(i) = sum( y(n-width+1+i:n) .* ...
                         weights(1:width-i)/sum(weights(1:width-i)) );
    end
    c = [cbegin;c(width:end);cend];
elseif ~xreps % I am not sure what this does
    % with no x repeats, can take ratio of two smoothed sequences
    yy = y;
    yy(ynan) = 0;
    nn = double(~ynan);
    ynum = moving121(x,yy,span);
    yden = moving121(x,nn,span);
    c = ynum ./ yden;
else % I am not sure what this does either
    % with some x repeats, loop
    notnan = ~ynan;
    yy = y;
    yy(ynan) = 0;
    c = zeros(n,1,'like',y);
    for i=1:n
        if i>1 && x(i)==x(i-1)
            c(i) = c(i-1);
            continue;
        end
        R = i;                                 % find rightmost value with same x
        while(R<n && x(R+1)==x(R))
            R = R+1;
        end
        hf = ceil(max(0,(span - (R-i+1))/2));  % need this many more on each side
        hf = min(min(hf,(i-1)), (n-R));
        L = i-hf;                              % find leftmost point needed
        while(L>1 && x(L)==x(L-1))
            L = L-1;
        end
        R = R+hf;                              % find rightmost point needed
        while(R<n && x(R)==x(R+1))
            R = R+1;
        end
        c(i) = sum(yy(L:R)) / sum(notnan(L:R));
    end
end
%--------------------------------------------------------------------
function isuniform = uniformx(diffx,x,y)
%ISUNIFORM True if x is of the form a:b:c

if any(isnan(y)) || any(isnan(x))
    isuniform = false;
else
    isuniform = all(abs(diff(diffx)) <= eps*max(abs([x(1),x(end)])));
end
