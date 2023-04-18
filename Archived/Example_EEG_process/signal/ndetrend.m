%
% y = ndetrend(x,dim,order,breaks)
%
% NDETREND removes a piecewise polynomial trend from data array x along 
% dimension dim.  The trend will be continuous between pieces unless the
% order of the trend is 0.
%
%       x: array to be detrended
%     dim: (optional) dimension to detrend; default is first nonsingleton
%   order: (optional) order of the piecewise polynomial; default is 1
%  breaks: (optional) breakpoint indices for each piece of polynomial
%
%       y: detrended array
%
% See also POLYTREND, REREF, DETREND, MEAN
%
% Created by Bill Winter, May 2006: permute and matrix multiplication
% based on DETREND.M, Copyright 1984-2004 The MathWorks, Inc.
% Modified by Bill Winter November 2012: reshape, matrix multiplication,
%    and power

% rather than slash operator, use qr and slash?
function x = ndetrend(x,dim,o,b)
siz = size(x);                                      % array size
if nargin<2 || isempty(dim), dim=find(siz>1,1); end % default: nonsingleton
if isempty(dim), x = zeros(siz); return; end        % perfect fit

if nargin < 3, o = 1; end
o = max(0,floor(o));                                % positive integers

N = siz(dim);                                       % dimension length
if nargin < 4, b = [1 N+1];l = 1;                   % default: one piece
else
    b = unique([1; b(:); N+1]);                     % breaks unique
    b(b < 1 | b > N+1) = [];                        % breaks within array
    l = length(b) - 1;                              % number of pieces
end

if o == 1                                           % piecewise linear
    M = N - b;                                      % length of segments
    a = [zeros(N,l) ones(N,1)];                     % preallocate linear
    for k = 1:l, a(1+b(k):N,k) = (1:M(k))/M(k); end
elseif o                                            % piecewise poly
    M = N - b;                                      % length of segments
    a = [zeros(N,l*o) ones(N,1)];                   % preallocate poly
    for k = 1:l
        a(1+b(k):N,o*(k-1)+1:o*k) = bsxfun(@power,(1:M(k)).'/M(k),1:o);
    end
elseif l > 1                                        % piecewise constant
    M = diff(b);                                    % length of segments
    for k = l:-1:1, a(b(k):b(k+1)-1,k) = 1/M(k); end% constant pieces
else                                                % single constant
    x = bsxfun(@minus,x,mean(x,dim));
    return;
end

N = prod(siz(dim+1:end));
if dim == 1
    try                                             % high-memory
        x(:,:) = x(:,:) - a*(a\x(:,:));
    catch                                           %  low-memory
        for k = 1:N, x(:,k) = x(:,k) - a*(a\x(:,k)); end
        % if second dimension isn't too long, could try:
        % for k = 1:prod(siz(3:end)), x(:,:,k)=x(:,:,k)-a*(a\x(:,:,k)); end
    end
else
    a = a';
    if dim ~= 2, x = reshape(x,[prod(siz(1:dim-1)) siz(dim) N]); end
    for k = 1:N, x(:,:,k) = x(:,:,k) - (x(:,:,k)/a)*a; end
    % consider indexing first dimension and reshaping the rest if small
    if dim ~= 2, x = reshape(x,siz); end
end

% assumes we only want unique integers 
function x = unique(x)
x = sort(floor(x));
x = x([logical(diff(x));true]);