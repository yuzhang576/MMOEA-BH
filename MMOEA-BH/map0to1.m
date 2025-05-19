function varargout = map0to1(x, varargin)
    if nargin == 1
        xmin = min(x, [], 1);
        xmax = max(x, [], 1);
    elseif nargin == 3
        xmin = varargin{1};
        xmax = varargin{2};
    else
        error('input error')
    end
    x1 = (x-xmin)./(xmax-xmin);
    if any(isnan(x1),'all')
        x1 = x;
    end
    varargout = {x1, xmin, xmax};
end