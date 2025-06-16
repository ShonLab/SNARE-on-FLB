function h = imshow3(varargin)
img = varargin{1};
if nargin >= 3
    n = varargin{2};
    varargin = [img, varargin(3:end)];
elseif nargin == 2
    n = varargin{2};
    varargin = {img};
end

img = single(img);
low = median(img(:))-median(n*std(img));
high = median(img(:))+median(n*std(img));

if low == high
    high = high+1; low = low-1;
end

varargin = [varargin,'DisplayRange',[low,high]];
hh = imshow2(varargin{:});
if (nargout > 0)
% Only return handle if caller requested it.
h = hh;
end