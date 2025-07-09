function mn = getgracemean(gsm, s, e)

% GETGRACEMEAN computes the mean of the dataset given a start and end
% date.
%
% mn = getgracemean(gsm)
% mn = getgracemean(gsm, s, e)
%
% INPUT
% gsm   - The 10-column internal format. See READGRC
% s     - Start date [year month]
% e     - End date   [year month]
%
% OUTPUT
% mn    - Mean of the GRACE monthly values
%
%-----------------------------------------------------------------------

% Created on: 4 August 2016, Hannover
% Author: Balaji Devaraju
%-----------------------------------------------------------------------

narginchk(1, 3)

if iscell(gsm)
    if ~ismember(10, size(gsm, 2))
        error('Unexpected format. Expecting a 10-column cell-array.')
    end
else
    error(sprintf(['Unexpected variable type: %s\n', ...
                    'Expecting a 10-column cell-array.'], class(gsm)))
end

if exist('s', 'var') == 1
    if numel(s) ~= 2
        error(sprintf(['Unexpected number of elements in the start date: %g\n', ...
                        'Expecting a [2 x 1] or [1 x 2] vector as input'], numel(s)))
    end
else
    s = [gsm{1,4} gsm{1,5}];
end

if exist('e', 'var') == 1
    if numel(e) ~= 2
        error(sprintf(['Unexpected number of elements in the start date: %g\n', ...
                        'Expecting a [2 x 1] or [1 x 2] vector as input'], numel(s)))
    end
else
    e = [gsm{end,4} gsm{end,5}];
end


t   = cell2mat(gsm(:,4:6));

k1  = find((t(:,1) == s(1)) & (t(:,2) == s(2)));
if isempty(k1)
    error(sprintf('Start date %g/%g does not exist', s))
end

k2  = find((t(:,1) == e(1)) & (t(:,2) == e(2)));
if isempty(k2)
    error(sprintf('End date %g/%g does not exist', e))
end

mn  = zeros(size(gsm{1,9}));
for k = k1:k2
    mn = mn + gsm{k,9}/(k2-k1+1);
end
mn(1:2,1:2)     = 0;
mn(1,1)         = 1;
