function klm = remgracemonths(gsm, varargin)

% REMGRACEMONTHS removes troublesome months in the GRACE time-series.
% The troublesome months are those that have been
%   a. regularized
%   b. the data used for computing the monthly solution comes from two
%      adjacent months
%
% INPUT
% gsm   - 10-column cell-array (internal) format.
% OPTIONAL VALUES
% 'regularized' - Removes regularized months (def: false)    [logical]
% 'spans2'      - Removes months whose solution spans two months
%                                            (def: false)    [logical]
%
% OUTPUT
% klm   - Time-series with the troublesome months removed given out as
%         10-column cell-array (internal) format.
%
%----------------------------------------------------------------------


% Checking the number of inputs
narginchk(1,Inf)

% Checking the format of inputs
if iscell(gsm)
    if size(gsm, 2) ~= 10
        error(sprintf(['Unexpected dimensions of the input GRACE data: [%g x %g]\n', ...
                        'Expecting a [n x 10] array'], size(gsm)))
    end
else
    error(sprintf(['Unexpected variable type of the input GRACE data: %s\n', ...
                    'Expecting a cell-array'], class(gsm)))
end

% Default values of optional input variables
defpar = {  'regularized',  false; ...
            'spans2',       false};
params = getopt(defpar, false, varargin);

% Removing regularized months
if params.regularized
    v   = cell2mat(gsm(:,7));
    gsm(v==2,:) = [];
end

% Removing entries that spans two months
if params.spans2
    t   = cell2mat(gsm(:,4:6));
    idx = (t(:,3)>t(:,4));
    gsm(idx,:)  = [];
    
    t   = cell2mat(gsm(:,4:6));
    m1  = doy2cal(t(:,1), t(:,3));
    m2  = doy2cal(t(:,1), t(:,4));
    gsm(m1~=m2,:)  = [];
end

klm = gsm;
