function cmap = fillcatch(val,cindx)

% FILLCATCH fills the given ctachment-wise data in all the pixels of the
% corresponding catchment. 
%
% cmap = fillcatch(val,cindx)
%
% INPUT
% val   - a two column vector having the ctahcment index in the first column
%         and the value to be filled in the second column.
% cindx - map of the catchments containing the catchment index in the
%         corresponding pixels.
%
% OUTPUT
% cmap  - catchment map filled with the values provided.
%
%--------------------------------------------------------------------------

% Balaji Devaraju. 6 September 2012, Stuttgart

if nargin < 2
    error('Insufficient input arguments')
end

cmap = ones(size(cindx))*-9999;

for k = 1:size(val,1)
    idx = (cindx==val(k,1));
    cmap(idx) = val(k,2);
end
