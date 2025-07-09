function otpt = findindx(inpt,cindx,sid)

% FINDINDX retrieves the time-series of particular catchments from the output 
% of CATCHTSRS or similar matrices.
%
% otpt = findindx(inpt,cindx)
% otpt = findindx(inpt,cindx,sid)
%
% INPUT
% inpt  - Time-series matrix [time x data]. The first row of the matrix must 
%         contain the indices of the catchments.
% cindx - Indices of catchments whose time-series have to be retrieved.
% sid   - The function assumes that the time-series starts from the third 
%         column. If this is not the case the column number must be provided.
%
% OUTPUT
% otpt  - Retrieved time-series of the catchments.
%------------------------------------------------------------------------------

% Christof Lorenz. Garmisch-Partenkirchen
%------------------------------------------------------------------------------

if nargin < 2
    error('Insufficient input arguments')
elseif nargin == 2
    sid = 3;
end


otpt(:, 1:sid-1) = inpt(:, 1:sid-1);

for k = 1:length(cindx)
    indx = (inpt(1, :) == cindx(k));
    otpt(:, k+sid) = inpt(:, indx);
end
