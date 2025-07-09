function cm = catchagg(fld,cindx,indx)

% CATCHAGG aggregates storage changes of the catchments into area-weighted
% averages, either for the specified catchments or all the catchments.
%
% cm = catchagg(fld,cindx)
% cm = catchagg(fld,cindx,indx)
%
% INPUT
% fld   -   Field of values whose values have to be aggregated. If multiple
%           fields are to be provided, they can provided in a cell array.
% cindx -   Map in which the individual pixels are assigned to the
%           respective catchments.
% indx  -   Index of catchments whose aggregated value is desired.
% Note: 'fld' and 'cindx' must have the same longitudinal ordering either
% [0:dlam:360] or [-180:dlam:180]
%
% OUTPUT 
% cm    -   Catchment aggregated values, with the first column bearing the
%           catchment indices.
%--------------------------------------------------------------------------

% Balaji Devaraju. 6 August 2012, Stuttgart.

if nargin < 2
    error('Insufficient input arguments')
elseif nargin == 2
    if iscell(fld)
        if ~isequal(size(fld{1}),size(cindx))
            error('Dimension mismatch between the field and catchment index matrices')
        end
        cll = true;
        l   = size(fld(:),1);
    else
        if ~isequal(size(fld),size(cindx))
            error('Dimension mismatch between the field and catchment index matrices')
        end
        cll = false;
    end
    indx = unique(cindx(:));
    indx(indx==-9999) = [];
elseif nargin == 3
    if iscell(fld)
        if ~isequal(size(fld{1}),size(cindx))
            error('Dimension mismatch between the field and catchment index matrices')
        end
        cll = true;
        l   = size(fld(:),1);
    else
        if ~isequal(size(fld{1}),size(cindx))
            error('Dimension mismatch between the field and catchment index matrices')
        end
        cll = false;
    end
    if min(size(indx)) ~= 1
        error('The index values must be given as a scalar or a vector')
    end
    idx = unique(cindx(:));
    idx = idx(idx>-1);
    t   = zeros(size(indx));
    for k = 1:length(indx)
        t(k) = ~isempty(find(indx(k) == idx));
    end
    t = logical(t);
    if ~isempty(indx(~t))
        fprintf('Warning: Some invalid catchment IDs were detected: \n %g \n',indx(~t)')
    end
    indx = indx(t);
end

[m,n] = size(cindx);
dlam  = pi/m;
dt    = dlam;
th    = (0:dt:pi)';

constants_grace % Loading the constants specific to GRACE.
                % We will only need the semi-major axis of
                % reference ellipsoid.

Area = ae^2*dlam*dt*(cos(th(1:end-1)) - cos(th(2:end)))*ones(1,n);
fA   = fld.*Area;

if cll
    cm = [indx(:) zeros(length(indx),l)];
    for k = 1:length(indx)
        t = (indx(k)==cindx);
        for q = 1:l
            cm(k,q) = sum(fA{q}(t))/sum(Area(t));
        end
    end
else
    cm = [indx(:) zeros(length(indx),1)];
    for k = 1:length(indx)
        t = (indx(k)==cindx);
        cm(k,2) = sum(fA(t))/sum(Area(t));
    end
end
