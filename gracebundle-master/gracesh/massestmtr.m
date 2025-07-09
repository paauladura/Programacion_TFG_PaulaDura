function [mest,t,l] = massestmtr(gsm,syn,mn,ds,cindx,indxnr)

% MASSESTMTR estimates the mass changes for every month in either surface
% mass density (kg/m^2) or water-equivalent-height (m)
%
% [mest,t,l] = massestmtr(gsm,syn,mntype,ds,cindx,indxnr)
%
% INPUT
% gsm       - SH co-efficient cells for every release
% syn       - Structure variable containing inputs for GSHS_RAD function.
%             Inputs are quant(ity) ('geoid' m, 'smd' kg/m^2 or 'water' m),
%             grid, n (grid size for 0.5° n = 360), h (altitude at which the 
%             function has to be evaluated, fil (filter in CS, SC or Colombo 
%             ordering formats, where a unique filter for every month can 
%             also be provided as a one-column cell array) & j (toggles 
%             between removing the normal field (1) and not removing it (0)).
% mn        - Long-term mean of the gravity field time-series. This has to 
%             be provided in CS-format.
% ds        - 'false' (no destriping) or 'true' (destriping with default 
%             values) or a vector of inputs for the DSTRPNGMTRX function.
%             Inputs are start degree, end degree, start order, end order,
%             and polynomial degree in that order.
% cindx     - Map containing the catchment indices. The size and grid must
%             be the same as the field that will be generated.
% indxnr    - If the values of only select catchments are desired, then their
%             indices must be provided here. [c,1]
%
% OUTPUT
% mest - Mass estimates for each catchment given in the index. This is
% given out as a cell array, with the rows of cell indicating the
% different months of the data
%
%--------------------------------------------------------------------------
% USES  SHbundle/gshs, sc2cs, cssc2clm
%       FilterBundle/dstrpngmtrx
%       GRACE/prpfltr
%--------------------------------------------------------------------------

% Created on 14 June 2007
% Author: Balaji Devaraju, Stuttgart
%--------------------------------------------------------------------------

m       = size(gsm,1);
if isscalar(gsm{1,8})
    lmax    = gsm{1,8};
else
    lmax    = size(gsm{1,9}) - 1;
end

if nargin == 1
    syn     = struct('quant','none','grid','block','n',180,'h',0,'fil',[],'j',0);
    mn      = 0;
    ds      = false;
    cindx   = [];
elseif nargin == 2
    mn      = 0;   
    ds      = false;
    cindx   = [];
elseif nargin == 3
    mn      = prpfltr(mn,lmax);
    ds      = false;
    cindx   = [];
elseif nargin == 4
    mn      = prpfltr(mn,lmax);
    if ~islogical(ds)
        dstrp   = ds;
        ds      = true;
    elseif ds
        dstrp   = [8 lmax 8 lmax 2];
    end
    cindx   = [];
elseif nargin == 5
    mn      = prpfltr(mn,lmax);
    if ~islogical(ds)
        dstrp   = ds;
        ds      = true;
    elseif ds
        dstrp   = [8 lmax 8 lmax 2];
    end
    indxnr = unique(cindx(:));
    indxnr = indxnr(indxnr ~= -9999);
elseif nargin == 6
    mn      = prpfltr(mn,lmax);
    if ~islogical(ds)
        dstrp   = ds;
        ds      = true;
    elseif ds
        dstrp   = [8 lmax 8 lmax 2];
    end
end

constants_grace % Loading the constants specific to GRACE.
                % We will only need the semi-major axis of
                % reference ellipsoid.

% Computes the monthly fields and the catchment aggregated values
if ~isempty(cindx) && ~isempty(indxnr)
    % Computing cell's Area
    dlam    = pi/syn.n;
    dt      = dlam;

    theta   = (0:dt:pi)';

    Area    = ae^2*dlam*dt*(cos(theta(1:end-1)) - cos(theta(2:end)))*ones(1,syn.n*2);

    c       = length(indxnr);
    ctch    = cell(c,1);
    for i = 1:c
        if ~any(find(cindx,indxnr(i)))
            indxnr(i)
        end
        ctch{i,1} = find(cindx == indxnr(i));
    end
    mest        = cell(m,6);
    mest(:,1:4) = gsm(:,4:7);

    if ~isempty(syn.fil) && ~iscell(syn.fil)
        fil = prpfltr(syn.fil,lmax);
    end

    for i = 1:m
        lmax    = size(gsm{1,9},1) - 1;
        if ds
            scfld = dstrpngmtrx(gsm{i,9}-mn,lmax,dstrp(5),dstrp(1),dstrp(2),dstrp(3),dstrp(4));
        else
            scfld = gsm{i,9} - mn;
        end
        if ~isempty(syn.fil)
            if iscell(syn.fil)
                fil = prpfltr(syn.fil{i},lmax);
                if ~isempty(fil)
                    display('Monthly filter successfully prepared')
                end
            end
            if ~isempty(fil)
                scfld = scfld.*fil;
                display('Filter successfully applied')
            end
        end
        [f,t,l]     = gshs(scfld,syn.quant,syn.grid,syn.n,syn.h,syn.j);
        mest{i,5}   = f;
        fArea       = f.*Area;
        mest{i,6}   = [indxnr,indxnr];
        for j = 1:c
            mest{i,6}(j,2)  = sum(fArea(ctch{j}))/sum(Area(ctch{j}));
        end
    end
% Computes only the monthly fields.
else
    mest        = cell(m,5);
    mest(:,1:4) = gsm(:,4:7);

    if ~isempty(syn.fil) && ~iscell(syn.fil)
        fil = prpfltr(syn.fil,lmax);
    end

    for i = 1:m
        if ds
            scfld = dstrpngmtrx(gsm{i,9}-mn,lmax,dstrp(5),dstrp(1),dstrp(2),dstrp(3),dstrp(4));
        else
            scfld = gsm{i,9} - mn;
        end
        if ~isempty(syn.fil)
            if iscell(syn.fil)
                fil = prpfltr(syn.fil{i},lmax);
                if ~isempty(fil)
                    display('Monthly filter successfully prepared')
                end
            end
            if ~isempty(fil)
                scfld = scfld.*fil;
                display('Filter successfully applied')
            end
        end
        [f,t,l] = gshs(scfld,syn.quant,syn.grid,syn.n,syn.h,syn.j);
        mest{i,5} = f;
    end
end
