function dv = degvar(spctrm,lmax,cum,quant,h,type)

% DEGVAR(SPCTRM,LMAX) computes the spherical harmonic spectral degree-
% variances.
%
% dv = degvar(spctrm,lmax)
% dv = degvar(spctrm,lmax,cum,quant,h,type)
%
% INPUT
% spctrm    -   Co-efficients of the spectral field given in sc- or
%               cs-formats or in [l m Clm Slm] format.
% lmax    	-   Maximum degree of the co-efficients.
% cum       -   Cumulative degree variances are computed if cum=1, else
%               degree variances are alone calculated.
% quant     -   String that defines the quantity to be calculated.
%               'none' - dimensionless, 'geoid' [m], 'smd' [kg/m^2],
%               'water' [m], 'potential' [m^2/s^2], 'dg' or 'gravity'
%               (gravity anomalies), 'tr' (gravity disturbances),'trr'
%               (d^2/dr^2), 'slope' (slope of the surface gradient).
%               def - ['none']
% h         -   height of the point of evaluation [m]
% type      -   'champ', 'grace', 'none'
%
% OUTPUT
% dv        -   Degree variances (power spectrum) of the spectral field.
%--------------------------------------------------------------------------

% Created on 29 May 2007
% Modified on 9 January 2009
%   - added additional inputs for applying eigengrav function
%   - removed FOR loops
% Author: Balaji Devaraju, Stuttgart
%--------------------------------------------------------------------------

% dv = [(0:lmax)',zeros(lmax+1,1)];
[rows, cols] = size(spctrm);

if cols == 4 || ((2*rows - 1) == cols) || (rows == cols)

    if nargin == 2
        cum = 0;
        quant = 'none';
        h = 0;
        type = 'none';
    elseif nargin == 3
        quant = 'none';
        h = 0;
        type = 'none';
    elseif nargin == 4
        h = 0;
        type = 'none';
    elseif nargin == 5
        type = 'none';
    end

    cum = logical(cum);

    switch type
    case 'champ'
        constants_champ
        type = [GM ae];
    case 'grace'
        constants_grace
        type = [GM ae];
    case 'none'
        type = [];
    otherwise
        error('Unknown TYPE for constants')
    end


    if cols == 4
        if length(spctrm) == sum(1:lmax+1)
            spctrm = sc2cs(gcoef2sc(cssc2clm(spctrm,120)));
        else
            error('Incomplete set of co-efficients')
        end
    elseif (2*rows-1 == cols)
        spctrm = sc2cs(spctrm);
    end

    if size(spctrm,1) == size(spctrm,2)
        dv 		= tril(spctrm);
        spctrm 	= spctrm - dv;
        dv 		= (sum(dv'.^2) + sum(spctrm.^2))';
        if nargin > 3
            dv = dv.*(eigengrav((0:lmax)',quant,h,type)).^2;
        end
        if cum
            dv = cumsum(dv);
        end
        dv = [(0:lmax)' dv];
    end
else
    error('The format is not compatible with the function')
end
