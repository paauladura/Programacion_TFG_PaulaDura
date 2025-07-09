function [P, dP, ddP,DegreeOrder] = legendreP(kmax, thetaRAD)
% IN:
%    kmax ........ maximum degree of the Legendre functions        [1 x 1]
%    thetaRAD .... co-latitude                                [rad][N x 1]
%
% OUT:
%    P ........... normalized Legendre functions                   [N x M]                             
%    dP .......... first order derivatives of Legendre functions   [N x M]              
%    ddP ......... second order derivatives of Legendre functions  [N x M]               
%    DegreeOrder.. matrix of [degree|order] per row                [N x 2]
%
% REMARKS:
%     1. P, dP, dPP will be 2D matrices
%     2. all 2D matrix indices are [ N(length of thetaRAD) * M(kmax+1) ]
%     3. if dP, dPP are not specified as outputs, then won't be computed
%     4. if thetaRAD is not Nx1, it will be reshaped

% -------------------------------------------------------------------------
% project: SHBundle 
% -------------------------------------------------------------------------
% authors:
%    Hailong FU (HF), Uni Stuttgart
%    Matthias ROTH (MR), GIS, Uni Stuttgart
%    <bundle@gis.uni-stuttgart.de>
% -------------------------------------------------------------------------
% revision history:
%    2014-10-27: MA, added matrix of [degree|order]
%    2013-06-11: MR, source code beautification, brush up help text, small
%                    speed-up (~10%)
%    2013-03-12: HF, initial version
% -------------------------------------------------------------------------
% license:
%    This program is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the  Free  Software  Foundation; either version 3 of the License, or
%    (at your option) any later version.
%  
%    This  program is distributed in the hope that it will be useful, but 
%    WITHOUT   ANY   WARRANTY;  without  even  the  implied  warranty  of 
%    MERCHANTABILITY  or  FITNESS  FOR  A  PARTICULAR  PURPOSE.  See  the
%    GNU General Public License for more details.
%  
%    You  should  have  received a copy of the GNU General Public License
%    along with Octave; see the file COPYING.  
%    If not, see <http://www.gnu.org/licenses/>.
% -------------------------------------------------------------------------

%% check number of output parameters
if nargout >= 2
    test_dP = true;
else
    test_dP = false;
end

if nargout >= 3
    test_ddP = true;
else
    test_ddP = false;
end


%% 2D parallel recursion
dmax             = kmax + 1;
Wlm_1_p(1, dmax) = 0; % faster than "Wlm_1_p = zeros(1, dmax);"
Wlm_1(1, dmax)   = 0; % faster than "Wlm_1 = zeros(1, dmax);"
Wlm_2(1, dmax)   = 0; % dito
Wlm_1_p(2)       = sqrt(3);
Wlm_1_p(1)       = sqrt(3);

siz      = size(thetaRAD); % reshape co-latitude to be one column
lenLat   = siz(1) * siz(2);
thetaRAD = reshape(thetaRAD, lenLat, 1);

if nargout >= 2
    slowMode = false; % there is a fast but singular solution, and also a slow non-singular solution
    if min(thetaRAD) < deg2rad(1) || max(thetaRAD) > deg2rad(179)
        slowMode = true;
    end
end

coslat = cos(thetaRAD);
sinlat = sin(thetaRAD); 
cotlat = cot(thetaRAD);

P(lenLat, (1 + dmax) * dmax / 2) = 0.0; % faster than "P = zeros(lenLat, (1 + dmax) * dmax / 2);"
P(:, 1) = 1;
P(:, 2) = sqrt(3) * coslat;
P(:, 3) = sqrt(3) * sinlat;

if test_dP
    dP(lenLat, (1 + dmax) * dmax / 2) = 0.0; % faster than "dP = zeros(lenLat, (1 + dmax) * dmax / 2);"
    dP(:, 2) = -P(:, 3);
    dP(:, 3) =  P(:, 2);
    if test_ddP
        ddP(lenLat, (1 + dmax) * dmax / 2) = 0.0; % faster than "ddP = zeros(lenLat, (1 + dmax) * dmax / 2);"
        ddP(:, 2) = -P(:, 2);
        ddP(:, 3) = -P(:, 3);
    end    
end

preCols = 3; % previously finished columns

for l = 2:kmax % if l >= 2
    lp1 = l + 1;  
    m = 0:(l - 1);  
    mp1 = m + 1;
    Wlm_1(mp1)    = sqrt((4 * l^2 - 1) ./ (l^2 - m.^2));
    Wlm_1(lp1)    = sqrt(1 + 0.5 / l);  
    mp1_          = 1:(l - 1);  
    Wlm_2(mp1_)   = Wlm_1(mp1_) ./ Wlm_1_p(mp1_); 
    absMp1        = preCols + mp1;    
    absMp1_pre    = absMp1 - l;   
    absMp1_prepre = absMp1_pre - (l - 1);  
    
    % Legendre
    P(:, absMp1)        = bsxfun(@times, bsxfun(@times, P(:, absMp1_pre), coslat),  Wlm_1(mp1)) ...
                          - bsxfun(@times, P(:, absMp1_prepre), Wlm_2(mp1));         
    P(:, preCols + lp1) = Wlm_1(lp1) * bsxfun(@times, P(:, preCols), sinlat);

    if test_dP % first order directive
        if slowMode
            dP(:, absMp1)        = bsxfun(@times, bsxfun(@times, dP(:, absMp1_pre), coslat) - bsxfun(@times, P(:, absMp1_pre), sinlat),  Wlm_1(mp1)) ...
                                 - bsxfun(@times, dP(:, absMp1_prepre), Wlm_2(mp1));          
            dP(:, preCols + lp1) = Wlm_1(lp1) * (bsxfun(@times, dP(:, preCols), sinlat) + bsxfun(@times, P(:, preCols), coslat));        
 
            if test_ddP % second order directive
                ddP(:, absMp1)        = bsxfun(@times, bsxfun(@times, (ddP(:, absMp1_pre) - P(:, absMp1_pre)), coslat) - bsxfun(@times, dP(:, absMp1_pre), sinlat * 2),  Wlm_1(mp1) )...
                                      - bsxfun(@times, ddP(:, absMp1_prepre), Wlm_2(mp1)); 
                ddP(:, preCols + lp1) = Wlm_1(lp1) * ( bsxfun(@times, (ddP(:, preCols) - P(:, preCols)), sinlat) + bsxfun(@times, dP(:, preCols), coslat * 2) );        
            end             
        else
            absM        = preCols + (1:l+1);
            absM_pre    = absM - l;
            e           = sqrt((2*l+1) / (2*l-1) * (l*l - (0:l) .* (0:l))); 
            dP(:, absM) = bsxfun(@times, P(:, absM), l*cotlat) - (1./sinlat*e) .* P(:, absM_pre);
            
            if test_ddP % second order directive
                ddP(:, absM) = bsxfun(@times, dP(:, absM), -cotlat) + (1./sinlat./sinlat*((0:l).*(0:l))-l*(l+1)).*P(:, absM);
            end
        end
    end   

    Wlm_1_p = Wlm_1;
    preCols = preCols + lp1;
end

%% transpose that output is the same as Markus'
P = P';
if test_dP
    dP = dP';
    if test_ddP
        ddP = ddP';
    end
end


if nargout>3
    % create a matrix of [degree|order] for this kind of sorting 
    E = triu(ones(kmax+1));
    D = diag([0:kmax]);
    n = E*D;
    m = D*E;
    n = n + tril(NaN+ ones(kmax+1),-1);

    n = n(:);
    m = m(:);
    kk = isnan(n);
    n(kk) = [];
    m(kk) = [];

    DegreeOrder = [n,m];
end
