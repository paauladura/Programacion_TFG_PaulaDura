function [dsfld,Wds] = dstrpngmtrx(fld,lmax,ply,sord,eord,sdeg,edeg)

% DSTRPNGMTRX computes the destriping matrix of the destriping filter
% proposed by Swenson and Wahr (2006).
%
% [dsfld,Wds] = dstrpngmtrx(fld,lmax)
% [dsfld,Wds] = dstrpngmtrx(fld,lmax,ply)
% [dsfld,Wds] = dstrpngmtrx(fld,lmax,ply,sord,eord,sdeg,edeg)
%
% INPUT
% fld   -   Co-efficients from the time-variable GRACE gravity-field. The
%           co-efficients must be in either CS, SC, or [l m Clm Slm]
%           formats.
% lmax  -   Maximum degree of spherical harmonic expansion.
% ply   -   Polynomial degree for the polynomial that will be fit to the
%           coefficients during destriping.
%           [optional; default: 2]
% sord  -   Order starting from which the destriping must be applied.
%           [optional; default order: 8]
% eord  -   Order until which the destriping must be applied.
%           [optional; default order: lmax-8]
% sdeg  -   Degree starting from which the destriping must be applied.
%           [optional; default degree: sord]
% edeg  -   Degree until which the destriping must be applied.
%           [optional; default degree: lmax]
%
% OUTPUT
% Wds   -   Destriping matrix
% dsfld -   Destriped coefficients provided in the same format as the 
%           input. However with [l m Clm Slm] format as input the output is
%           same as input, but the output is Colombo-ordered.
%
%--------------------------------------------------------------------------

% Created on: 15 September 2008, Stuttgart
% Author: Balaji Devaraju
%--------------------------------------------------------------------------

[rows,cols] = size(fld);
if isequal(cols,4) && isequal(rows, sum(1:lmax+1))
    fld = colombo(fld);
    fld = clm2sc(fld);
elseif isequal(rows,lmax+1) && isequal(rows,cols)
    fld = cs2sc(fld);
elseif ~isequal(rows,lmax+1) && ~isequal(cols,(2*lmax+1))
    error('Error matrix not in the SC, CS, or Colombo ordering formats')
end

if nargin == 2 
    sord = 8;
    eord = lmax-10;
    sdeg = sord;
    edeg = lmax;
    ply  = 2;
elseif nargin == 3
    sord = 8;
    eord = lmax-10;
    sdeg = sord;
    edeg = lmax;
elseif nargin == 4
    eord = lmax-10;
    sdeg = sord;
    edeg = lmax;
elseif nargin == 5
    sdeg = sord;
    edeg = lmax;
elseif nargin == 6
    edeg = lmax;
elseif nargin < 2
    error('Insufficient input arguments')
end

if isempty(ply), ply = 2; end
if isempty(edeg), edeg = lmax; end
if isempty(sord) || (sord == 0), sord = 8; end
if isempty(eord) || (eord > (edeg-10)), eord = edeg-10; end
if isempty(sdeg) || (sdeg < sord), sdeg = sord; end
% if edeg <= eord, edeg = eord + 10; end

Wds = speye((lmax+1)^2);
Wds(Wds==1) = 0;
dsfld = fld;
cindx2 = cumsum(fliplr(1:lmax+1))';
cindx1 = [1;cindx2(1:end-1)+1];
sindx2 = cumsum(fliplr(1:lmax))' + cindx2(end);
sindx1 = [cindx2(end)+1;sindx2(1:end-1)+1];


if sord > 0
    Wds(cindx1(1):cindx2(1),cindx1(1):cindx2(1)) = eye(lmax+1);

    for m=1:sord-1
        Wds(cindx1(m+1):cindx2(m+1),cindx1(m+1):cindx2(m+1)) = eye(lmax+1-m);
        Wds(sindx1(m):sindx2(m),sindx1(m):sindx2(m)) = eye(lmax+1-m);
    end
end

if eord < lmax
    for m = eord:lmax
        Wds(cindx1(m+1):cindx2(m+1),cindx1(m+1):cindx2(m+1)) = eye(lmax+1-m);
        Wds(sindx1(m):sindx2(m),sindx1(m):sindx2(m)) = eye(lmax+1-m);
    end
end

for m = sord:eord
    Word = zeros(lmax+1-m);
    
    nsmooth = round(30*exp(-m/10)+1);
    if nsmooth < 5
        nsmooth = 5;
    end
    hnsmooth = floor(nsmooth/2);
    % nsmooth = nsmooth-1;
    
    if sdeg>m
        lvec = sdeg:edeg;
    else
        lvec = m:edeg;
    end
    
    for l = 1:length(lvec)
        if mod(nsmooth,2) == 0
            if lvec(l)-(2*(hnsmooth-1)) < lvec(1)
                degvec = lvec(l) + 2*(-(hnsmooth-1):hnsmooth);
                len = length(find(degvec < lvec(1)));
                degvec = lvec(l) + 2*(-(hnsmooth-1-len):hnsmooth+len);
            elseif lvec(l)+(2*hnsmooth) > edeg
                degvec = lvec(l) + 2*(-(hnsmooth-1):hnsmooth);
                len = length(find(degvec > edeg));
                degvec = lvec(l) + 2*(-(hnsmooth-1+len):hnsmooth-len);
            else
                degvec = lvec(l)+2*(-(hnsmooth-1):hnsmooth);
            end
        else
            if lvec(l)-(2*(hnsmooth)) < lvec(1)
                degvec = lvec(l) + 2*(-hnsmooth:hnsmooth);
                len = length(find(degvec < lvec(1)));
                degvec = lvec(l) + 2*(-(hnsmooth-len):hnsmooth+len);
            elseif lvec(l)+(2*hnsmooth) > edeg
                degvec = lvec(l) + 2*(-hnsmooth:hnsmooth);
                len = length(find(degvec > edeg));
                degvec = lvec(l) + 2*(-(hnsmooth+len):hnsmooth-len);
            else
                degvec = lvec(l)+2*(-hnsmooth:hnsmooth);
            end
        end
        degvec = degvec(degvec>=sdeg & degvec<=edeg & degvec>=m);
        V = vander(degvec);
        V = fliplr(V);
        V = V(:,1:ply+1);
        vlm = V((degvec==lvec(l)),:);
        M = zeros(length(degvec),lmax+1-m);
        for i = 1:length(degvec)
            M(i,degvec(i)+1-m) = 1;
        end
        Word(lvec(l)+1-m,:) = vlm*(V\M);
    end
    Word = (eye(size(Word))-Word);
    dsfld(m+1:end,lmax+1-m) = Word*dsfld(m+1:end,lmax+1-m);
    dsfld(m+1:end,lmax+1+m) = Word*dsfld(m+1:end,lmax+1+m);
    Wds(cindx1(m+1):cindx2(m+1),cindx1(m+1):cindx2(m+1)) = Word;
    Wds(sindx1(m):sindx2(m),sindx1(m):sindx2(m)) = Word;
end

save('tempmat.mat','Wds')
Wds = covord('tempmat.mat','./',lmax,1);
delete('tempmat.mat')

if isequal(cols,4) && isequal(rows, sum(1:lmax+1))
    dsfld = cssc2clm(dsfld,lmax);
elseif isequal(rows,lmax+1) && isequal(rows,cols)
    dsfld = sc2cs(dsfld);
end
