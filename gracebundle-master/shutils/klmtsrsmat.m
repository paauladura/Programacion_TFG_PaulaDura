function [nmat,ordrng] = klmtsrsmat(gsm,lmax,jflag,mfix,dg)

% KLMTSRSMAT K_{lm} time-series matrix
% Converts the GSM cells of monthly/weekly gravity field data from GRACE to
% a matrix with the time-series of each coefficient occupying the columns
% of the matrix.
%
% nmat = klmtsrsmat(gsm,lmax)
% [nmat,ordrng] = klmtsrsmat(gsm,lmax,jflag,mfix,dg)
%
% INPUT
% gsm   -   Cells of monthly/weekly gravity field data as output by READGRC
%           function
% lmax  -   Maximum degree of expansion of the coefficients
% jflag -   A '0' flag indicates normal-field removed/only conversion
%           '1' flag indicates normal-field must be removed before
%           conversion. [def. '0']
% mfix  -   If '1' adds rows for missing months or weeks and sets the
%           coefficients to 'NaN' [def. '0']
% dg    -   If '1' arranges the coefficients in degree-leading order. 
%
% OUTPUT
% nmat  -   Matrix with the dimensions [time x k_{lm}]
% ordrng-   Vector with the ordering of the coefficients [deg ord]
%--------------------------------------------------------------------------
% USES cssc2clm, normalklm, monthfix
%--------------------------------------------------------------------------

% Created on: 15 February 2009, Stuttgart
% Author: Balaji Devaraju
%--------------------------------------------------------------------------

if nargin < 2,
    error('Insufficient inputs')
end

if nargin == 2,
    jflag = 0; mfix = 0; dg = 0;
elseif nargin == 3,
    mfix = 0; dg = 0;
elseif nargin == 4,
    dg = 0;
end

r = size(gsm,1);

nmat = [cell2mat(gsm(:,4:6)) zeros(r,(lmax+1)^2)];
nmat(:,5) = [];

if jflag==1,
    nrml = normalklm(lmax);
    if dg==1,
        temp = cssc2clm(full(nrml),lmax);
        [~,idx] = sort(temp(:,1));
        ordrng = temp(idx,1:2);
        ordrng = [ordrng; ordrng(ordrng(:,2)~=0,:)];
        for i = 1:r,
            temp = cssc2clm(gsm{i,9}-nrml,lmax);
            temp = temp(idx,:);
            nmat(i,5:end) = [temp(:,3); temp(temp(:,2)~=0,4)]';
        end
    else
        temp = cssc2clm(full(nrml),lmax);
        ordrng = [temp(:,1:2); temp(lmax+2:end,1:2)];
        for i = 1:r,
            temp = cssc2clm(gsm{i,9}-nrml,lmax);
            nmat(i,5:end) = [temp(:,3); temp(lmax+2:end,4)]';
        end
    end
else
    if dg==1,
        temp = cssc2clm(gsm{1,9},lmax);
        [~,idx] = sort(temp(:,1));
        ordrng = temp(idx,1:2);
        ordrng = [ordrng; ordrng(ordrng(:,2)~=0,:)];
        temp = temp(idx,:);
        nmat(1,5:end) = [temp(:,3); temp(temp(:,2)~=0,4)]';
        for i = 2:r,
            temp = cssc2clm(gsm{i,9},lmax);
            temp = temp(idx,:);
            nmat(i,5:end) = [temp(:,3); temp(temp(:,2)~=0,4)]';
        end
    else
        temp = cssc2clm(gsm{1,9},lmax);
        ordrng = [temp(:,1:2); temp(lmax+2:end,1:2)];
        nmat(1,5:end) = [temp(:,3); temp(lmax+2:end,4)]';
        for i = 2:r,
            temp = cssc2clm(gsm{i,9},lmax);
            nmat(i,5:end) = [temp(:,3); temp(lmax+2:end,4)]';
        end
    end
end

if mfix==1
    nmat = monthfix(nmat);
end
