function varargout = klm2tsrs(gsm, varargin)

% KLMTSRSMAT K_{lm} time-series matrix
% Converts the GSM cells of monthly/weekly gravity field data from GRACE to
% a matrix with the time-series of each coefficient occupying the columns
% of the matrix.
%
% klm_ts = klm2tsrs(gsm)
% [klm_ts, ordrng] = klm2tsrs(gsm, OPTIONS, VALUES)
% [klm_ts, dklm_ts, ordrng] = klm2tsrs(gsm, OPTIONS, VALUES)
%
% INPUT
% gsm   -   Cells of monthly/weekly gravity field data as output by READGRC
%           function
% OPTIONAL ARGUMENTS
%--------------------
% 'max_degree'  -   Maximum degree of expansion of the coefficients.
% 'normalfield' -   Toggle for subtracting normal field
%                       true    - do not subtract normal-field 
%                       false   - subtract normal-field [def. false]
% 'monthfix'    -   Toggle switch to add rows for missing months or weeks
%                   and sets the coefficients to 'NaN' [def. false]
% 'degreewise'  -   Arranges the coefficients in degree-leading order. 
%                   [C00 C10 C11 C20 C21 ... CLL S11 S21 S22 ... SLL]
%                   [def. false]
% 'error'       -   Toggle switch for converting the errors also as a
%                   time-series [def. false]
%
% OUTPUT
% klm_ts        -   Coefficients arranged as a time-series matrix with the 
% dklm_ts           dimensions [time x k_{lm}]
%                   [year month start_day end_day klm]
% ordrng        -   Vector with the ordering of the coefficients [deg ord]
%--------------------------------------------------------------------------
% USES cssc2clm, normalklm, monthfix
%--------------------------------------------------------------------------

% Created on: 18 September 2017, Hannover
% Author: Balaji Devaraju
%--------------------------------------------------------------------------

narginchk(1, inf)

% Default values for optional arguments
defpar = { 'max_degree',    gsm{1,8}; ...
            'normalfield',  false; ...
            'monthfix',     false; ...
            'degreewise',   false; ...
            'errors',       false  };

params = getopt(defpar, false, varargin);


r    = size(gsm,1);
lmax = params.max_degree;

% Initializing output
klm_ts      = [cell2mat(gsm(:,4:6)) zeros(r,(lmax+1)^2)];
klm_ts(:,5) = [];
if params.errors
    dklm_ts = klm_ts;
end

if params.normalfield
    nrml = normalklm(lmax);
    if params.degreewise
        temp    = cssc2clm(full(nrml),lmax);
        [~,idx] = sort(temp(:,1));
        ordrng  = temp(idx,1:2);
        ordrng  = [ordrng; ordrng(ordrng(:,2)~=0,:)];
        for i = 1:r,
            temp = cssc2clm(gsm{i,9}-nrml,lmax);
            temp = temp(idx,:);
            klm_ts(i,5:end) = [temp(:,3); temp(temp(:,2)~=0,4)]';

            if params.errors
                temp = cssc2clm(gsm{i,10},lmax);
                temp = temp(idx,:);
                dklm_ts(i,5:end) = [temp(:,3); temp(temp(:,2)~=0,4)]';
            end
        end
    else
        temp    = cssc2clm(full(nrml),lmax);
        ordrng  = [temp(:,1:2); temp(lmax+2:end,1:2)];
        for i = 1:r,
            temp = cssc2clm(gsm{i,9}-nrml,lmax);
            klm_ts(i,5:end) = [temp(:,3); temp(lmax+2:end,4)]';

            if params.errors
                temp = cssc2clm(gsm{i,10},lmax);
                dklm_ts(i,5:end) = [temp(:,3); temp(lmax+2:end,4)]';
            end
        end
    end
else
    if params.degreewise
        temp    = cssc2clm(gsm{1,9},lmax);
        [~,idx] = sort(temp(:,1));
        ordrng  = temp(idx,1:2);
        ordrng  = [ordrng; ordrng(ordrng(:,2)~=0,:)];
        temp    = temp(idx,:);
        klm_ts(1,5:end) = [temp(:,3); temp(temp(:,2)~=0,4)]';
        for i = 1:r,
            temp = cssc2clm(gsm{i,9},lmax);
            temp = temp(idx,:);
            klm_ts(i,5:end) = [temp(:,3); temp(temp(:,2)~=0,4)]';

            if params.errors
                temp = cssc2clm(gsm{i,10},lmax);
                temp = temp(idx,:);
                dklm_ts(i,5:end) = [temp(:,3); temp(temp(:,2)~=0,4)]';
            end
        end
    else
        temp = cssc2clm(gsm{1,9},lmax);
        ordrng = [temp(:,1:2); temp(lmax+2:end,1:2)];
        klm_ts(1,5:end) = [temp(:,3); temp(lmax+2:end,4)]';
        for i = 1:r,
            temp = cssc2clm(gsm{i,9},lmax);
            klm_ts(i,5:end) = [temp(:,3); temp(lmax+2:end,4)]';

            if params.errors
                temp = cssc2clm(gsm{i,10},lmax);
                dklm_ts(i,5:end) = [temp(:,3); temp(lmax+2:end,4)]';
            end
        end
    end
end

if params.monthfix
    klm_ts = monthfix(klm_ts);
end

if nargout == 1
    varargout{1} = klm_ts;
elseif nargout == 2
    varargout{1} = klm_ts;
    varargout{2} = ordrng;
elseif nargout == 3
    varargout{1} = klm_ts;
    if params.errors
        varargout{2} = dklm_ts;
    else
        error('Requesting DKLM_TS without switching on the option ERRORS')
    end
    varargout{3} = ordrng;
end
