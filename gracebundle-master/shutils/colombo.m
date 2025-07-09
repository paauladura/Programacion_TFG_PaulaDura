function [mord,ind] = colombo(lord,lmax)

% COLOMBO performs the colombo ordering for a given spherical harmonic
% coefficients in the form of [l m Clm Slm] vector sorted along degree (l).
% Colombo ordering refers to the ordering of the coefficients by sorting
% along order (m).
%
% mord = colombo(lord,lmax)
%
% INPUT
% lord      -   [l m Clm Slm] or [l m Clm; l m Slm]
% lmax    	-   maximum degree of the spherical harmonic development [optional]
%               This is more for checking the integrity of the SH spectrum 
%               provided, and also for sorting in case the spectrum is given in
%               three columns.
%
% OUTPUT
% mord      -   [l m Clm Slm] sorted along 'm' column.
% ind 		- 	Index of the sorting
%--------------------------------------------------------------------------

% Created on 17 September 2007, Stuttgart
% Author: Balaji Devaraju


% Revision:
% BD    25/07/2014  Fixed a bug due to improper input argument checking.
%--------------------------------------------------------------------------

[r,c] = size(lord);

if c < 2
    error('LORD must have atleast two columns. Please verify your input.')
end

if nargin == 2
	degsum  = sum(1:lmax+1);
	totelem = (lmax+1)^2;
    if r == degsum
        [msort,ind] = sort(lord(:,2));
        if msort == lord(:,2)
            mord = lord;
        else
            mord = lord(ind,:); 
        end
    elseif r == totelem
        mord 					= [lord(1:degsum,:), zeros(degsum,1)];
        mord((mord(:,2)~=0),4) 	= lord(degsum+1:end,3);
        [msort,ind] 			= sort(mord(:,2));
        mord 					= mord(ind,:);
    else
        error('Please verify the inputs.')
    end
elseif nargin == 1
    [msort,ind] = sort(lord(:,2));
    if msort == lord(:,2)
        mord = lord;
    else
        mord = lord(ind,:);
    end
else
	error('Please check the number of inputs.')
end

