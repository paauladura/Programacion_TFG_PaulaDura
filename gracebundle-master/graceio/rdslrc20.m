function c20 = rdslrc20(rfile,sfile)

% RDSLRC20 reads the C_{20} values from SLR provided in Technical note 7, 
% which can be used to replace the GRACE C_{20} coefficients.
%
% c20 = rdslrc20(rfile)
% c20 = rdslrc20(rfile,sfile)
%
% INPUT
% rfile - File (including path) from which the C20 coefficients have to be read
% sfile - File (including path) to which the resulting matrix have to be saved
%
% OUTPUT
% c20   - A six column matrix containing all the data in the technical note
%------------------------------------------------------------------------------
% USES: GRACEBundle/isleap
%------------------------------------------------------------------------------

% Balaji Devaraju, Stuttgart, 10 August 2011
%------------------------------------------------------------------------------

if nargin == 1
    if ~ischar(rfile)
        error('Filename is not a character array')
    end
    sv = false;
elseif nargin == 2
    if ~ischar(rfile)
        error('Read filename is not a character array')
    end
    if ~ischar(sfile)
        error('Save filename is not a character array')
    end
    sv = true;
end


fprintf('\n Reading SLR C20 coefficients ... ')

fid     = fopen(rfile,'r');
k       = 0;

n       = 300;
c20     = zeros(n,6);

while(~feof(fid))
    k       = k+1;
    tmpdat  = fgets(fid);
    if strncmp(tmpdat,'PRODUCT',7);
        n = k;
    end
    if k>n
        c20(k-n,1)  = str2double(tmpdat(1:5));
        c20(k-n,2)  = str2double(tmpdat(9:17));
        c20(k-n,3)  = str2double(tmpdat(20:39));
        c20(k-n,4)  = str2double(tmpdat(41:47))*1e-10;
        c20(k-n,5)  = str2double(tmpdat(50:55))*1e-10;
        if length(tmpdat)>56
            if strcmp('*',tmpdat(57))
                c20(k-n,6)  = 2;
            end
        else
            c20(k-n,6)  = 0;
        end
    end
end
c20(k-n+1:end,:)    = [];
c20                 = [c20(:,1:2) zeros(size(c20,1),2) c20(:,3:end)];
c20(:,3)            = floor(c20(:,2));
leapind             = isleap(c20(:,3));
c20(:,4)            = round((c20(:,2)-c20(:,3)) .* (double(leapind) + 365)) + 1;

if sv, save(sfile, '-v7', 'c20'); end

fprintf('done \n\n')
