function Q = readgrccov(in_file, out_file)

% READGRCCOV reads GRACE covariance matrices from text files as
% provided by the analysis centers.
%
% Q = readgrccov(in_file)
% Q = readgrccov(in_file, out_file)
%
% INPUT
% in_file   -   Filename with absolute path of the GACE covariance text
%               file
% out_file  -   Optional output filename to save the matrix as a Matlab 
%               MAT file.
%
% OUTPUT
% Q         -   GRACE covariance matrix
%----------------------------------------------------------------------

if nargin == 1
    if ~ischar(in_file)
        error('Input filename is not a string!')
    end
    savecov = false;
elseif nargin == 2
    if ~ischar(in_file)
        error('Input filename is not a string!')
    end
    if ischar(out_file)
        savecov = true;
    else
        error('Output filename is not a string!')
    end
end

fid = fopen(in_file, 'r+');


fprintf('Reading header data ...\n')
loop    = true;
readdata= false;
chkl    = true;
nhead   = 0;
nlines  = 0;
j       = 1;
J       = 0;
L       = [];

while loop
    str     = fgets(fid);
    if str ~= -1
        nlines  = nlines + 1;
        str = str(1:end-1);
        if readdata
            N = length(str)/12;
            J = J + N;
            % keyboard
            tmp = str2double(reshape(str', 12, N)');

            Q(idx(j:j+N-1)) = tmp;
            j = j + N;
            % fprintf('\b\b\b\b\b\b\b\b\b\b\b\b%12d', j)
%            if j > sum(1:(lmax+1)^2)
%                keyboard
%                fprintf('\n%s\n', str)
%            end
        else
            if chkl
                L       = regexp(str, 'ORDER');
            end
            lmax_str = num2str(lmax);
            Sk      = regexp(str, 'S');
            k       = regexp(str, [lmax_str, ' ', lmax_str]);
            nhead   = nhead + 1;
            if chkl && ~isempty(L)
                lmax    = sqrt(strread(str, '%*s %f')) - 1;
                nelem   = (lmax+1)^2;
                idx     = nelem * (0:nelem-1);
                idx     = bsxfun(@plus, (1:nelem)', idx);
                idx     = tril(idx)';
                idx     = idx(:);
                idx     = idx(idx~=0);

                Q       = zeros((lmax+1)^2);
                fprintf('Maximum number of spherical harmonic degrees = %8f\n', lmax)
                chkl    = false;
            end
            if ~isempty(Sk) && ~isempty(k)
                readdata = true;
                fprintf('Read %12d header lines ...\n', nhead)
                %fprintf('Reading the covariance matrices: 000000000000')
            end
        end
    else
        fprintf('Last record in the file %12d\n', str)
        fprintf('Read %12d lines\n', nlines)
        fprintf('Read %12d covariance matrix elements\n', J)
        loop     = false;
        readdata = false;
    end
end
fclose(fid);

% fprintf('Reading the covariance matrices ...\n')
% k_lm_nk = textread(fnm, '%12.6f', 'headerlines', nhead);
