function coeff = readgeocentercoeff(filename)
  
% READGEOCENTERCOEFF reads the degree 1 spherical harmonic coefficients
% time-series that describe the geocenter motion.
%
% coeff = readgeocentercoeff(fname)
%
% INPUT
% filename - Filename with absolute path
%-----------------------------------------------------------------------



% 
%  Project: GRACE Bundle
%  Copyright Balaji Devaraju (BD) 
%  devaraju at ife dot uni-hannover dot de
% 
%  License: GNU GPLv3 or later
%  You should have received a copy of the GNU General Public License
%  along with EGRAFS;  If not, see <http://www.gnu.org/licenses/>.
%  
%  Authors: Balaji Devaraju
% 
%  Version control
%  Auhtor   YYYY:MM:DD  Comment
%   BD      2019:07:12  Initial version
%----------------------------------------------------------------------


narginchk(1, 1)

%% check parameters
if ~exist(filename, 'file')
    error('File "%s" not found...', filename);
end
fprintf('Parsing file "%s"...\n', filename);

% read header
fid       = fopen(filename, 'r');
hasErrors = true;
hlines    = 0;
while 1
    one_line = fgets(fid);
    if ~ischar(one_line), break, end
    hlines  = hlines + 1;
    if isempty(one_line)
        continue
    else
        keyword = strread(one_line,'%s');

        idx1 = strcmp('end', keyword);
        idx2 = strcmp('header', keyword);
        idx3 = (any(idx1) && any(idx2));
        
        if idx3
            fprintf('Found %s at line #%g\n', sprintf('%s %s %s', keyword{1}, keyword{2}, keyword{3}), hlines);
            break
        else
            continue
        end
    end
end

if isoctave
    [code, l, m, C, S, dC, dS, t0, t_end] = textread(filename, '%s %f %f %f %f %f %f %s %s', 'headerlines', hlines);
else
    [code, l, m, C, S, dC, dS, t0, t_end] = textread(filename, '%s %f %f %f %f %f %f %s %s %*[^\n]', 'headerlines', hlines);
end

coeff = cell(length(m)/length(unique(m)), 10);
 
for k = 1:2:length(l)
    n = (k+1)/2;
    coeff{n, 1} = 'TN-13';
    coeff{n, 2} = 'GSM';
    coeff{n, 3} = 6;
    
    if strcmp(t0{k}, t0{k+1}) && strcmp(t_end{k}, t_end{k+1})
        coeff{n, 4}   = str2num(t0{k}(1:4));
        coeff{n, 5}   = str2num(t0{k}(5:6));
        d0            = dayofyear(coeff{n, 4}, coeff{n, 5}, str2num(t0{k}(7:8)));
        d_end         = dayofyear(str2num(t_end{k}(1:4)), str2num(t_end{k}(5:6)), str2num(t_end{k}(7:8)));
        coeff{n, 6}   = [d0, d_end, d_end-d0];
        coeff{n, 7}   = 0;
        coeff{n, 8}   = max(l);
        coeff{n, 9}   = sc2cs(clm2sc([l(k:k+1) m(k:k+1) C(k:k+1) S(k:k+1)]));
        coeff{n, 10}  = sc2cs(clm2sc([l(k:k+1) m(k:k+1) dC(k:k+1) dS(k:k+1)]));
    else
        fprintf('Found mismatch between %s , %s , %s and %s in iteration %g\n', t0{k}, t_end{k}, t_0{k+1}, t_end{k+1}, k);
        error('Please verify the input file for possible errors')
    end
end
