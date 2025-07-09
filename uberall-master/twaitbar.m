function hand_out = twaitbar(perc, hand_in, tag)

% TWAITBAR displays an ASCII waitbar in the command window. 
% It is only possible to use a single twaitbar at the same time. Also,
% avoid output to the command window as long as a twaitbar is active --
% output would interfere with twaitbar! 
% IN:
%   perc ...... initialization: 'init', finalization; 'close' 
%               or progress ratio, i.e. a number between 0 and 1
%   hand_in ... waitbar handle; 
%               initialization: takes the length of the waitbar or an empty 
%               array (default: 30)
%   tag ....... initialization: takes a text for the waiting bar
%               (default: '')
% OUT:
%   hand_out .. waitbar handle (stores some information for waitbar writing 
%               process)  
% EXAMPLE:
%   N = 1000; % number of data samples
%   wait_h = twaitbar('init', 40, 'computation'); % create waitbar of length 40
%   for i = 1:N
%       wait_h = twaitbar(i / N, wait_h);
%       % ...
%   end
%   twaitbar('close', wait_h);

% -------------------------------------------------------------------------
% project: SHBundle 
% -------------------------------------------------------------------------
% authors:
%    Matthias Roth (MR), GI, Uni Stuttgart               
%    <bundle@gis.uni-stuttgart.de>
% -------------------------------------------------------------------------
% revision history:
%    2015-05-20: MR, initial version
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

switch perc
    case 'init' % initialize handle
        % set default waitbar length
        if ~exist('hand_in', 'var') || isempty(hand_in) || ~isnumeric(hand_in)
            hand_out.barlen = 30;
        else 
            hand_out.barlen = hand_in;
        end
        hand_out.nchar = hand_out.barlen + 7; % add length of additional letters
        hand_out.k = 1;
        hand_out.p = 1;
        % print tag
        if exist('tag', 'var')
            fprintf('%s ... ', tag);
        end
        % print empty bar
        fprintf('%3d%% [%s]', 0, repmat(' ', 1, hand_out.barlen));
        return
    case 'close' % finalize handle
        hand_out = hand_in;
        % delete waitbar line
        fprintf([repmat('\b', 1, hand_out.nchar) 'done.\n']);
        hand_out.nchar = 0;
        return
    otherwise
        hand_out = hand_in;
end

% write waitbar line
while hand_out.p < fix(perc * 100)
    fprintf([repmat('\b', 1, hand_out.nchar) '%3d%% [%s%s]'], fix(perc * 100), ...
        repmat('=', 1, hand_out.k), repmat(' ', 1, hand_out.barlen - hand_out.k));
    while hand_out.k < fix(perc * hand_out.barlen)
        hand_out.k = hand_out.k + 1;
    end
    hand_out.p = hand_out.p + 1;
end
