function m_names = monthnames(stle)

% MONTHNAMES provides the names of months in different formats for plotting 
% purposes.
%
% m_names = monthnames(style)
%
% INPUT
% style  - Any of 'long', 'short', 'vshort' and 'ssnl'. ['short']
%
% OUTPUT
% m_names - Names of months in the required format given as a cell-array.
%------------------------------------------------------------------------------

% Christof Lorenz. Garmisch-Partenkirchen.
%------------------------------------------------------------------------------


if nargin < 1
    stle = 'short';
end

if strcmp(stle, 'long')
    m_names{1,1} = 'January';
    m_names{2,1} = 'February';
    m_names{3,1} = 'March';
    m_names{4,1} = 'April';
    m_names{5,1} = 'May';
    m_names{6,1} = 'June';
    m_names{7,1} = 'July';
    m_names{8,1} = 'August';
    m_names{9,1} = 'September';
    m_names{10,1} = 'October';
    m_names{11,1} = 'November';
    m_names{12,1} = 'December';
elseif strcmp(stle, 'short')
    m_names{1,1} = 'Jan';
    m_names{2,1} = 'Feb';
    m_names{3,1} = 'Mar';
    m_names{4,1} = 'Apr';
    m_names{5,1} = 'May';
    m_names{6,1} = 'Jun';
    m_names{7,1} = 'Jul';
    m_names{8,1} = 'Aug';
    m_names{9,1} = 'Sep';
    m_names{10,1} = 'Oct';
    m_names{11,1} = 'Nov';
    m_names{12,1} = 'Dec';
elseif strcmp(stle, 'vshort')
    m_names{1,1} = 'J';
    m_names{2,1} = 'F';
    m_names{3,1} = 'M';
    m_names{4,1} = 'A';
    m_names{5,1} = 'M';
    m_names{6,1} = 'J';
    m_names{7,1} = 'J';
    m_names{8,1} = 'A';
    m_names{9,1} = 'S';
    m_names{10,1} = 'O';
    m_names{11,1} = 'N';
    m_names{12,1} = 'D';
elseif strcmp(stle, 'ssnl')
    m_names{1,1} = 'DJF';
    m_names{2,1} = 'MAM';
    m_names{3,1} = 'JJA';
    m_names{4,1} = 'SON';
end


    
    
