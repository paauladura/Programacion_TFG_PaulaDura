function [f, varargout] = tsplot(nrplots, tscale, datatype, indxs, varargin)

% TSPLOT plots the time-series of multiple datasets at multiple locations
% using the SUBPLOT function. 
%
% INPUT
% nrplots   - number of plots that should be used in the subplot
% tscale    - Time period specification [from_month from_year to_month to_year].
% 

if strcmp(datatype, 'full')
    sind = 2;
    if isequal(tscale,[1 12])
        tmeinfo = 1;
    else
        tmeinfo = 3;
    end
    if length(indxs) == 1
        ccol = zeros(length(varargin),1);
        for i = 1:length(varargin)
            ccol(i) = find(varargin{i}(1,:) == indxs);
        end
    else
        ccol = zeros(length(varargin),length(indxs));
        for i = 1:length(varargin)
            for j = 1:length(indxs)
                ccol(i,j) = find(varargin{i}(1,:) == indxs(j));
            end
        end
    end
elseif strcmp(datatype, 'normal')
    sind = 1;
    tmeinfo = 1;
    ccol = ones(length(varargin),1)*2;
    % tlte = [ ];
end

if isequal(tscale,[1 12])
    mnths = mnthnms('vshort');
end



clr = [   0    0    0;
         30   60  255;
        250   60   60;
          0  220    0;
        240  130   40;
          0  200  200;
        230  220   50;
        160    0  200;
        160  230   50;
          0  160  255;
        240    0  130;
        230  175   45;
          0  210  140;
        130    0  220]/255; 
          
        
scrsz = get(0,'ScreenSize');
f = figure('OuterPosition',[1 scrsz(4)/2 scrsz(3)/3 scrsz(4)/2]);

if nrplots == 1
    h = subplot(1,1,1);
    for i = 1:length(varargin)
        
        plot(varargin{i}(sind:end, tmeinfo), varargin{i}(sind:end, ccol(i)), 'Color', clr(i,:), 'Linewidth', 1);
        hold on
    end
    datetick('x')
    if  ~isequal(tscale,[1 12])
        set(gca, 'xlim',  [datenum(tscale(1), 1, 15) datenum(tscale(2), 12, 15)]);
    else
        set(gca, 'xtick',  1:1:12);
        set(gca, 'xticklabel', mnths);
        set(gca, 'xlim', [0.5 12.5]);
    end
    
    
else
    h = zeros(nrplots(1)*nrplots(2),1);
    for i = 1:nrplots(1)*nrplots(2)
        h(i) = subplot(nrplots(1), nrplots(2), i);
        for j = 1:length(varargin)
            plot(varargin{j}(sind:end, tmeinfo), varargin{j}(sind:end, ccol(j, i)), 'Color', clr(j,:), 'Linewidth', 1);
            hold on
        end
        datetick('x');
        if  ~isequal(tscale,[1 12])
            set(gca, 'xlim',  [datenum(tscale(1), 1, 15) datenum(tscale(2), 12, 15)]);
        else
            set(gca, 'xlim',  [1 12]);
            set(gca, 'xticklabel', mnths);
        end
    end
    
end

varargout{1} = h;
