function [otpt,cont,agg] = intannum(ts,type,y,clms)

% INTANNUM computes the desired statistic (mean or RMS or sum) from a 
% monthly time-series dataset with a minimum number of continuous 
% datapoints. This is mainly useful in computing the statistic of a 
% time-series only for integer number of years.
%
% [otpt,cont,agg] = intannum(ts,type,y,clms)
% 
% INPUT
% ts 	- 	Monthly time-series. The data gaps must be filled with NaN.
% type 	- 	'full' where the first row contains information on the 
% 			columns, or 'normal' if the time-series starts from the 
%			first row. [default: 'full']
% y 	- 	integer number of years. Must be a scalar value. [default:3]
% clms 	- 	In general the ordering expected is 
% 			[month year data], which is [1 2 3], but if the 
% 			ordering is [year month data] then [2 1 3]. The third 
% 			number indicates the column from which the data starts. 
% 			This must be given if the data does not start from the 
% 			third column. [default: 1 2 3]
%
% OUTPUT
% otpt 	- 	The statistics of the dataset.
% cont 	- 	Provides a logical matrix which indicates the datagaps.
% agg 	- 	Aggregate of the individual continuous parts.
%-----------------------------------------------------------------------

% Balaji Devaraju. 10 October 2012, Stuttgart.
%-----------------------------------------------------------------------

if nargin == 1
	type 	= 'normal';
	y		= 36;
	clms 	= [1 2 3];
elseif nargin == 2
	y 		= 36;
	clms 	= [1 2 3];
	if strcmp(type,'full')
		cinfo 	= ts(1,:);
		ts 		= ts(2:end,:);
	end
elseif nargin == 3
	clms 	= [1 2 3];
	y 		= y*12;
	if strcmp(type,'full')
		cinfo 	= ts(1,:);
		ts 		= ts(2:end,:);
	end
elseif nargin == 4
	y 		= y*12;
	if strcmp(type,'full')
		cinfo 	= ts(1,:);
		ts 		= ts(2:end,:);
	end
end

% Finding contiguous parts of the dataset
cont 	= double(~isnan(ts(:,clms(3):end)));
agg		= cont;
for k = 2:size(cont,1)
	cndtn 			= cont(k,:) - cont(k-1,:);
	agg(k,cndtn==0) = agg(k,cndtn==0) + agg(k-1,cndtn==0);
	agg(k,cndtn==1) = 1;
	agg(k,cndtn<0) 	= 0;
end

% Computing mean runoff from the dataset
id 				= find(max(agg) >= y);
m 				= max(agg);
otpt.mn			= [cinfo(clms(3):end); ts(1,clms(3):end)];
otpt.mn(2,:) 	= NaN;
otpt.sum 		= otpt.mn;
otpt.std 		= otpt.mn;
otpt.rms 		= otpt.mn;

for k = 1:length(id)
	maxid 				= find(agg(:,id(k))==m(id(k)));
	temp 				= ts(:,id(k)+clms(3)-1);
	maxid 				= (maxid - m(id(k)) + 1):(maxid  - mod(m(id(k)),12));
	otpt.mn(2,id(k)) 	= mean(temp(maxid));
	otpt.sum(2,id(k)) 	= sum(temp(maxid));
	otpt.std(2,id(k)) 	= sqrt(mean((temp(maxid)-mean(temp(maxid))).^2));
	otpt.rms(2,id(k))   = sqrt(mean(temp(maxid).^2));
	cont(:,id(k)) 		= 0;
	cont(maxid,id(k)) 	= 1;
end

cont 	= [0 0 cinfo(clms(3):end); ts(:,clms(1:2)) cont];
agg 	= [0 0 cinfo(clms(3):end); ts(:,clms(1:2)) agg];
