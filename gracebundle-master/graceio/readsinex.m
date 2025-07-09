function N = readsinex(fname, var_)


% N = readsinex(fname, var_)
%

narginchk(2, 2)

if exist(fname, 'file')==2
    fid       = fopen(filename, 'r');
    while 1
        one_line = fgets(fid);
        if ~ischar(one_line), break, end
        hlines  = hlines + 1;
        if isempty(one_line)
            continue
        else
            keyword = textscan(one_line,'%s',1);
            if isempty(keyword{1})
                continue
            else
                keyword = lower(keyword{1}{1});
            end
        end
        
        
    end
else
    error('File does not exist. Check input')
end
