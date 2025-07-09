function [M] = blockinv(a,d,c,b,inopt,outopt)

% BLOCKINV computes the block inverse of a matrix that can be
% compartmentalised into four blocks [A B;C D]
%
% M = blockinv(a,d,c,b)
% Calculates the block inverse of a non-symmetric matrix, whose 
% blocks are [a, b; c, d]
%
% M = blockinv(a,d,c)
% calculates the block inverse of a symmetric matrix, whose lower 
% triangular block alone is provided apart from the diagonal blocks.
%
% M = blockinv(a,d,b,'upper')
% The option 'upper' specifies upper triangular block is provided 
% instead of lower triangular block in a symmetric matrix.
%
% M = blockinv(...,'struct')
% The output is provided as a MATLAB structure instead of a full 
% matrix. If only three blocks were provided, only three blocks of 
% the symmetric matrix are provided in the structure. The output
% off-diagonal block is always the lower triangular block of the 
% inverse.
%
% Possible combinations:
% M = blockinv(a,d,c,b)
% M = blockinv(a,d,c,b,'struct')
% M = blockinv(a,d,c)
% M = blockinv(a,d,c,'struct')
% M = blockinv(a,d,b,'upper')
% M = blockinv(a,d,b,'upper','struct')
%
%--------------------------------------------------------------------------

% Created on: 4 February 2008, Stuttgart
% Author: Balaji Devaraju
%--------------------------------------------------------------------------

if nargin == 5
    if strcmp('struct',inopt)
        outopt = 'struct';
        if ischar(b)
            if strcmp('upper',b)
                inopt = 'upper';
                c = c';
                b = [];
                outopt = [];
            elseif strcmp('struct',b)
                outopt = 'struct';
                b = [];
                inopt = [];
            end
            if ~isequal(a,a') || ~isequal(d,d')
                error('A or D is not symmetric')
            end
            if ~isequal(size(a,2),size(c,2)) || ~isequal(size(d,1),size(c,1))
                error('Size of C does not the match the sizes of A and/or D')
            end
            sym = 1;
        else
            if ~isequal(size(a,1),size(b,1)) || ~isequal(size(a,2),size(c,2)) || ~isequal(size(d,1),size(c,1)) || ~isequal(size(d,2),size(b,2))
                error('Dimensions of the matrices do not match')
            end
            sym = 0;
            inopt = [];
        end
    else
        error('Unknown string provided for the output option')
    end
elseif nargin == 4
    if ischar(b)
        if strcmp('upper',b)
            inopt = 'upper';
            c = c';
            b = [];
            outopt = [];
        elseif strcmp('struct',b)
            outopt = 'struct';
            b = [];
            inopt = [];
        end
        if ~isequal(a,a') || ~isequal(d,d')
            error('A or D is not symmetric')
        end
        if ~isequal(size(a,2),size(c,2)) || ~isequal(size(d,1),size(c,1))
            error('Size of C does not the match the sizes of A and/or D')
        end
        sym = 1;
    else
        if ~isequal(size(a,1),size(b,1)) || ~isequal(size(a,2),size(c,2)) || ~isequal(size(d,1),size(c,1)) || ~isequal(size(d,2),size(b,2))
            error('Dimensions of the matrices do not match')
        end
        sym = 0;
        inopt = [];
        outopt = [];
    end
elseif nargin == 3
    if ~isequal(size(a,2),size(c,2)) || ~isequal(size(d,1),size(c,1))
        error('Size of C does not the match the sizes of A and/or D')
    end
    sym = 1;
    inopt = [];
    outopt = [];
end

if sym==0
    d = d\eye(size(d));
    if strcmp('struct',outopt)
        M = struct('NW',[],'NE',[],'SW',[],'SE',[]);
        b = b*d;
        M.NW = (a - b*c);
        M.NW = M.NW\eye(size(M.NW));
        c = d*c;
        M.NE = -M.NW*b;
        M.SW = -c*M.NW;
        M.SE = d - M.SW*b;
    elseif isempty(outopt)
        b = b*d;
        Q = a - b*c;
        Q = Q\eye(size(Q));
        c = d*c;
        E = -Q*b;
        W = -c*Q;
        S = d - W*b;
        M = [Q E; W S];
    end
elseif sym==1
    d = d\eye(size(d));
    if strcmp('struct',outopt)
        M = struct('NW',[],'SW',[],'SE',[]);
        M.NW = d*c;
        M.NW = (a - c'*M.NW);
        M.NW = M.NW\eye(size(M.NW));
        c = d*c;
        M.SW = -c*M.NW;
        M.SE = d - M.SW*c';
    elseif isempty(outopt)
        Q = d*c;
        Q = a - c'*Q;
        Q = Q\eye(size(Q));
        c = d*c;
        W = -c*Q;
        M = [Q W'; W (d - W*c')];
    end
end
