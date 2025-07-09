% This function calls the mex compiler for all C code of the SHbundle.
% At the moment, this is only Legendre_mex.
%
% Please note, that you also need additionally an appropriate C compiler
% for your system architecture.

mex ./src/Legendre_mex.c ./src/legendre_od.c ./src/legendre_plm.c ./src/x_numbers.c