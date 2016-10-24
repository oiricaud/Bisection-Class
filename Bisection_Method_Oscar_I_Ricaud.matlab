% Oscar I. Ricaud Math 4329: Numerical Analysis

% function bisect(a0, b0, ep,max_iterate)
% This is the bisection method for solving an equation f(x) = 0
% 
%The function f is defined below by the user. The function f is to be
% to be continuous on the interval [a0, b0], and it is to be of opposite
% signs at a0 and b0. The quantity ep is the rror tolerance. The parameter
% max_iterate is an upper limit on the number of iterates to be computed.
%
%
% This program gurantees ep as the error bound for the computed root
% provided: (1) the restrictions on the given function f and the inition
% [a0, b0] are satisfied; (2) ep is not too small when the machine epsilon
% is taken into account; and (3) the number of iterates computed is at most
% max_iterate. Only some of these conditions are checked in the program!
% 
% For the given function f(x), an example of calling sequence might be the
% following: 
%           [root, error_bound] = bisect(1,1.5,1.0E-6,10)
%
% The following is printed for the each iteration the values of count, a,b
% ,c ,f(c), (b-a)/2
% with c the ucrrent iterate and (b-a)/2 the error bound for c.
% The variable count is the index of the current iterate. Tap the carriage
% return to continue with the iteration. 
function [root, error_bound]=Bisection_Method_Oscar_I_Ricaud(a0,b0, ep, max_iterate)

fprintf('\n it_count,    a,     b,      c,     f(c), b-c(error), rel_error \n')

if(a0 >= b0)
    disp('a0 < b0 is not true. Stop!')
    return
end

format short e 
a = a0; b = b0;
fa = f(a); fb = f(b);

if sign(fa)*sign(fb) >0
    disp('f(a0) and f(b0) are of the same sign. Stop!')
    return
end

c = (a+b)/2;
old = c;
it_count = 0;

while b-c > ep && it_count < max_iterate
    it_count = it_count + 1;
    fc = f(c);
    iteration = [it_count a b c fc b-c (b-c/c)]
    if sign(fb)* sign(fc) <= 0
        a = c;
        fa = fc;
    else
        b = c;
        fb = fc;
    end 
    c = (a+b)/2;
    
end
format long 
root = c;
format short e
error_bound = b-c;
disp(iteration);



function value = f(x)
value = x.^3-3*x.^2+3*x-1;