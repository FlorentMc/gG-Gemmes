function [ x, y ] = euler ( f, xRange, yInitial, numSteps )
% [ x, y ] = euler ( f, xRange, yInitial, numSteps ) uses
% Euler's explicit method to solve a system of first-order ODEs
% y'=f(x,y).
% f = name of an m-file with signature
%    yprime = f(x,y)
% to compute the right side of the ODE as a row vector
% xRange = [x1,x2] where the solution is sought on x1<=x<=x2
% yInitial = column vector of initial values for y at x1
% numSteps = number of equally-sized steps to take from x1 to x2
% x = row vector of values of x
% y = matrix whose n-th row is the approximate solution at x(n).

% M. M. Sussman
% $Id: euler.m,v 1.2 2009-02-04 14:30:11 mike Exp $

x=zeros(numSteps+1,1);
x(1) = xRange(1);
h = ( xRange(2) - xRange(1) ) / numSteps;
y(1,:) = transpose(yInitial);
for k = 1 : numSteps
  x(k+1) = x(k) + h;
  y(k+1,:) = y(k,:) + h * transpose(feval( f, x(k), transpose(y(k,:)) ));
end
