function dy = central_difference(y,x)
%%-----------------------------------------------------------------------
% Function central_difference to calculate the derivatives dy/dx
% x and y are vector

Tflag=0;
if size(y,1)==1    % Treat row vector as a column vector
   y=y.';
   Tflag=1;
end;

% Forward difference at left end
dy(1)=(y(2)-y(1))/(x(2)-x(1));

% Backward difference at right end
dy(length(y))=(y(end)-y(end-1))/(x(end)-x(end-1));

% Central Difference in interior
h=diff(x); 
h_i=h(1:end-1); 
h_ip1=h(2:end);
dy(2:end-1)=(-(h_ip1./h_i).*y(1:end-2)+(h_i./h_ip1).*y(3:end))./(h_i+h_ip1)+(1./h_i-1./h_ip1).*y(2:end-1);

% end function central_difference
%%-----------------------------------------------------------------------