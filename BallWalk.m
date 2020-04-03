X = [0, 0, 0];
calls = 0;
pos = [X]; 
out_pos = [X];
function [retval] = oracle(x)
  retval = (x <= 1) && sum(x < 0) == 0;
endfunction

delta = .989;
figure(); 
for j =1:10000;
 for i = 1:100;
  tmpx = X + unifrnd(0, delta, 1, 3);
  if(oracle(tmpx))
   calls++;
   X = tmpx;
   pos = [pos; X];
   plot3(pos(:, 1), pos(:,2), pos(:,3), 'b-');
   hold on; 
   %plot(X, 'b-'); hold on;
  else
   break;
  end
 end
 X = [0, 0, 0];
 pos = X;
 out_pos = X;
end
pos;
out_pos;
calls
%{
numjumps = 200; %number of steps
R = 0.5; %probability of step to right
x = rand (numjumps, 2);
step = 2 * (x >= R) - 1;
position = [0, 0; cumsum(step)]
plot (position(:, 1), position (:, 2))
%}