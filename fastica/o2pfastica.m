function [R] = o2pfastica(S)
  [A, ~] = fastica(S, 'verbose', 'off');
  % normalize columns before return
  R = A ./ repmat(sqrt(sum((A .^ 2),1)),size(A,1),1);
end
