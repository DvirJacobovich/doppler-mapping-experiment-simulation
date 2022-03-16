function out_mat = fht_c(in_mat)
% the function operates on matrix columns

N = size(in_mat,1);

out_mat = hadamardc(in_mat)/N;

% out_mat = zeros(size(in_mat));
% 
% for i=1:size(in_mat,2)
%     out_mat(:,i) = hadamardc(in_mat(:,i))/N;
% end
