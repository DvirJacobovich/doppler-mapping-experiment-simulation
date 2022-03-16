function out_mat = ifht_c(in_mat)
% the function operates on matrix columns

out_mat = hadamardc(in_mat);

% out_mat = zeros(size(in_mat));
% 
% for i=1:size(in_mat,2)
%     out_mat(:,i) = hadamardc(in_mat(:,i));
% end
