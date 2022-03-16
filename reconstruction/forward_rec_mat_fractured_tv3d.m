function v_arr=forward_rec_mat_fractured_tv3d(x_arr,zero_idx_list,n,m,l)
% v = A*w
% A = S*W were S is the sensing matrix and W is the wavelet sparsifier
% w is the vector in wavelet space
% x=W*w
N_had=n*m;

x_mat = reshape(x_arr,[m*n,l]);


if isa(x_arr,'gpuArray')
    %v_mat = gpuArray(fht_c(gather(x_mat)));
    
    temp_cube = reshape(x_mat,m,n,l);
    v_cube = fast_separable_hadamard_image_transform(temp_cube);
    v_mat = reshape(v_cube,m*n,l);
else
    v_mat = fht_c(x_mat);
end

v_mat=v_mat(zero_idx_list,:)*sqrt(N_had);

v_arr = v_mat(:);
