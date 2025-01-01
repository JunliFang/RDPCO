function g_inv=generalize_inv(sigma,D,E)
[Q_S_o, L_S_o] = eig(D*E*E'*D'+ D*diag(sigma)*D');
l_s_o = diag(L_S_o);
nnz_idx = find(abs(l_s_o)>1e-8);
l_s_o_inv = l_s_o;
l_s_o_inv(nnz_idx) = 1./l_s_o(nnz_idx);
g_inv = real(Q_S_o*diag(l_s_o_inv)*Q_S_o');
end