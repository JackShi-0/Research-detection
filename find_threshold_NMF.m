function threshold = find_threshold_NMF(data,Pfa,bin_n,RunNum)
%Pfa = 1e-4;
%bin_n = 1e3;
%RunNum = 1e6;
nPfa = length(Pfa);
Pfa_n = round(Pfa*RunNum);
[n_count,bin_c] = hist(data,bin_n);
sum_mt = triu(ones(bin_n));
sum_n = sum_mt*n_count';
ind = zeros(nPfa,1);
for kk = 1:nPfa
    ind_n = find(sum_n < Pfa_n(kk), 1, 'first');

    if ind_n > 1
        [tmp, ind_sel] = min([abs(sum_n(ind_n -1) - Pfa_n(kk)), abs(sum_n(ind_n) - Pfa_n(kk))]);
        if ind_sel == 1
            ind(kk) = ind_n -1;
        else
            ind(kk) = ind_n;
        end
    else
        ind(kk) = ind_n;
    end
end

threshold = bin_c(ind);

end