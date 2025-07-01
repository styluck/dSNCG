% compute purity and entropy of given classification results
function [pur, ent] = comp_pe(X,true_ans,eps)
[n, k] = size(X);
tes = zeros(k,k);
for i=1:n
    for j=1:k
        if(X(i,j)>eps)
            tes(true_ans(i),j) = tes(true_ans(i),j)+1;
        end
    end
end
pur = sum(max(tes))/n;
tes_log = log2(tes./sum(tes));
tes2 = tes.*tes_log;
for i=1:k
    for j=1:k
        if(tes(i,j)==0)
            tes2(i,j) = 0;
        end
    end
end
ent = (sum(sum(tes2)))/(-n*log2(k));
end