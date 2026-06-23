function Dist = BrayCurtisDistance(Tot_Mat_abund)
% X: N x p (N samples, sum of replicate split into two group 1 or 2 (SynCom20 or SynCom21), p species), nonnegative = identical, 1 completly different
N = size(Tot_Mat_abund, 1);
Dist = zeros(N, N);
for i = 1:N
    xi = Tot_Mat_abund(i,:);
    for j = i+1:N
        xj = Tot_Mat_abund(j, :);
        numerator = sum(abs(xi - xj));
        denominator = sum(xi + xj);
        if denominator == 0
            d = 0; % both all-zero
        else
            d = numerator/denominator;
        end
        Dist(i,j) = d;
        Dist(j,i) = d;
    end
end
end