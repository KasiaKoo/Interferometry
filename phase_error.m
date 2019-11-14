function pe = phase_error(alpha, n, epsilon)
    if n == 0
        pe = alpha;
    elseif n == 1
        pe = alpha*(1+epsilon);
    elseif n == 2
        pe = alpha.*(1+epsilon*alpha);
    elseif n == -1
        pe = alpha.*(1+epsilon.*alpha - epsilon);
    elseif n == 'v'
        pe = alpha.*(1 + rand*alpha);
    end
    
end