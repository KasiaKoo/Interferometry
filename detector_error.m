function de = detector_error(I, n, epsilon)
    if n == 0
        de = I;
    elseif n == 1
        de = I.*(1+epsilon);
    elseif n == 2
        de = I.*(1+epsilon.*I);
    elseif n == 'g'
        de = I.^epsilon;
    elseif n=='p'
        de = I + sqrt(I);
    elseif n=='t'
        de = I.*(1 + rand);
    end
end