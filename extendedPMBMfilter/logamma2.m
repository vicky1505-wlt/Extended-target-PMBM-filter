function y = logamma2(x)
    y = 1/2*log(pi)+gammaln(x)+gammaln(x-1/2);
end