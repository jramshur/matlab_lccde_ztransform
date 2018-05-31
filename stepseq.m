function [x,n]=stepseq(n0,n1,n2)
%Generate  x(n)=u(n-n0);n1<=n<=n2
%[x,n]=stepseq(n0,n1,n2)

    STEP = 1;

    n = n1:STEP:n2;
    x = n > n0;
end