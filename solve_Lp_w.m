function   x   =  solve_Lp_w( y, W, p, J, lambda )

% tau is generalized thresholding vector
% lambda = 0.6; W \in (M*1)
% tau   =  (2*lambda*W.*(1-p)).^(1/(2-p)) + lambda*p*W.*(2*lambda*(1-p)*W).^((p-1)/(2-p));
tau   =  (2*W.*(1-p)).^(1/(2-p)) + p*W.*(2*(1-p)*W).^((p-1)/(2-p));

x     =   zeros( size(y) );

i0    =   find( abs(y)>tau );

if length(i0)>=1
    y0    =   y(i0);    
    t     =   abs(y0);
    W0    =   W(i0);
    for j = 1:J		% 论文中 K
        t =  abs(y0) - p*W0.*(t).^(p-1);
    end
    x(i0) =  sign(y0).*t;
end
return;

% 计算tau时多了lambda=0.6，源代码的lambda改名为W(向量)，其余一样