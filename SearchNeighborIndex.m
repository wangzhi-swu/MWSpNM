function  par  =  SearchNeighborIndex( par )
% This Function Precompute the all the patch indexes in the Searching window
% -NeighborIndex is the array of neighbor patch indexes for each keypatch
% -NumIndex is array of the effective neighbor patch numbers for each keypatch
% -SelfIndex is the index of keypatches in the total patch index array
par.maxr  =  par.h - par.ps + 1;      % patch 行数
par.maxc  =  par.w - par.ps + 1;      % patch 列数
r         =  1:par.step:par.maxr;     % 1,6,11,16,…,(长度为 maxr/step +1 )
par.r     =  [r r(end) + 1:par.maxr]; % 1,6,11,16,…,余下的以1为步长
c         =  1:par.step:par.maxc;     
par.c     =  [c c(end) + 1:par.maxc];
par.lenr  =  length(par.r);           % key patch 行数
par.lenc  =  length(par.c);           % key patch 列数
par.ps2   =  par.ps^2;     % patch面积
par.ps2ch =  par.ps2 * par.ch; % patch面积 * 3通道

% Total number of patches in the test image
par.maxrc   =  par.maxr * par.maxc;   % patch 总数
% Total number of seed patches being processed (seed patch ~ key patches)
par.lenrc   =  par.lenr * par.lenc;   % key patch 总数
% index of each patch in image
par.Index = (1:par.maxrc);
par.Index = reshape(par.Index, par.maxr, par.maxc); % 所有patch列优先排列
% preset variables for all the patch indexes in the Searching window
par.NeighborIndex =   int32(zeros(4 * par.win^2, par.lenrc));
par.NumIndex      =   int32(zeros(1, par.lenrc));
par.SelfIndex     =   int32(zeros(1, par.lenrc));

for  i  =  1 : par.lenr
    for  j  =  1 : par.lenc % for each key patch
        row = par.r(i);  % 当前key patch 所在行(原始图像)
        col = par.c(j);  % 当前key patch 所在列(原始图像)
        off = (col-1) * par.maxr + row; % 当前key patch 的索引(all patch)
        off1 = (j-1) * par.lenr + i;    % 当前key patch 的索引(key patch)
        
       % the range indexes of the window for searching the similar patches
        rmin    =   max( row - par.win, 1 );
        rmax    =   min( row + par.win, par.maxr );
        cmin    =   max( col - par.win, 1 );
        cmax    =   min( col + par.win, par.maxc );
        
        idx     =   par.Index(rmin:rmax, cmin:cmax); % sw 覆盖的all patch
        idx     =   idx(:);  % 拉伸成列向量
        
        par.NumIndex(off1)  =  length(idx);
                    % NumIndex(i) := 第i号key patch 的search window 面积(i是key patch下标)
        par.NeighborIndex(1:par.NumIndex(off1), off1)  =  idx;
                    % NeighborIndex(1:NumIndex(i), i) := 第i号key patch的 sw
                    % 究竟覆盖了哪些 all patch (all patch下标)
        par.SelfIndex(off1) = off;
                    % SelfIndex(i) := 第i号key patch 在 all patch 的索引
    end
end
