%%
function [C, time, H, sigma, B, omiga, iter]=DEpD_Homo(Y1,X1,  Y2,X2,  anneling, P_thr, gamma0, Constraints, Bilateral, DscDis, h)
%% 换X点

global Sift_Ratio  
            
[N, D] = size(X2);

%% initialization
AffineCPD_Initial;    
if ~exist('Constraints','var');      Constraints = 'L2-Length';     end;
if ~exist('gamma0','var');      gamma0 = 10;     end;
if exist('DscDis', 'var')   
    if ~isempty(DscDis)
        alpha_l = 2e-4;      alpha_p = 0.5; 
        Omega = ComputeOmega(DscDis, alpha_l, alpha_p);     %hist(Omega.lines)
    end
end
if ~exist('Omega','var');       Omega = [opt.outliers, opt.outliers];   end;
if ~exist('h','var');           h = eye(3);  	h = h(1:8)';   end; 

%% default mean and scaling
normal.xd=0; normal.yd=0; normal.xscale=1; normal.yscale=1;


%% transfer to multi-points
X = [];         Y = []; 
for i = 1:4:D-3
    if all(abs(X2(:, 4))<10)
        [X2temp, Y2temp]=LAF_to_pt3x3(X2(:,i:i+3), Y2(:,i:i+3)); % 
    else
        X2temp = X2;        Y2temp = Y2;
    end
    X = [X,  X2temp(:, 1:2),  X2temp(:, 3:4)];              Y = [Y,   Y2temp(:, 1:2),  Y2temp(:, 3:4)];
end
% if ~isempty(X1)   
%     [Xtemp,Ytemp,normal]=Homo_normalize([X(:,1:2); X(:,3:4); X1(:,1:2)],[Y(:,1:2); Y(:,3:4); Y1(:,1:2)],1); 
%     X= [Xtemp(1:N, 1:2), Xtemp(N+1:2*N, 1:2)];
%     Y= [Ytemp(1:N, 1:2), Ytemp(N+1:2*N, 1:2)];
%     X1(:,1:2)= Xtemp(2*N+1:end, 1:2);
%     Y1(:,1:2)= Ytemp(2*N+1:end, 1:2);
%     X1(:,3) = X1(:,3)/normal.xscale;  
%     Y1(:,3) = Y1(:,3)/normal.yscale;
% end



%% 归一化
if isempty(X1)   
%     [X,Y,normal]=Homo_normalize(X,Y,1);  % flag: 按2列归一化，但保留原size
    [X,Y,normal]=AffineCPD_normalize(X,Y,1);  % flag: 按2列归一化，但保留原size
else
    [X,Y,normal]=AffineCPD_normalize(X,Y,1);  % flag: 按2列归一化，但保留原size
    [X1(:,1:2),Y1(:,1:2),~]=AffineCPD_normalize(X1(:,1:2),Y1(:,1:2),1,normal); 
    X1(:,3) = X1(:,3)/normal.xscale;  
    Y1(:,3) = Y1(:,3)/normal.yscale;
end


%% 归一化之后重新计算尺度和方向
Linelength = ones(size(X,1), 1);       theta = ones(size(X,1), 1); 
for i = 1:N
    for j = 1:4:D-3
        [theta(i,j), Linelength(i,j)] = cart2pol(X(i, j+2)-X(i, j),  X(i, j+3)-X(i, j+1));
    end
end
s = 1;          
S = getS(s, theta);
S.LinelengthX = Linelength;
S.LinelengthX2 = Linelength.^2;

for i = 1:N
    for j = 1:4:D-3
        [theta(i,j), Linelength(i,j)] = cart2pol(Y(i, j+2)-Y(i, j),  Y(i, j+3)-Y(i, j+1));
    end
end
S.LinelengthY = Linelength;
S.LinelengthY2 = Linelength.^2;

%% 找邻近点0
if ~isempty(X1)  
    K = 2; % 最近5个点中，距离最小的
    if strcmp( Constraints, 'L2-Clear-ptsScale')
        K = 1;
    end
    [Sel_Pts_X, Sel_idx] = findClosestPoint2(X, X1, K);
    % 选 重投影距离（点到直线）的最小
    [Reprojected_error, Relative_Reprojected_error_horizentoal, Relative_Reprojected_error_Vertical, X1_sel, Y1_sel, Sel_Pts_idx] = find_Reprojection_Closest(X1, Y1, X, Y,  S, Sel_idx);
    % 选 旋转角度最接近 的
%     [Reprojected_error, Reprojected_error_horizentoal, Reprojected_error_Vertical, X1_sel, Y1_sel, Sel_Pts_idx] = find_Frame_Closest(X1, Y1, X, Y,  S, Sel_idx);
    X1 = X1(:, 1:2);
    Y1 = Y1(:, 1:2);
else
    Relative_Reprojected_error_horizentoal = []; Reprojected_error = []; Relative_Reprojected_error_Vertical = []; X1=[];Y1=[];
end


%% 自适应初始化H
if exist('DscDis', 'var')    
    [h] = InitializationViaDsc(X, Y, DscDis.line, S, 10);
end

%%
tic
[C,  h,  sigma, B, omiga, iter] = Matching(X1, Y1, X, Y, anneling, opt.max_it, opt.tol, opt.viz, Omega, S, h, gamma0, P_thr, Constraints, Bilateral, Reprojected_error, Relative_Reprojected_error_horizentoal, Relative_Reprojected_error_Vertical);
time=toc;

%% 因为采用x-Rb 的计算方法直接对x进行漂移，而x的尺度归一化为normal.xscale，因此
B = B * normal.xscale;
% B(1:2:end, :) = B(1:2:end, :) * normal.xscale;
% B(2:2:end, :) = B(2:2:end, :) * normal.yscale;
% B=B*normal.xscale/normal.yscale;        t=t*normal.xscale-B*normal.yd(1:2)'+normal.xd(1:2)';  
H = reshape([h; 1], 3,3)';
T1 = diag([normal.xscale,normal.xscale,1]);    T1(1:2,3) = normal.xd(1:2)';
T2 = diag([normal.yscale,normal.yscale,1]);    T2(1:2,3) = normal.yd(1:2)';
H = T1 * H / T2; %T2 / (T1 * H); %inv(T1 * H / T2);
H = H/H(3,3);
H = H([2,1,3],[2,1,3]);
save normal normal
end
        
%% 算法核心程序
function [C, h, sigma, B, omega, iter] = Matching(X1, Y1, X, Y, anneling, max_it,  tol,  viz, omega,  S, h, gamma0, P_thr, Constraints, Bilateral, Rep_err, Rep_err_h, Rep_err_v)

[Nl, D] = size(X);             [Np, Dp] = size(X1);           

%% 判断是否有特征点 辅助
if exist('Rep_err', 'var') && ~isempty(Rep_err)
    SigGammaError = 0;
else
    Rep_err = zeros(Nl, 4);
    SigGammaError = 0;
end

%% R 矩阵
R = [];
for j = 1 : Nl
    temp = blkdiag(S.St(:,:,j), S.St(:,:,j), eye(D-4));
    temp = temp(:, [1,3]);
    R = blkdiag(R, temp);
end

%% 提前添加局部估计的偏移
gammaMax = 1e4;
switch Constraints
    case '0'
        gamma = 0;
    case 'L2'
        gamma = mean(gamma0(1)./S.LinelengthX2); 
    % 根据"线段长度"
    case 'L2-Length'
        gamma = gamma0(1)./S.LinelengthX2; 
%         gamma = gamma0(1)./max(S.LinelengthX2, S.LinelengthY2); 
        gamma = gamma*mean(gamma)/median(gamma);
    % 根据"线段长度差"
    case 'L2-Ratio'
        S.ratio = abs(S.LinelengthY2 - S.LinelengthX2);
        gamma = gamma0(1)./S.ratio;
        gamma = gamma*mean(gamma0(1)./S.LinelengthX2)/median(gamma);
    % 根据"重投影" 误差 【无特征点】
    case 'L2-Clear-ptsScale'
        X1 = [];     Y1 = [];       Np = 0;     
        if isempty(Rep_err_h)
            gamma = gamma0(1)./max(S.LinelengthX2, S.LinelengthY2); 
        else
            Estimated_Frac_Lenth2 = abs(Rep_err_h');
            gamma = gamma0(1)  ./  (Estimated_Frac_Lenth2(:)+eps)/5;
%             Estimated_Ver_Dis = kron(sum(Rep_err_v.^2, 2), [1,1]');
%             gamma = Estimated_Ver_Dis  .*  gamma0(1)  ./  (Estimated_Frac_Lenth2(:)+eps);
%             gamma = Estimated_Ver_Dis.^4  .*  gamma0(1)  ./  (Estimated_Frac_Lenth2(:).*kron(S.LinelengthY2, [1,1]'));
        end
    case 'L2-pts'
        gamma = gamma0(1)./max(S.LinelengthX2, S.LinelengthY2); 
    case 'L2-ptsScale'
        if isempty(Rep_err_h)
            gamma = gamma0(1)./max(S.LinelengthX2, S.LinelengthY2); 
        else
            Estimated_Frac_Lenth2 = abs(Rep_err_h');
            gamma = gamma0(1)  ./  (Estimated_Frac_Lenth2(:)+eps)/5;
%             Estimated_Ver_Dis = kron(sum(Rep_err_v.^2, 2), [1,1]');
%             gamma = Estimated_Ver_Dis  .*  gamma0(1)  ./  (Estimated_Frac_Lenth2(:)+eps);
%             gamma = Estimated_Ver_Dis.^4  .*  gamma0(1)  ./  (Estimated_Frac_Lenth2(:).*kron(S.LinelengthY2, [1,1]'));
        end
end
gamma = min(gamma, gammaMax);



%%
S1 = X1;               T1 = Y1;             S2 = X;             T2 = Y; 

%%
Zl = zeros(D*Nl, 8);           Zp = zeros(2*Np, 8);
for i = 1 : Nl
    % line
    Z2_1_1 = [kron(eye(2), [Y(i,1:2),  1]),        -X(i,1:2)'*Y(i,1:2)];
    Z2_1_2 = [kron(eye(2), [Y(i,3:4),  1]),        -X(i,3:4)'*Y(i,3:4)];
    temp = [Z2_1_1; Z2_1_2];
%     if PtNum==3
%         disp('还没有实现');% 第三点
%     end
    Zl(D*(i-1)+1: D*i, :) = temp; 
end
Zlh = Zl * h;  %8M, 1
X = X';      %X=X(:);% 4, 2*N

for i = 1 : Np
    % pts
    Z2_1 = [kron(eye(2), [Y1(i,:),  1]),        X1(i,1:2)'*Y1(i,:)];
    Zp(2*(i-1)+1: 2*i, :) = Z2_1;
end
Zph = Zp * h;
X1 = X1';

%% 自适应确定初始sigma   
sigma = (trace(X*X')   -   2*X(:)'*Zlh   +   Zlh'*Zlh   +   ...
         trace(X1*X1')   -   2*X1(:)'*Zph   +   Zph'*Zph) ...
         / (D*Nl + Dp*Np);
% sigma = 1;


saveGIF = 0;               
if saveGIF; 
    filename = 'results1';
%     global filename
end
%% 换X.点
iter=0; ntol=tol+10; L=1;   L_old = 0;  sigma_save = zeros(max_it, 1);  flag = 0;
ntol_thr = 1e-5;    sigma_thr = 1e-5;
u = (max(X(:))-min(X(:)))/2 ;
while ((iter<max_it) && (sigma > sigma_thr)  && (ntol > ntol_thr))
    E(iter+1) = L;
%     a_Save(:, iter+1) = h;
    sigma_save(iter+1) = sigma;
    L_old=L;

    [Npp, Npl, Pp, Pl, Ep, El] = AffineCPD_Pts_Line_Omega(  S1,T1,  S2,T2,  sigma, omega,u);
%         [Npp, Npl, Pp, Pl, Ep, El] = AffineCPD_Pts_Line(  S1,T1,  S2,T2,  sigma, omiga);
     
    L = Ep + El;
    
    Pl = max(Pl, 1e-100);           Pp = max(Pp, 1e-100); 
    Pl4 = kron(Pl, eye(D));         Pp2 = kron(Pp, eye(2));       
    
    if isnan(L);
        pause(1)
    end
    ntol = abs((L-L_old)/L);
%     disp([' CPD Affine ' st ' : dL= ' num2str(ntol) ', iter= ' num2str(iter) ' sigma2= ' num2str(sigma)]);

    %% 估计b 
    Xvec = X(:);
    if length(gamma)==Nl
        gamma = kron(gamma, [1,1]');
    end
    coe = diag(1 ./ (1 + 2*sigma*gamma));
    if ~Bilateral
        % (x 单边漂移)
        b = coe * R'  * ( Xvec - Zlh + 2*sigma*SigGammaError(:));%
        Rb = R*b;
        X_B = X - reshape(Rb, D, []);
        X_Z = Zl;
    else
        % (x-y 双边漂移)
        bx = coe * R' * Xvec;%
        bz = coe * R' * Zl;%
        b = bx - bz*h;
        Rb = R*b;
        X_B = X - reshape(R*bx, D, []);
        X_Z = Zl - R*bz;
    end

    
    % Solve for h
    h1 = X_Z' * Pl4 * X_Z   +   Zp' * Pp2 * Zp;
    h2 = X_B(:)' * Pl4 * X_Z   +   X1(:)' * Pp2 * Zp;
    h = h1 \ h2'; 

    % update Y postions 
    Zlh = X_Z * h;                                  
    T2 =  reshape(Zlh, D, Nl)';                   
    Zph = Zp * h;   %2, 4M
    T1 =  reshape(Zph, Dp, Np)';                     
    
    % update individual translation                                  
    S2 = reshape(X_B, D, Nl)';
    

   %% 更新sigma
    sigma_trial = (trace(X_B * Pl  * X_B')   -   2*X_B(:)' * Pl4 * Zlh   +   Zlh' * Pl4 * Zlh   +   ...
             trace(X1*Pp*X1')   -   2*X1(:)'*Pp2*Zph   +   Zph'*Pp2*Zph) ...
             / (D*Npl + Dp*Npp);
   %% 确定性退火
   if anneling
       ntol_thr = -inf;
       sigma_thr = 5e-4;
       if iter>=1
    %         iter_Flag = abs((sigma_save(iter) - sigma_save(iter+1)) / sigma_save(iter+1));
            if (ntol < 1e-5) || flag == 1    % && sigma_trial>1e-3)
                flag = 1;
                alpha = 0.5;
                sigma_trial = alpha*sigma;
            end
       end
   end
    sigma = sigma_trial;
%     omega(1) = 1 - Npp/Np;
%     omega(2) = 1 - Npl/Nl;
%     omega = max(omega, 0.02);
%     omega = min(omega, 0.98);
    
     iter=iter+1;
    if viz, cpd_plot_iter(reshape(S2, 2, [])', reshape(T2, 2, [])'); end;
    
    B = b;
    b_save(:,iter) = Rb;
    if viz, 
        figure(1)
        set(gcf,'units','normalized','outerposition',[0 0 1 1])
        cpd_plot_iter_our(S2, T2, iter); 
        frame = getframe(1);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        if iter==1
            imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
        else
            imwrite(imind,cm,filename,'gif','WriteMode','append');
        end
        pause(0.1)
    end;
%     close 
end

save b_save b_save  % 找出对的，分析误差e，画图显示
[~, ~, Pp, Pl, ~, ~] = AffineCPD_Pts_Line_Omega(  S1,T1,  S2,T2,  sigma, omega,u);
Pp = diag(Pp);
Pl = diag(Pl);
Index_p = find(Pp>P_thr); %行最大
Index_l = find(Pl>P_thr); %行最大
% C.points = [Index_p,Index_p,Pp(Index_p)];
% C.line = [Index_l,Index_l,Pl(Index_l)];
C.points = [Index_p,Index_p];
C.line = [Index_l,Index_l];

% % 检查对正确点的约束是否有效
%     Rep_Error = Reprojected_error(Index_l,:);
%     Error = X' - T2;    Error = Error(Index_l,:);

end

%% 计算S
function [Linelength, theta] = verifyParameter(Linelength, theta, NL, NNN)
    if length(Linelength) == NNN(2)
        LinelengthTemp = Linelength;
        Linelength = zeros(NL(end), 1);
        Linelength(NL(1)+1:NL(3)) = [LinelengthTemp; LinelengthTemp];
        thetaTemp = theta;
        theta = zeros(NL(end), 1);
        theta(NL(1)+1:NL(3)) = [thetaTemp; thetaTemp];
    end
end



%% 计算XY123
function [X1,Y1,  X2_1,Y2_1,  X2_2,Y2_2,  X3_1,Y3_1,  X3_2,Y3_2,  X3_3,Y3_3] = getXY123(X, Y, NL)
    X1 = X(1:NL(1), :);     
    X2_1 = X(NL(1)+1:NL(2), :);  
    X2_2 = X(NL(2)+1:NL(3), :);
    X3_1 = X(NL(3)+1:NL(4), :); 
    X3_2 = X(NL(4)+1:NL(5), :); 
    X3_3 = X(NL(5)+1:NL(6), :);
    Y1 = Y(1:NL(1), :);     
    Y2_1 = Y(NL(1)+1:NL(2), :); 
    Y2_2 = Y(NL(2)+1:NL(3), :);     
    Y3_1 = Y(NL(3)+1:NL(4), :);
    Y3_2 = Y(NL(4)+1:NL(5), :);
    Y3_3 = Y(NL(5)+1:NL(6), :);
end
%% 计算S
function S = getS(s, theta)
    % 计算线段组数量
    [N, Line_Number] = size(theta);
    S.S = zeros(2, 2, N);
    S.invS = zeros(2, 2, N);
    for j = 1:Line_Number
        for i=1:N   
%             if i>N && i<=N 
                if numel(s) == 1;     s_i=s;    else;    s_i = s(i);        end;
                scaling = [s_i, 1];
                rot = [cos(theta(i)), sin(theta(i));  -sin(theta(i)), cos(theta(i))];
                S.S(:,:,i,j) = diag(scaling) * rot;
                S.St(:,:,i,j) = S.S(:,:,i)';
                S.StS(:,:,i,j) = S.St(:,:,i)*S.S(:,:,i);
                S.invS(:,:,i,j) = diag(1./scaling) * rot';
%             else
%                 S.S(:,:,i,j) = eye(2);
%                 S.St(:,:,i,j) = eye(2);
%                 S.StS(:,:,i,j) = eye(2);
%                 S.invS(:,:,i,j) = eye(2);
%             end 
        end
    end
end

% 选 重投影距离（点到直线）的最小
function [Reprojected_error, Relative_Reprojected_error_horizentoal, Relative_Reprojected_error_Vertical, X1_sel, Y1_sel, Sel_Pts_idx] = find_Reprojection_Closest(X1, Y1, X2, Y2, S, Sel_idx)

N = size(Y2, 1);    
K = size(Sel_idx, 2);       
Sel_Pts_idx = zeros(N, 1);
Reprojected_2pts = cell(N, K);

Reprojected_error = zeros(N, 4);
Relative_Reprojected_error_horizentoal = zeros(N, 2);
Relative_Reprojected_error_Vertical = zeros(N, 2);
Error = 1e5*ones(4,K);
Error_Vertical = 1e5*ones(2,K);
Error_Horizentoal = 1e5*ones(2,K);
Scaling = X1(:,3)./Y1(:,3);   Theta = X1(:,4) - Y1(:,4);

for i = 1 : N
    % 计算Y2 K次重投影后的端点
    for k = 1 : K
        idx = Sel_idx(i, k);
        ptsX = X1(idx, [1,2]);
        ptsY = Y1(idx, [1,2]);
        LineY_2pts = reshape(Y2(i,:), 2,2) - ptsY';
        LineX_2pts = reshape(X2(i,:), 2,2) - ptsX';
        Rotation = [cos(Theta(idx)), sin(Theta(idx)); -sin(Theta(idx)), cos(Theta(idx))];   % Rotation
        Reprojected_2pts{i, k} = Scaling(idx) .* (Rotation * LineY_2pts);
        Error_temp = LineX_2pts - Reprojected_2pts{i, k};                      % Error
        Rot = blkdiag(S.St(:,:,i), S.St(:,:,i));
        Rot = Rot(:, [2,4]);
        Error_Vertical(:,k) = Rot'*Error_temp(:);
        Rot = blkdiag(S.St(:,:,i), S.St(:,:,i));
        Rot = Rot(:, [1,3]);
        Error_Horizentoal(:,k) = Rot'*Error_temp(:);
        Error(:,k) = Error_temp(:);
    end
%     if i==113
%         pause(0.1)
%     end
    % 取最小值为最终y_H(y_r)，并输出
    [~, k_min] = min(sum(Error_Vertical.^2));
    Sel_Pts_idx(i) = Sel_idx(i, k_min);
    X1_sel = Y1(Sel_Pts_idx(i), :);
    Y1_sel = X1(Sel_Pts_idx(i), :);
%     Reprojected_pts(i, :) = Reprojected_2pts{i, k_min}(:)';
    Lenth_Reprojected_Y = norm(Reprojected_2pts{i, k}(:,2)-Reprojected_2pts{i, k}(:,1));
    Relative_Reprojected_error_horizentoal(i, :) = Error_Horizentoal(:,k_min)';% * Lenth_Reprojected_Y.^1.2;
    % 如果不重合，横向偏移远超过 线段本身长度，且是同一个方向偏移； 则需要去掉
    if prod(Error_Horizentoal(:,k_min))>0 && all(abs(Error_Horizentoal(:,k_min)) > 2*Lenth_Reprojected_Y)
        Relative_Reprojected_error_horizentoal(i, :) = 1e-5;
    end
    Relative_Reprojected_error_Vertical(i, :) = Error_Vertical(:,k_min)';
    % 如果垂直误差很大或其相对线长很大，则需要去掉
    Error_Length_ratio = norm(Error_Vertical(:,k_min)) / Lenth_Reprojected_Y;
    if Error_Length_ratio > 0.5 || norm(Error_Vertical(:,k_min))>0.1
%         Relative_Reprojected_error_Vertical(i, :) = 1e3;
        Relative_Reprojected_error_horizentoal(i, :) = 1e-5;
    end
    Reprojected_error(i, :) = Error(:,k_min)';
end

end

function [H_initial] = InitializationViaDsc(X2, Y2, DscDis, S, K)

H_initial=eye(3); H_initial=H_initial(:); H_initial=H_initial(1:8);

alpha = 0.05;
deta = 0.1;
[Sorted_DscDis, Dsc_Sel] = sort(DscDis);
if Sorted_DscDis(1)>30;     return;     end
Idx_sel_25 = find(Sorted_DscDis<25);
K = min(max(K, length(Idx_sel_25)), 40);
X2 = X2(Dsc_Sel(1:K), :);
Y2 = Y2(Dsc_Sel(1:K), :); 

N = size(Y2, 1);           
Dis = zeros(N, 3);
Reprojected_2pts = cell(N, K);
Error = 1e5*ones(4,K);
Error_Vertical = 1e5*ones(2,K);
Error_Horizentoal = 1e5*ones(2,K);
Scaling = X2(:,3)./Y2(:,3);   Theta = X2(:,4) - Y2(:,4);

for i = 1 : K
    % 计算Y2 K次重投影后的端点
    for k = 1 : K
        if k==i;    continue;       end;
        LineY_2pts = reshape(Y2(i,:), 2,2);
        LineX_2pts = reshape(X2(i,:), 2,2);
        Rotation = [cos(Theta(k)), sin(Theta(k)); -sin(Theta(k)), cos(Theta(k))];   % Rotation
        Reprojected_2pts{i, k} = Scaling(k) .* (Rotation * LineY_2pts);
        Error_temp = LineX_2pts - Reprojected_2pts{i, k};                      % Error
        Rot = blkdiag(S.St(:,:,i), S.St(:,:,i));
        Rot = Rot(:, [2,4]);
        Error_Vertical(:,k) = Rot'*Error_temp(:);
        Rot = blkdiag(S.St(:,:,i), S.St(:,:,i));
        Rot = Rot(:, [1,3]);
        Error_Horizentoal(:,k) = Rot'*Error_temp(:);
        Error(:,k) = Error_temp(:);
    end
    Dis(i, 1) = sum(sqrt(sum(Error_Vertical.^2))<deta);
    Dis(i, 2) = sum(exp(-alpha*sum(Error_Horizentoal.^2)));
    Dis(i, 3) = sum(sqrt(sum(Error.^2))<deta);
%     Dis(i, 1) = sum(exp(-alpha*sum(Error_Vertical.^2)));
%     Dis(i, 2) = sum(exp(-alpha*sum(Error_Horizentoal.^2)));
%     Dis(i, 3) = sum(exp(-alpha*sum(Error.^2)));

end
Dis_vertical_sorted = sort(Dis(:, 1), 'descend');
thr = max([2, Dis_vertical_sorted(6), multithresh(Dis_vertical_sorted)]);
idx_vertical = find(Dis(:, 1)>=thr);
if isempty(idx_vertical);    return;    end

[~, sel_min] = min(Dis(idx_vertical, 2));
Idx_min = idx_vertical(sel_min);
Rotation = [cos(Theta(Idx_min)), sin(Theta(Idx_min)); -sin(Theta(Idx_min)), cos(Theta(Idx_min))];   % Rotation
H_initial = eye(3);
H_initial(1:2, 3) = mean(LineX_2pts - LineY_2pts, 2);
H_initial(1:2, 1:2) = Scaling(Idx_min) .* Rotation;
H_initial(3, :) = [0,0,1];
H_initial=H_initial(:); H_initial=H_initial(1:8);
end

function [Selected_Points, ind] = findClosestPoint2(Lines_2Pts, Points, K)
    
N_lines = size(Lines_2Pts, 1);
M_points = size(Points, 1);

dis_son = 1e5*ones(M_points, 1);
Selected_Points = cell(N_lines, 1);
ind = zeros(N_lines, K);

for i = 1:N_lines

    for j = 1:M_points
        dis_son(j) = norm(Points(j, 1:2) - Lines_2Pts(i, 1:2))  +  norm(Points(j, 1:2) - Lines_2Pts(i, 3:4));
    end
    
    [~, ind_temp] = sort(dis_son);
    ind(i, :) = ind_temp(1:K);
    Selected_Points{i} = Points(ind(i, :), :);
end

end

function proj_point = ProjPoint( point,line_p )
x1 = line_p(1);
y1 = line_p(2);
x2 = line_p(3);
y2 = line_p(4);

x3 = point(1);
y3 = point(2);

yk = ((x3-x2)*(x1-x2)*(y1-y2) + y3*(y1-y2)^2 + y2*(x1-x2)^2) / (norm([x1-x2,y1-y2])^2);
xk = ((x1-x2)*x2*(y1-y2) + (x1-x2)*(x1-x2)*(yk-y2)) / ((x1-x2)*(y1-y2));


if x1 == x2
    xk = x1;
end

if y1 == y2
    xk = x3;
end

proj_point = [xk,yk];

end