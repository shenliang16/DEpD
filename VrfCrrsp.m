function [Crr_GT,d1,d2, crsp, Lengthratio, C, Dis,bn] = VrfCrrsp(B_opt,  X1, X2, size1, size2, threshold, NNN, thr_Line,Im1, Im2, C_putative, LSM_GT)
bn = []; crsp=[]; Lengthratio=[]; C=[];   Dis = [];

%% 形变场
    if isfield(B_opt, 'dfield')  ||  (size(B_opt,1)>10 || size(B_opt, 3)>1)
        if isfield(B_opt, 'dfield')
            B_opt_new = inv(B_opt.tform)';
            T = B_opt_new([2,1,3],[2,1,3]); 
            N=size(X1,1);        X3=[X1(:,1:2),ones(N,1)];     
            XT=(T*X3')';         X1=XT(:,1:2)./XT(:,3);
            B_opt = B_opt.dfield;
        end

        h = fspecial('average', 5);
        B_opt(:,:,1) = imfilter(B_opt(:,:,1), h);
        B_opt(:,:,2) = imfilter(B_opt(:,:,2), h);

        moveX = interp2(B_opt(:,:,1), X2(:,2), X2(:,1));
        moveY = interp2(B_opt(:,:,2), X2(:,2), X2(:,1));
        T_X2 = X2(:, 1:2) - [moveX, moveY];

        d1 = X1(:,1:2)-T_X2;
        d2=sum(d1.^2, 2);
        Crr_GT=(sqrt(d2)<threshold);
        return;
    end



if exist('NNN','var')
    if isfield(NNN, 'flag')
        flag = NNN.flag;
    else
        flag=B_opt(end,end); 
    end
else
    flag=B_opt(end,end); 
end

%         F_gt=B_opt;    % threshold=max(0.0025,threshold);
%         [Crr_GT1,~,~] = DeterminInlier(F_gt,  X2_1, X2_2, size1, size2, threshold);
%                         [X4_1, X4_2]=LAF_to_pt3x3(X4_1(:,1:4), X4_2(:,1:4)); % 
%         [Crr_GT2,~,~] = DeterminInlierLine(F_gt,  X4_1, X4_2, size1, size2, threshold);
%         [Crr_GT2,d1,d2] = DeterminInlier(F_gt,  X6_1, X6_2, size1, size2, threshold);

if abs(flag+1)<1e-99
    mean=size1/2;

    % Make center of the image coordinates 0,0
    xd=X1(:,1)-mean(1); 
    yd=X1(:,2)-mean(2);

    % Calculate the Transformed coordinates
    Tlocalx = mean(1) + B_opt(1,1) * xd + B_opt(1,2) *yd + B_opt(1,3) * 1;
    Tlocaly = mean(2) + B_opt(2,1) * xd + B_opt(2,2) *yd + B_opt(2,3) * 1;
    d1=[]; 
    threshold=max(6,threshold);
    d2=sqrt(sum((X2(:,1:2)-[Tlocalx,Tlocaly]).^2,2));

    Crr_GT=(d2<threshold);
    
elseif abs(flag-1)<1e-99
%     if isfield(threshold, 'Exausted')
%         threshold = threshold.threshold;
%         
%     end
    T=B_opt([2,1,3],[2,1,3]);     %T=inv(T); T = T/T(3,3);
    threshold=max(5,threshold);
    [Crr_GT, d1, crsp, d2]=DetGTusingT(X1,X2,[],[],threshold,T);
    
%         threshold = 8*max(size1)/512;
%         T = B_opt;
%         [Crr_GT_trial,d1_trial, crsp_trial]=DetGTusingT(X1,X2,[],[],threshold,T);
%         if sum(Crr_GT_trial)/(sum(Crr_GT)+eps)>2
%             Crr_GT =Crr_GT_trial;
%             d1=d1_trial;
%             crsp=crsp_trial;
%         end

    % 线特征Ground Truth
    if exist('NNN','var')
        [Crr_GT, Lengthratio, C, Dis,bn] = get_LineGroundTruth(B_opt,  X1, X2, size1, size2, threshold,  Crr_GT,  NNN, thr_Line);
    end
elseif flag==0
    % 运行 exp1:  getGroundTruth 估计真值
    d2=[];          threshold=max(5,threshold);
    [Crr_GT, d1, crsp]=DetGTusingT(X1,X2,Im1,Im2,threshold);
    
    
else
    F_gt=B_opt;    threshold=min(0.003,threshold);
    if isempty(X1)
        Crr_GT = []; 
        d1 = [];
        d2 = [];
        return;
    end
    [Crr_GT,d1,d2] = DeterminInlier(F_gt,  X1, X2, size1, size2, threshold);
%     if exist('NNN','var')
%         Crr_GT2 = LSVFDR_verify (X1, X2, NNN);
%         Crr_GT = Crr_GT | Crr_GT2;
%     end
    Crr_GT2 = false(length(Crr_GT), 1);
    Crr_GT2_idx = find(Crr_GT);
    if exist('NNN','var') && length(Crr_GT2_idx)>9
        N_tested = sum(Crr_GT);
        X1_v = [X1(Crr_GT,1:2); NNN.X1];
        X2_v = [X2(Crr_GT,1:2); NNN.X2];
        [Clpm, ~]=LPM4test(X1_v, X2_v);
        %[Clpm, ~] = SIR_Method_4test(X1(Crr_GT,1:2),X2(Crr_GT,1:2), NNN.siftratio); %NNN.tested_idx
        Clpm = Clpm(Clpm<N_tested);
        Crr_GT2(Crr_GT2_idx(Clpm)) = true; 
        Crr_GT = Crr_GT2;
%         tested_idx = NNN.tested_idx;
%         Crr_GT(tested_idx) = Crr_GT2(tested_idx);
    end
end


%% 给定数据集的GroundTryth
if exist('LSM_GT', 'var')
    if ~isempty(LSM_GT)
        Crr_GT_41 = zeros(size(C_putative,1), 1);
        Crrsp1 = LSM_GT{1};
        Crrsp2 = LSM_GT{2};
        % 计算homography
        for i = 1:size(C_putative,1)
            Contain_flag1 = any(C_putative(i,1) == Crrsp1,2);
            if any(Contain_flag1)
                Index_in_1 = find(Contain_flag1);
                Contain_flag2 = any(C_putative(i,2) == Crrsp2(Index_in_1(1), :),2);
                if Contain_flag2
                    Crr_GT_41(i) = 1;
                end
            end
        end
        if size(C_putative,1) == size(Crr_GT,1) %% 如果只有直线特征
            Crr_GT = Crr_GT_41;
        else
            flag = find(NNN.Flag == 41)+1;  %[6, 41, 4, 2] 
            N = [0, cumsum(NNN.N)];
            if ~isempty(flag) 
                Crr_GT(N(flag-1)+1:N(flag)) = Crr_GT_41;
            end
        end
    end
end


Crr_GT = Crr_GT==1;

end

function [Crr_GT, Lengthratio, C, DisAll,bn] = get_LineGroundTruth(B_opt,  X1, X2, size1, size2, threshold,  Crr_GT,  NNN, thr_Line)
if isfield(NNN, 'Thr_times');   Thr_times = NNN.Thr_times;   else;    Thr_times = 20;    end;
% if isfield(NNN, '');        else; 0.93;
% 线特征Ground Truth
    T=B_opt([2,1,3],[2,1,3]);           Lengthratio_all=[];         threshold=min(5,threshold);         %DisAll = [];
    Lengthratio=-ones(size(X1,1), 3);   cosError=-ones(size(X1,1), size(X2,1));
    Lengthratio_temp = -ones(size(X1,1), size(X2,1));       Length_temp = zeros(size(X1,1), size(X2,1));
    flag = find(NNN.Flag == 41)+1;  %[6, 41, 4, 2] 
    N = [0, cumsum(NNN.N)];
    if ~isempty(flag)           
        X1_4 = X1(N(flag-1)+1:N(flag),1:4);               X2_4 = X2(N(flag-1)+1:N(flag),1:4);     
        
        if any(X1_4(:,4)>2*pi+1)
            X2temp = X1_4;
            Y2temp = X2_4;
        else
            [X2temp, Y2temp]=LAF_to_pt3x3(X1_4, X2_4); 
        end
        
        %% 自建数据集（自己估计的真值变换）
%         if max(size1)>2500
%             Y2temp = Y2temp(:,[2,1,4,3]);
%             X2temp = X2temp(:,[2,1,4,3]);
            threshold = Thr_times*max(size1)/512;
%         end

        bn = zeros(size(X1_4, 1), 2);    
        %% 计算 X 的变换
        for i = 1 : size(X2_4, 1)
            X1_32 = [reshape(X2temp(i, 1:4), 2,2); 1, 1];   
            TX1_32{i} = T * X1_32;     
            TX1_32{i} = TX1_32{i}./TX1_32{i}(3, :);
        end
        %% 计算 Y 的变换
        for j = 1 : size(X2_4, 1)
            scaling = [1, 1];
            [theta, ~] = cart2pol(Y2temp(j, 4)-Y2temp(j, 2),  Y2temp(j, 3)-Y2temp(j, 1));
%             theta = X2_4(j, 4);
            rot = [cos(theta), -sin(theta);  -sin(theta), cos(theta)];
            SR{j} = blkdiag(diag(scaling) * rot, 1);
            Y1_32{j} = [reshape(Y2temp(j, 1:4), 2,2); 1, 1];  
            SR_Y1_32_PreCal{j} = SR{j} * Y1_32{j}; 
        end
        
        %% 画图，双点对应变换后的空间关系图
        % idx = 8;scatter(TX1_32{idx}(1, :)', TX1_32{idx}(2, :)'); hold on;   scatter(Y2temp(idx, [1,3])', Y2temp(idx, [2,4])');  axis equal
        % idx = 8;scatter(TX1_32{idx}(1, :)', TX1_32{idx}(2, :)'); hold on;   scatter(Y1_32{idx}(1, :)', Y1_32{idx}(2, :)');  axis equal
        
        Crr_GT_1 = [];   Crr_GT_2 = [];    
        for i=1:size(X1_4, 1)  
            %% 判断是否需要 穷尽搜索 所有潜在匹配
            Seach_range = [i, i];                
            if isfield(NNN, 'Exausted')
                if NNN.Exausted==1
                    Seach_range = [1, size(X2_4, 1)];      
                end
            end
            %% 开始验证
            for j = Seach_range(1) : Seach_range(2)
                
                SR_TX1_32 = SR{j} * TX1_32{i};              SR_Y1_32 = SR_Y1_32_PreCal{j}; 
                SR_TX1_32 = SR_TX1_32([2,1,3],:);           SR_Y1_32 = SR_Y1_32([2,1,3],:);
                Dis = (SR_Y1_32 - SR_TX1_32);               DisAll(:,:,i) = Dis(1:2,1:2); %DisAll = [DisAll, Dis(1:2,1:2)];
                Dis = abs(Dis);
                %% 垂直距离约束
                Crr_GT_1(i,j) = all(Dis(2,:) < threshold);
                %% 角度约束
%                 [theta1, ~] = cart2pol(TX1_32{i}(2, 1)-TX1_32{i}(2, 2),  TX1_32{i}(1, 2)-TX1_32{i}(1, 1));
%                 [theta2, ~] = cart2pol(Y1_32{i}(2, 1)-Y1_32{i}(2, 2),  Y1_32{i}(1, 2)-Y1_32{i}(1, 1));
%                 Crr_GT_1(i,j) = Crr_GT_1(i,j) &&  (abs(wrapToPi(theta2 - theta1)) <pi/9);
                cos_line1 = [TX1_32{i}(2, 2)-TX1_32{i}(2, 1),  TX1_32{i}(1, 2)-TX1_32{i}(1, 1)];
                cos_line2 = [Y1_32{i}(2, 2)-Y1_32{i}(2, 1),  Y1_32{i}(1, 2)-Y1_32{i}(1, 1)];
                cos_2line = (cos_line1*cos_line2')/(norm(cos_line1)*norm(cos_line2));
                Crr_GT_1(i,j) = Crr_GT_1(i,j) &&  (abs(cos_2line) > 0.92);%0.93
                %% 重合度约束
%                 if i==30;
%                     pause(0.1);
%                 end
                Max_length = max(norm(SR_Y1_32(1:2,2)-SR_Y1_32(1:2,1)), norm(SR_TX1_32(1:2,2)-SR_TX1_32(1:2,1)));
                if Crr_GT_1(i,j)
                    aveDrift = mean(abs(SR_Y1_32(1,1:2) - SR_TX1_32(1,1:2)));
                    Lengthratio_temp(i,j) = aveDrift / abs(diff(SR_Y1_32(1,1:2)));
                    Length_temp(i,j) = abs(diff(SR_Y1_32(1,1:2)));
                    cosError(i,j) = abs(cos_2line);
                    if 0 %min(SR_Y1_32(1,:)) - max(SR_TX1_32(1,:)) >Max_length/6 || min(SR_TX1_32(1,:)) - max(SR_Y1_32(1,:)) >Max_length/6
                        %Crr_GT_1(i,j) = all(Dis(2,:) < threshold/4);
%                         Crr_GT_2(i,j) = 0;
%                         Lengthratio_temp(i,j) = 0; % 不重合
%                     elseif (SR_Y1_32(1,2) - SR_Y1_32(1,1))  *  (SR_TX1_32(1,2) - SR_TX1_32(1,1))  <  0
%                         Crr_GT_2(i,j) = 0;
                    else
                        S4_points_x = [SR_TX1_32(1,:), SR_Y1_32(1,:)];
                        %SR_TX1_32:    以直线方向为x轴的坐标，第一行为垂直坐标;
                        Bn1 = SR_Y1_32(1,:) - SR_TX1_32(1,:);
                        Bn2 = SR_Y1_32(1,:) - SR_TX1_32(1,[2,1]);
                        if sum(abs(Bn1)) < sum(abs(Bn2))
                            bn(i,:) = Bn1;
                        else
                            bn(i,:) = Bn2;
                        end
                        Points_4_Sorted = sort(S4_points_x);
                        % 断裂程度：端点距离均值 ÷ 线段长度
                        ratio = norm(Points_4_Sorted([2,3])) / norm(Points_4_Sorted([1,4]));
                        switch length(thr_Line) 
                            case 1      % ground truth
                                Crr_GT_2(i,j) =  ratio<thr_Line  &  ratio>(1/thr_Line);
                            case 2      % picking fracture
                                Crr_GT_2(i,j) =  ratio<=50  &  ratio>=(1/50);
                                if Crr_GT_2(i,j)
                                    lx = norm(SR_TX1_32(1:2,1)-SR_TX1_32(1:2,2));
                                    ly = norm(SR_Y1_32(1:2,1)-SR_Y1_32(1:2,2));
                                    ratio = max(lx,ly)/min(lx,ly);
                                    Crr_GT_2(i,j) =  (ratio>=thr_Line(1)  &&  ratio<=thr_Line(2)) || ((ratio>=1/thr_Line(2)  &&  ratio<=1/thr_Line(1)));
                                end
                         end
                    end
                end
            end
        end
    end
%% 是否考虑 《线段重合度》
    if thr_Line==0
        temp = Crr_GT_1;
        Crr_GT(N(flag-1)+1:N(flag)) = diag(Crr_GT_1);
    else
        temp = Crr_GT_1 & Crr_GT_2;
        Crr_GT(N(flag-1)+1:N(flag)) = diag(Crr_GT_1) & diag(Crr_GT_2);
    end
    
    C = [];
    for i=1:size(X1_4, 1) 
        if any(temp(i,:))
            ind = find(temp(i,:));
            % 如果原对应点为正确对应之一，则选择此对应
            if any(ind == i);       
                ind = i;        
            end
            C = [C;    i, ind(1)];
            Lengthratio(i,:) = [Lengthratio_temp(i, ind(1)), Length_temp(i, ind(1)), cosError(i, ind(1))];
        end
    end
end

%% 9.15
% function [Crr_GT,d1,d2, crsp, Lengthratio, C] = VrfCrrsp(B_opt,  X1, X2, size1, size2, threshold, NNN, thr_Line,Im1, Im2, C_putative, LSM_GT)
% 
% flag=B_opt(end,end);   crsp=[]; Lengthratio=[]; C=[];
% %         F_gt=B_opt;    % threshold=max(0.0025,threshold);
% %         [Crr_GT1,~,~] = DeterminInlier(F_gt,  X2_1, X2_2, size1, size2, threshold);
% %                         [X4_1, X4_2]=LAF_to_pt3x3(X4_1(:,1:4), X4_2(:,1:4)); % 
% %         [Crr_GT2,~,~] = DeterminInlierLine(F_gt,  X4_1, X4_2, size1, size2, threshold);
% %         [Crr_GT2,d1,d2] = DeterminInlier(F_gt,  X6_1, X6_2, size1, size2, threshold);
% 
% if abs(flag+1)<1e-99
%     mean=size1/2;
% 
%     % Make center of the image coordinates 0,0
%     xd=X1(:,1)-mean(1); 
%     yd=X1(:,2)-mean(2);
% 
%     % Calculate the Transformed coordinates
%     Tlocalx = mean(1) + B_opt(1,1) * xd + B_opt(1,2) *yd + B_opt(1,3) * 1;
%     Tlocaly = mean(2) + B_opt(2,1) * xd + B_opt(2,2) *yd + B_opt(2,3) * 1;
%     d2=[]; 
%     threshold=max(6,threshold);
%     d1=sqrt(sum((X2(:,1:2)-[Tlocalx,Tlocaly]).^2,2));
% 
%     Crr_GT=(d1<threshold);
%     
% elseif abs(flag-1)<1e-99
%     
%     T=B_opt([2,1,3],[2,1,3]);     %T=inv(T);
%     d2=[]; 
%     threshold=max(5,threshold);
%     [Crr_GT,d1, crsp]=DetGTusingT(X1,X2,[],[],threshold,T);
% %         threshold = 8*max(size1)/512;
% %         T = B_opt;
% %         [Crr_GT_trial,d1_trial, crsp_trial]=DetGTusingT(X1,X2,[],[],threshold,T);
% %         if sum(Crr_GT_trial)/(sum(Crr_GT)+eps)>2
% %             Crr_GT =Crr_GT_trial;
% %             d1=d1_trial;
% %             crsp=crsp_trial;
% %         end
%     % 线特征Ground Truth
%     if exist('NNN','var')
%         [Crr_GT, Lengthratio, C] = get_LineGroundTruth(B_opt,  X1, X2, size1, size2, threshold,  Crr_GT,  NNN, thr_Line);
%     end
% elseif flag==0
%     % 运行 exp1:  getGroundTruth 估计真值
%     d2=[];          threshold=max(5,threshold);
%     [Crr_GT, d1, crsp]=DetGTusingT(X1,X2,Im1,Im2,threshold);
%     
%     
% else
%     F_gt=B_opt;    % threshold=max(0.0025,threshold);
%     [Crr_GT,d1,d2] = DeterminInlier(F_gt,  X1, X2, size1, size2, threshold);
%     if exist('NNN','var')
%         Crr_GT2 = LSVFDR_verify (X1, X2, NNN);
%         Crr_GT = Crr_GT | Crr_GT2;
%     end
% end
% 
% 
% %% 给定数据集的GroundTryth
% if exist('LSM_GT', 'var')
%     if ~isempty(LSM_GT)
%         Crr_GT_41 = zeros(size(C_putative,1), 1);
%         Crrsp1 = LSM_GT{1};
%         Crrsp2 = LSM_GT{2};
%         % 计算homography
%         for i = 1:size(C_putative,1)
%             Contain_flag1 = any(C_putative(i,1) == Crrsp1,2);
%             if any(Contain_flag1)
%                 Index_in_1 = find(Contain_flag1);
%                 Contain_flag2 = any(C_putative(i,2) == Crrsp2(Index_in_1(1), :),2);
%                 if Contain_flag2
%                     Crr_GT_41(i) = 1;
%                 end
%             end
%         end
%         if size(C_putative,1) == size(Crr_GT,1) %% 如果只有直线特征
%             Crr_GT = Crr_GT_41;
%         else
%             flag = find(NNN.Flag == 41)+1;  %[6, 41, 4, 2] 
%             N = [0, cumsum(NNN.N)];
%             if ~isempty(flag) 
%                 Crr_GT(N(flag-1)+1:N(flag)) = Crr_GT_41;
%             end
%         end
%     end
% end
% 
% 
% Crr_GT = Crr_GT==1;
% 
% end
% 
% function [Crr_GT, Lengthratio, C] = get_LineGroundTruth(B_opt,  X1, X2, size1, size2, threshold,  Crr_GT,  NNN, thr_Line)
% if isfield(NNN, 'Thr_times');   Thr_times = NNN.Thr_times;   else;    Thr_times = 20;    end;
% % 线特征Ground Truth
%     T=B_opt([2,1,3],[2,1,3]);           Lengthratio_all=[];         threshold=min(5,threshold); 
%     Lengthratio=-ones(size(X1,1), 1);   
%     Lengthratio_temp = -ones(size(X1,1), size(X2,1));       Length_temp = zeros(size(X1,1), size(X2,1));
%     flag = find(NNN.Flag == 41)+1;  %[6, 41, 4, 2] 
%     N = [0, cumsum(NNN.N)];
%     if ~isempty(flag)           
%         X1_4 = X1(N(flag-1)+1:N(flag),1:4);               X2_4 = X2(N(flag-1)+1:N(flag),1:4);     
%         
%         if any(X1_4(:,4)>2*pi+1)
%             X2temp = X1_4;
%             Y2temp = X2_4;
%         else
%             [X2temp, Y2temp]=LAF_to_pt3x3(X1_4, X2_4); 
%         end
%         
%         %% 自建数据集（自己估计的真值变换）
% %         if max(size1)>2500
% %             Y2temp = Y2temp(:,[2,1,4,3]);
% %             X2temp = X2temp(:,[2,1,4,3]);
%             threshold = Thr_times*max(size1)/512;
% %         end
%         
%         %% 计算 X 的变换
%         for i = 1 : size(X2_4, 1)
%             X1_32 = [reshape(X2temp(i, :), 2,2); 1, 1];   
%             TX1_32{i} = T * X1_32;     
%             TX1_32{i} = TX1_32{i}./TX1_32{i}(3, :);
%         end
%         %% 计算 Y 的变换
%         for j = 1 : size(X2_4, 1)
%             scaling = [1, 1];
% %             [theta, ~] = cart2pol(Y2temp(j, 4)-Y2temp(j, 2),  Y2temp(j, 3)-Y2temp(j, 1));
%             theta = X2_4(j, 4);
%             rot = [cos(theta), -sin(theta);  -sin(theta), cos(theta)];
%             SR{j} = blkdiag(diag(scaling) * rot, 1);
%             Y1_32{j} = [reshape(Y2temp(j, :), 2,2); 1, 1];  
%             SR_Y1_32_PreCal{j} = SR{j} * Y1_32{j}; 
%         end
%         
%         %% 画图，双点对应变换后的空间关系图
%         % idx = 8;scatter(TX1_32{idx}(1, :)', TX1_32{idx}(2, :)'); hold on;   scatter(Y2temp(idx, [1,3])', Y2temp(idx, [2,4])');  axis equal
%         % idx = 8;scatter(TX1_32{idx}(1, :)', TX1_32{idx}(2, :)'); hold on;   scatter(Y1_32{idx}(1, :)', Y1_32{idx}(2, :)');  axis equal
%         
%         Crr_GT_1 = [];   Crr_GT_2 = [];
%         for i=1:size(X1_4, 1)  
%             %% 判断是否需要 穷尽搜索 所有潜在匹配
%             Seach_range = [i, i];                
%             if isfield(NNN, 'Exausted')
%                 if NNN.Exausted==1
%                     Seach_range = [1, size(X2_4, 1)];      
%                 end
%             end
%             %% 开始验证
%             for j = Seach_range(1) : Seach_range(2)
%                 
%                 SR_TX1_32 = SR{j} * TX1_32{i};              SR_Y1_32 = SR_Y1_32_PreCal{j}; 
%                 SR_TX1_32 = SR_TX1_32([2,1,3],:);           SR_Y1_32 = SR_Y1_32([2,1,3],:);
%                 Dis = abs(SR_Y1_32 - SR_TX1_32);
%                 %% 垂直距离约束
%                 Crr_GT_1(i,j) = all(Dis(2,:) < threshold);
%                 %% 角度约束
% %                 [theta1, ~] = cart2pol(TX1_32{i}(2, 1)-TX1_32{i}(2, 2),  TX1_32{i}(1, 2)-TX1_32{i}(1, 1));
% %                 [theta2, ~] = cart2pol(Y1_32{i}(2, 1)-Y1_32{i}(2, 2),  Y1_32{i}(1, 2)-Y1_32{i}(1, 1));
% %                 Crr_GT_1(i,j) = Crr_GT_1(i,j) &&  (abs(wrapToPi(theta2 - theta1)) <pi/9);
%                 cos_line1 = [TX1_32{i}(2, 1)-TX1_32{i}(2, 2),  TX1_32{i}(1, 2)-TX1_32{i}(1, 1)];
%                 cos_line2 = [Y1_32{i}(2, 1)-Y1_32{i}(2, 2),  Y1_32{i}(1, 2)-Y1_32{i}(1, 1)];
%                 cos_2line = (cos_line1*cos_line2')/(norm(cos_line1)*norm(cos_line2));
%                 Crr_GT_1(i,j) = Crr_GT_1(i,j) &&  (abs(cos_2line) > 0.94);
%                 %% 重合度约束
%                 if i==8;
%                     pause(0.1);
%                 end
%                 Max_length = max(norm(SR_Y1_32(1:2,2)-SR_Y1_32(1:2,1)), norm(SR_TX1_32(1:2,2)-SR_TX1_32(1:2,1)));
%                 if Crr_GT_1(i,j)
%                     aveDrift = mean(abs(SR_Y1_32(1,1:2) - SR_TX1_32(1,1:2)));
%                     Lengthratio_temp(i,j) = aveDrift / abs(diff(SR_Y1_32(1,1:2)));
%                     Length_temp(i,j) = abs(diff(SR_Y1_32(1,1:2)));
%                     if min(SR_Y1_32(1,:)) - max(SR_TX1_32(1,:)) >Max_length/6 || min(SR_TX1_32(1,:)) - max(SR_Y1_32(1,:)) >Max_length/6
%                         Crr_GT_1(i,j) = all(Dis(2,:) < threshold/4);
% %                         Crr_GT_2(i,j) = 0;
% %                         Lengthratio_temp(i,j) = 0; % 不重合
% %                     elseif (SR_Y1_32(1,2) - SR_Y1_32(1,1))  *  (SR_TX1_32(1,2) - SR_TX1_32(1,1))  <  0
% %                         Crr_GT_2(i,j) = 0;
%                     else
%                         S4_points_x = [SR_TX1_32(1,:), SR_Y1_32(1,:)];
%                         Points_4_Sorted = sort(S4_points_x);
%                         % 断裂程度：端点距离均值 ÷ 线段长度
%                         ratio = norm(Points_4_Sorted([2,3])) / norm(Points_4_Sorted([1,4]));
%                         % 断裂程度：变换后线段长度比
% %                         lx = norm(SR_TX1_32(1:2,1)-SR_TX1_32(1:2,2));
% %                         ly = norm(SR_Y1_32(1:2,1)-SR_Y1_32(1:2,2));
% %                         Lengthratio_temp(i,j) = max(lx,ly)/min(lx,ly);  % Fracture ratio
% 
%                         switch length(thr_Line) 
%                             case 1      % ground truth
%                                 Crr_GT_2(i,j) =  ratio<thr_Line  &  ratio>(1/thr_Line);
%                             case 2      % picking fracture
%                                 Crr_GT_2(i,j) =  ratio<=50  &  ratio>=(1/50);
%                                 if Crr_GT_2(i,j)
%                                     lx = norm(SR_TX1_32(1:2,1)-SR_TX1_32(1:2,2));
%                                     ly = norm(SR_Y1_32(1:2,1)-SR_Y1_32(1:2,2));
%                                     ratio = max(lx,ly)/min(lx,ly);
%                                     Lengthratio_temp(i,j) = ratio;
%                                     Crr_GT_2(i,j) =  (ratio>=thr_Line(1)  &&  ratio<=thr_Line(2)) || ((ratio>=1/thr_Line(2)  &&  ratio<=1/thr_Line(1)));
%                                 end
%                          end
%                     end
%                 end
%             end
%         end
%     end
% %% 是否考虑 《线段重合度》
%     if thr_Line==0
%         temp = Crr_GT_1;
%         Crr_GT(N(flag-1)+1:N(flag)) = diag(Crr_GT_1);
%     else
%         temp = Crr_GT_1 & Crr_GT_2;
%         Crr_GT(N(flag-1)+1:N(flag)) = diag(Crr_GT_1) & diag(Crr_GT_2);
%     end
%     
%     C = [];
%     for i=1:size(X1_4, 1) 
%         if any(Lengthratio_temp(i,:)~=0)
% %             Lengthratio(i) = -1;
%         end
%         if any(temp(i,:))
%             ind = find(temp(i,:));
%             % 如果原对应点为正确对应之一，则选择此对应
%             if any(ind == i);       
%                 ind = i;        
%             end
%             C = [C;    i, ind(1)];
%             Lengthratio(i) = Lengthratio_temp(i, ind(1));
% %             Lengthratio_all = [Lengthratio_all, Lengthratio_temp(i, find(temp(i,:)))];
%         end
%     end
%         
% %     Lengthratio = Lengthratio_all;
%     
% 
% end


%% 3,15
% function [Crr_GT,d1,d2, crsp, Lengthratio] = VrfCrrsp(B_opt,  X1, X2, size1, size2, threshold, NNN, thr_Line,Im1, Im2, C_putative, LSM_GT)
% 
% flag=B_opt(end,end);   crsp=[]; Lengthratio=[];
% %         F_gt=B_opt;    % threshold=max(0.0025,threshold);
% %         [Crr_GT1,~,~] = DeterminInlier(F_gt,  X2_1, X2_2, size1, size2, threshold);
% %                         [X4_1, X4_2]=LAF_to_pt3x3(X4_1(:,1:4), X4_2(:,1:4)); % 
% %         [Crr_GT2,~,~] = DeterminInlierLine(F_gt,  X4_1, X4_2, size1, size2, threshold);
% %         [Crr_GT2,d1,d2] = DeterminInlier(F_gt,  X6_1, X6_2, size1, size2, threshold);
% 
% if abs(flag+1)<1e-99
%     mean=size1/2;
% 
%     % Make center of the image coordinates 0,0
%     xd=X1(:,1)-mean(1); 
%     yd=X1(:,2)-mean(2);
% 
%     % Calculate the Transformed coordinates
%     Tlocalx = mean(1) + B_opt(1,1) * xd + B_opt(1,2) *yd + B_opt(1,3) * 1;
%     Tlocaly = mean(2) + B_opt(2,1) * xd + B_opt(2,2) *yd + B_opt(2,3) * 1;
%     d2=[]; 
%     threshold=max(6,threshold);
%     d1=sqrt(sum((X2(:,1:2)-[Tlocalx,Tlocaly]).^2,2));
% 
%     Crr_GT=(d1<threshold);
%     
% elseif abs(flag-1)<1e-99
%     
%     T=B_opt([2,1,3],[2,1,3]); 
%     d2=[]; 
%     threshold=max(5,threshold);
%     [Crr_GT,d1, crsp]=DetGTusingT(X1,X2,[],[],threshold,T);
% %         threshold = 8*max(size1)/512;
%         T = B_opt;
%         [Crr_GT_trial,d1_trial, crsp_trial]=DetGTusingT(X1,X2,[],[],threshold,T);
%         if sum(Crr_GT_trial)/(sum(Crr_GT)+eps)>2
%             Crr_GT =Crr_GT_trial;
%             d1=d1_trial;
%             crsp=crsp_trial;
%         end
%     % 线特征Ground Truth
%     if exist('NNN','var')
%         [Crr_GT, Lengthratio] = get_LineGroundTruth(B_opt,  X1, X2, size1, size2, threshold,  Crr_GT,  NNN, thr_Line);
%     end
% elseif flag==0
%     % 运行 exp1:  getGroundTruth 估计真值
%     d2=[];          threshold=max(5,threshold);
%     [Crr_GT, d1, crsp]=DetGTusingT(X1,X2,Im1,Im2,threshold);
%     
%     
% else
%     F_gt=B_opt;    % threshold=max(0.0025,threshold);
%     [Crr_GT1,d1,d2] = DeterminInlier(F_gt,  X1, X2, size1, size2, threshold);
%     if exist('NNN','var')
%         Crr_GT2 = LSVFDR_verify (X1, X2, NNN);
%     end
%     Crr_GT = Crr_GT1 | Crr_GT2;
% end
% 
% 
% %% 给定数据集的GroundTryth
% if exist('LSM_GT', 'var')
%     if ~isempty(LSM_GT)
%         Crr_GT_41 = zeros(size(C_putative,1), 1);
%         Crrsp1 = LSM_GT{1};
%         Crrsp2 = LSM_GT{2};
%         % 计算homography
%         for i = 1:size(C_putative,1)
%             Contain_flag1 = any(C_putative(i,1) == Crrsp1,2);
%             if any(Contain_flag1)
%                 Index_in_1 = find(Contain_flag1);
%                 Contain_flag2 = any(C_putative(i,2) == Crrsp2(Index_in_1(1), :),2);
%                 if Contain_flag2
%                     Crr_GT_41(i) = 1;
%                 end
%             end
%         end
%         if size(C_putative,1) == size(Crr_GT,1) %% 如果只有直线特征
%             Crr_GT = Crr_GT_41;
%         else
%             flag = find(NNN.Flag == 41)+1;  %[6, 41, 4, 2] 
%             N = [0, cumsum(NNN.N)];
%             if ~isempty(flag) 
%                 Crr_GT(N(flag-1)+1:N(flag)) = Crr_GT_41;
%             end
%         end
%     end
% end
% 
% 
% Crr_GT = Crr_GT==1;
% 
% end
% 
% function [Crr_GT Lengthratio] = get_LineGroundTruth(B_opt,  X1, X2, size1, size2, threshold,  Crr_GT,  NNN, thr_Line)
%     
% % 线特征Ground Truth
%     T=B_opt([2,1,3],[2,1,3]);           threshold=min(5,threshold);             Lengthratio=[];
%     flag = find(NNN.Flag == 41)+1;  %[6, 41, 4, 2] 
%     N = [0, cumsum(NNN.N)];
%     if ~isempty(flag)           
%         X1_4 = X1(N(flag-1)+1:N(flag),1:4);               X2_4 = X2(N(flag-1)+1:N(flag),1:4);     
%         [X2temp, Y2temp]=LAF_to_pt3x3(X1_4, X2_4); 
%         
%         %% 自建数据集（自己估计的真值变换）
% %         if max(size1)>2500
%             Y2temp = Y2temp(:,[2,1,4,3]);
%             X2temp = X2temp(:,[2,1,4,3]);
%             threshold = 8*max(size1)/512;
% %         end
%         
% %         theta = X2(:, 4);
%         for i=1:size(X1_4, 1)  
%             scaling = [1, 1];
%             [theta(i), ~] = cart2pol(Y2temp(i, 3)-Y2temp(i, 1),  Y2temp(i, 4)-Y2temp(i, 2));
%             rot = [cos(theta(i)), sin(theta(i));  -sin(theta(i)), cos(theta(i))];
%             SR = blkdiag(diag(scaling) * rot, 1);
%             X1_32 = [reshape(X2temp(i, :), 2,2); 1, 1];             Y1_32 = [reshape(Y2temp(i, :), 2,2); 1, 1];  
%             SR_TX1_32 = SR * T * X1_32;             SR_Y1_32 = SR * Y1_32; 
%             Dis = abs(SR_Y1_32 - SR_TX1_32);
%             Crr_GT_1(i) = all(Dis(2,:) < threshold);
% %             if prod(Dis(1,:)<0)
% %                 ratio = norm(SR_TX1_32(1,:)) / norm(SR_Y1_32(1,:));
% %                 Crr_GT_2(i) =  ratio<thr_Line  &  ratio>(1/thr_Line);
%             if min(SR_Y1_32(1,:)) - max(SR_TX1_32(1,:)) >0 || min(SR_TX1_32(1,:)) - max(SR_Y1_32(1,:)) >0
%                 Crr_GT_2(i) = 0;
%             elseif (SR_Y1_32(1,2) - SR_Y1_32(1,1))  *  (SR_TX1_32(1,2) - SR_TX1_32(1,1))  <  0
%                 Crr_GT_2(i) = 0;
%             else
%                 S4_points_x = [SR_TX1_32(1,:), SR_Y1_32(1,:)];
%                 Points_4_Sorted = sort(S4_points_x);
%                 ratio = norm(Points_4_Sorted([2,3])) / norm(Points_4_Sorted([1,4]));
%                 
%                 switch length(thr_Line) 
%                     case 1
%                     % ground truth
%                         Crr_GT_2(i) =  ratio<thr_Line  &  ratio>(1/thr_Line);
%                     % picking fracture
%                     case 2
%                         Crr_GT_2(i) =  ratio<5  &  ratio>(1/5);
%                         if Crr_GT_2(i)
%                             lx = norm(SR_TX1_32(1:2,1)-SR_TX1_32(1:2,2));
%                             ly = norm(SR_Y1_32(1:2,1)-SR_Y1_32(1:2,2));
%                             ratio = max(lx,ly)/min(lx,ly);
%                             Lengthratio(i) = ratio;
%                             Crr_GT_2(i) =  (ratio>=thr_Line(1)  &&  ratio<=thr_Line(2)) || ((ratio>=1/thr_Line(2)  &&  ratio<=1/thr_Line(1)));
%                         end
%                 end
%             end
%         end
%     end
%     Crr_GT(N(flag-1)+1:N(flag)) = Crr_GT_1 & Crr_GT_2;
% end




%% 
% % function [Crr_GT,d1,d2, crsp] = VrfCrrsp(B_opt,  X1, X2, size1, size2, threshold, NNN, thr_Line,Im1, Im2, C_putative, LSM_GT)
% % 
% % flag=B_opt(end,end);   crsp=[];
% % %         F_gt=B_opt;    % threshold=max(0.0025,threshold);
% % %         [Crr_GT1,~,~] = DeterminInlier(F_gt,  X2_1, X2_2, size1, size2, threshold);
% % %                         [X4_1, X4_2]=LAF_to_pt3x3(X4_1(:,1:4), X4_2(:,1:4)); % 
% % %         [Crr_GT2,~,~] = DeterminInlierLine(F_gt,  X4_1, X4_2, size1, size2, threshold);
% % %         [Crr_GT2,d1,d2] = DeterminInlier(F_gt,  X6_1, X6_2, size1, size2, threshold);
% % 
% % if abs(flag+1)<1e-99
% %     mean=size1/2;
% % 
% %     % Make center of the image coordinates 0,0
% %     xd=X1(:,1)-mean(1); 
% %     yd=X1(:,2)-mean(2);
% % 
% %     % Calculate the Transformed coordinates
% %     Tlocalx = mean(1) + B_opt(1,1) * xd + B_opt(1,2) *yd + B_opt(1,3) * 1;
% %     Tlocaly = mean(2) + B_opt(2,1) * xd + B_opt(2,2) *yd + B_opt(2,3) * 1;
% %     d2=[]; 
% %     threshold=max(6,threshold);
% %     d1=sqrt(sum((X2(:,1:2)-[Tlocalx,Tlocaly]).^2,2));
% % 
% %     Crr_GT=(d1<threshold);
% %     
% % elseif abs(flag-1)<1e-99
% %     
% %     T=B_opt([2,1,3],[2,1,3]); 
% %     d2=[]; threshold=max(5,threshold);
% %     [Crr_GT,d1, crsp]=DetGTusingT(X1,X2,[],[],threshold,T);
% %     % 线特征Ground Truth
% %     if exist('NNN','var')
% %         Crr_GT = get_LineGroundTruth(B_opt,  X1, X2, size1, size2, threshold,  Crr_GT,  NNN, thr_Line);
% %     end
% % elseif flag==0
% %     % 运行 exp1:  getGroundTruth 估计真值
% %     d2=[];          threshold=max(5,threshold);
% %     [Crr_GT, d1, crsp]=DetGTusingT(X1,X2,Im1,Im2,threshold);
% %     
% %     
% % else
% %     F_gt=B_opt;    % threshold=max(0.0025,threshold);
% %     [Crr_GT1,d1,d2] = DeterminInlier(F_gt,  X1, X2, size1, size2, threshold);
% %     if exist('NNN','var')
% %         Crr_GT2 = LSVFDR_verify (X1, X2, NNN);
% %     end
% %     Crr_GT = Crr_GT1 | Crr_GT2;
% % end
% % 
% % 
% % %% 给定数据集的GroundTryth
% % if exist('LSM_GT', 'var')
% %     if ~isempty(LSM_GT)
% %         Crr_GT_41 = zeros(size(C_putative,1), 1);
% %         Crrsp1 = LSM_GT{1};
% %         Crrsp2 = LSM_GT{2};
% %         % 计算homography
% %         for i = 1:size(C_putative,1)
% %             Contain_flag1 = any(C_putative(i,1) == Crrsp1,2);
% %             if any(Contain_flag1)
% %                 Index_in_1 = find(Contain_flag1);
% %                 Contain_flag2 = any(C_putative(i,2) == Crrsp2(Index_in_1(1), :),2);
% %                 if Contain_flag2
% %                     Crr_GT_41(i) = 1;
% %                 end
% %             end
% %         end
% %         if size(C_putative,1) == size(Crr_GT,1) %% 如果只有直线特征
% %             Crr_GT = Crr_GT_41;
% %         else
% %             flag = find(NNN.Flag == 41)+1;  %[6, 41, 4, 2] 
% %             N = [0, cumsum(NNN.N)];
% %             if ~isempty(flag) 
% %                 Crr_GT(N(flag-1)+1:N(flag)) = Crr_GT_41;
% %             end
% %         end
% %     end
% % end
% % 
% % end
% % 
% % function Crr_GT = get_LineGroundTruth(B_opt,  X1, X2, size1, size2, threshold,  Crr_GT,  NNN, thr_Line)
% %     
% % % 线特征Ground Truth
% %     T=B_opt([2,1,3],[2,1,3]);           threshold=min(5,threshold);
% %     flag = find(NNN.Flag == 41)+1;  %[6, 41, 4, 2] 
% %     N = [0, cumsum(NNN.N)];
% %     if ~isempty(flag)           
% %         X1_4 = X1(N(flag-1)+1:N(flag),1:4);               X2_4 = X2(N(flag-1)+1:N(flag),1:4);     
% %         [X2temp, Y2temp]=LAF_to_pt3x3(X1_4, X2_4); 
% %         
% %         %% 自建数据集（自己估计的真值变换）
% % %         if max(size1)>2500
% %             Y2temp = Y2temp(:,[2,1,4,3]);
% %             X2temp = X2temp(:,[2,1,4,3]);
% %             threshold = 8*max(size1)/512;
% % %         end
% %         
% % %         theta = X2(:, 4);
% %         for i=1:size(X1_4, 1)  
% %             scaling = [1, 1];
% %             [theta(i), ~] = cart2pol(Y2temp(i, 3)-Y2temp(i, 1),  Y2temp(i, 4)-Y2temp(i, 2));
% %             rot = [cos(theta(i)), sin(theta(i));  -sin(theta(i)), cos(theta(i))];
% %             SR = blkdiag(diag(scaling) * rot, 1);
% %             X1_32 = [reshape(X2temp(i, :), 2,2); 1, 1];             Y1_32 = [reshape(Y2temp(i, :), 2,2); 1, 1];  
% %             SR_TX1_32 = SR * T * X1_32;             SR_Y1_32 = SR * Y1_32; 
% %             Dis = abs(SR_Y1_32 - SR_TX1_32);
% %             Crr_GT_1(i) = all(Dis(2,:) < threshold);
% % %             if prod(Dis(1,:)<0)
% % %                 ratio = norm(SR_TX1_32(1,:)) / norm(SR_Y1_32(1,:));
% % %                 Crr_GT_2(i) =  ratio<thr_Line  &  ratio>(1/thr_Line);
% %             if min(SR_Y1_32(1,:)) - max(SR_TX1_32(1,:)) >0 || min(SR_TX1_32(1,:)) - max(SR_Y1_32(1,:)) >0
% %                 Crr_GT_2(i) = 0;
% %             elseif (SR_Y1_32(1,2) - SR_Y1_32(1,1))  *  (SR_TX1_32(1,2) - SR_TX1_32(1,1))  <  0
% %                 Crr_GT_2(i) = 0;
% %             else
% %                 S4_points_x = [SR_TX1_32(1,:), SR_Y1_32(1,:)];
% %                 Points_4_Sorted = sort(S4_points_x);
% %                 ratio = norm(Points_4_Sorted([2,3])) / norm(Points_4_Sorted([1,4]));
% %                 
% %                 switch length(thr_Line) 
% %                     case 1
% %                     % ground truth
% %                         Crr_GT_2(i) =  ratio<thr_Line  &  ratio>(1/thr_Line);
% %                     % picking fracture
% %                     case 2
% %                         Crr_GT_2(i) =  (ratio>thr_Line(1)  &&  ratio<thr_Line(2)) || ((ratio>1/thr_Line(2)  &&  ratio<1/thr_Line(1)));
% %                 end
% %             end
% %         end
% %     end
% %     Crr_GT(N(flag-1)+1:N(flag)) = Crr_GT_1 & Crr_GT_2;
% % 
% % end