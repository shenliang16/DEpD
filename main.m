
% Pts1 = X4p_1(:,1:4); 
% Pts2 = X4p_2(:,1:4); 
% Line2Pts1 = X4_1;
% Line2Pts2 = X4_2;
% NNDR_kps = Sift_Ratio_4p;
% save TestData I1 I2 NNDR_kps LBDdistance Pts1 Pts2 Line2Pts1 Line2Pts2 Crr_GT_4 bn_LBD thr_ratio_Line NNN_GT B_opt thre


load('TestData.mat')
        matcher = 'DEpD-pts-C-0.5';     % The method
        DscDis.point = NNDR_kps;        % Key points 1st-to-2ed descriptor distance ratior
        DscDis.line = LBDdistance;      % LBD descriptor distance
        Bilateral = 1;                  % Bilateral Drifting or lateral Drifting
        anneling = 0;                   
        P_thr = 0.75;                   % decision threshold
        Thr_kps = 0.5;
        size1 = size(I1);
        size2 = size(I2);

        if strcmp(matcher, 'DEpD0')||strcmp(matcher, ['DEpD0-', num2str(Thr_kps)])
            Pts1 = [];                  Pts2 = [];                  gamma=0;        Constraints = '0';
        elseif strcmp(matcher, 'DEpD-L2')||strcmp(matcher, ['DEpD-L2-', num2str(Thr_kps)])
            Pts1 = [];                  Pts2 = [];                  gamma=1;        Constraints = 'L2';
        elseif strcmp(matcher, 'DEpD-LenL2')||strcmp(matcher, ['DEpD-LenL2-', num2str(Thr_kps)])
            Pts1 = [];                  Pts2 = [];                  gamma=1;        Constraints = 'L2-Length';
        elseif strcmp(matcher, 'DEpD-Ratio')||strcmp(matcher, ['DEpD-Ratio-', num2str(Thr_kps)])
            Pts1 = [];                  Pts2 = [];                  gamma=1;        Constraints = 'L2-Ratio';
        elseif strcmp(matcher, 'DEpD-pts')||strcmp(matcher, ['DEpD-pts-', num2str(Thr_kps)])
            gamma=1;                    Constraints = 'L2-pts';
        else
            gamma=[1, 1, 1];            Constraints = 'L2-ptsScale';
        end
        %% The DEpD Method
        [C, time, H, ~, B, omiga2, iter] = DEpD_Homo(Pts1, Pts2,  Line2Pts1,Line2Pts2,  anneling, P_thr, gamma, Constraints, Bilateral, DscDis);
        %% Evaluation
            bn_GT = bn_LBD;
            [~, ~, ~, ~, ~,~,~,bn_est] = VrfCrrsp(H, Line2Pts1,Line2Pts2, size1(1:2), size2(1:2), thre, NNN_GT, thr_ratio_Line);
            [B_err, ratio] = cal_B_err(bn_est(Crr_GT_4,:), bn_LBD(Crr_GT_4,:));
        crspline = C.line;
        if length(crspline)>2
            NNN_GT.N=[size(crspline,1)];
            NNN_GT.Flag = [41];
            [Crr_line, ~, ~, ~, Lengthratio] = VrfCrrsp(B_opt, Line2Pts1(crspline(:,1),:), Line2Pts2(crspline(:,2),:), size1(1:2), size2(1:2), thre, NNN_GT, thr_ratio_Line);
            crspline(Lengthratio==-1,:)=[];  Crr_line(Lengthratio==-1,:)=[]; 
            Matching_Plot_PLR(I1,I2,crspline,Crr_line,[],Line2Pts1,Line2Pts2);    set(gca, 'LooseInset', [0,0,0,0])
            P=sum(Crr_line)/length(Crr_line);
            R=sum(Crr_line)/sum(Crr_GT_4);
            F= 2*R.*P./(P+R); 
            XX1=Line2Pts1(crspline(:,1),:); 
            XX2=Line2Pts2(crspline(:,2),:); 
        else; XX1=[];  XX2=[];  P=0;  R=0;  F= 0;
        end