function [B_error1, ratio] = cal_B_err(B_est, B_GT)
    B_error1 = B_est - B_GT;
    B_error1 = B_error1;
    B_error2 = B_est + B_GT;
    B_error2 = B_error2;
    idx_2 = (sum(abs(B_error1), 2)>sum(abs(B_error2), 2));
    B_error1(idx_2,:) = B_error2(idx_2,:);
    ratio = sqrt(sum(B_error1.^2, 2)./sum(B_GT.^2, 2));
            