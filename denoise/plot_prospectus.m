
n = 128;
M = load_image('duorou3', 256);
M = M(:, :, 1);
M = rescale(crop(M,n));
T = importdata('duorou3_T.mat');
options.sub = 5;

r = importdata('denoise_results_pr01.mat');
T_hat_10 = importdata('duorou3_T_hat_pr01.mat');
T_hat = T_hat_10{1};

% plot the noisy image
figure(1)
plot_tensor_field(T_hat, M, options);
title('Noisy tensor image','FontSize',25,'Interpreter','Latex');


dT = r.dT_euclid{1};
U_dnoise = perform_tensor_mapping(dT, +1);
U_dnoise(:,:,1) = perform_histogram_equalization(U_dnoise(:,:,1), 'linear');
T1_dnoise = perform_tensor_mapping(U_dnoise,-1);
figure(2)
plot_tensor_field(T1_dnoise, M, options);
title('Denoised tensor image','FontSize',25,'Interpreter','Latex');


dT = r.dT_logeuclid{1};
U_dnoise = perform_tensor_mapping(dT, +1);
U_dnoise(:,:,1) = perform_histogram_equalization(U_dnoise(:,:,1), 'linear');
T1_dnoise = perform_tensor_mapping(U_dnoise,-1);
figure(3)
plot_tensor_field(T1_dnoise, M, options);
title('Denoised tensor image','FontSize',25,'Interpreter','Latex');



dT = r.dT_symmldbreg{1};
U_dnoise = perform_tensor_mapping(dT, +1);
U_dnoise(:,:,1) = perform_histogram_equalization(U_dnoise(:,:,1), 'linear');
T1_dnoise = perform_tensor_mapping(U_dnoise,-1);
figure(4)
plot_tensor_field(T1_dnoise, M, options);
title('Denoised tensor image','FontSize',25,'Interpreter','Latex');



dT = r.dT_riemann{1};
for i = 1 : n
    for j = 1 : n
        flag = sum(sum(isnan(dT(i, j, :,:))));
        if(flag > 0)
            dT(i, j, :, :) = T_hat(i, j, :, :);
        end
    end
end
U_dnoise = perform_tensor_mapping(dT, +1);
U_dnoise(:,:,1) = perform_histogram_equalization(U_dnoise(:,:,1), 'linear');
T1_dnoise = perform_tensor_mapping(U_dnoise,-1);
figure(5)
plot_tensor_field(T1_dnoise, M, options);
title('Denoised tensor image','FontSize',25,'Interpreter','Latex');




dT = r.dT_medianriemann{1};
for i = 1 : n
    for j = 1 : n
        flag = sum(sum(isnan(dT(i, j, :,:))));
        if(flag > 0)
            dT(i, j, :, :) = T_hat(i, j, :, :);
        end
    end
end
U_dnoise = perform_tensor_mapping(dT, +1);
U_dnoise(:,:,1) = perform_histogram_equalization(U_dnoise(:,:,1), 'linear');
T1_dnoise = perform_tensor_mapping(U_dnoise,-1);
figure(6)
plot_tensor_field(T1_dnoise, M, options);
title('Denoised tensor image','FontSize',25,'Interpreter','Latex');




dT = r.dT_alpha0{1};
for i = 1 : n
    for j = 1 : n
        flag = sum(sum(isnan(dT(i, j, :,:))));
        if(flag > 0)
            dT(i, j, :, :) = T_hat(i, j, :, :);
        end
    end
end
U_dnoise = perform_tensor_mapping(dT, +1);
U_dnoise(:,:,1) = perform_histogram_equalization(U_dnoise(:,:,1), 'linear');
T1_dnoise = perform_tensor_mapping(U_dnoise,-1);
figure(7)
plot_tensor_field(T1_dnoise, M, options);
title('Denoised tensor image','FontSize',25,'Interpreter','Latex');




dT = r.dT_medianalpha1{1};
for i = 1 : n
    for j = 1 : n
        flag = sum(sum(isnan(dT(i, j, :,:))));
        if(flag > 0)
            dT(i, j, :, :) = T_hat(i, j, :, :);
        end
    end
end
U_dnoise = perform_tensor_mapping(dT, +1);
U_dnoise(:,:,1) = perform_histogram_equalization(U_dnoise(:,:,1), 'linear');
T1_dnoise = perform_tensor_mapping(U_dnoise,-1);
figure(8)
plot_tensor_field(T1_dnoise, M, options);
title('Denoised tensor image','FontSize',25,'Interpreter','Latex');