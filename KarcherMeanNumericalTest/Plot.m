function Plot

% specify the maximum time of the plots
max_time = 1e-2; 

% specify the maximum number of iterations in the plots
max_iter = 50;

% specify the lower and upper limit of the plots
y_lower = 1e-16; 
y_upper = 100;

% specify the data path
DataPath = sprintf('L2Results/k3n3well/');

TestMethods = {'RSDBB', 'RSDQR', 'LRBFGS', 'BINI', 'RBFGS'};
% TestMethods = {'RSDBB', 'LRBFGS'};

MemorySize = {0, 2, 4};

%% =========================== Plot iter vs distance ===========================
figure(1)
h = semilogy(1,1);
set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
hold on;
set(gca,'FontSize',18)
legendInfo = {};


% =========================== Plot Bini ===========================
if(ismember('BINI', TestMethods))
    BINI = importdata([DataPath 'AveBINI.mat']);
    plot(BINI.dis_iter, '^-', 'Markersize', 12, 'linewidth',2);
    legendInfo = union(legendInfo, 'RL Iteration', 'stable');
end


% =========================== Plot RSDQR ===========================
if(ismember('RSDQR', TestMethods))
    RSDQR = importdata([DataPath 'AveRSDQR.mat']);
    plot(RSDQR.dis_iter, '*-', 'Markersize', 12, 'linewidth',2);
    legendInfo = union(legendInfo, 'RSD-QR', 'stable');
end


% =========================== Plot RSDBB ===========================
if(ismember('RSDBB', TestMethods))
    RSDBB = importdata([DataPath 'AveLRBFGSm0.mat']);
    plot(RSDBB.dis_iter, 's-', 'Markersize', 12,'linewidth',2);
    legendInfo = union(legendInfo, 'RBB', 'stable');
end


% =========================== Plot LRBFGS ===========================
if(ismember('LRBFGS', TestMethods))
    for m = 2 : length(MemorySize)
        LRBFGS = importdata([DataPath 'AveLRBFGSm' int2str(MemorySize{m}) '.mat']);
        plot(LRBFGS.dis_iter, 'o-', 'Markersize', 12,'linewidth',2);
        legendInfo = union(legendInfo, ['LRBFGS: ' int2str(MemorySize{m})], 'stable');
    end
end


% =========================== Plot RBFGS ===========================
if(ismember('RBFGS', TestMethods))
    RBFGS = importdata([DataPath 'AveRBFGS.mat']);
    plot(RBFGS.dis_iter, 'v-', 'Markersize', 12, 'linewidth',2);
    legendInfo = union(legendInfo, 'RBFGS', 'stable');
end


l = legend(legendInfo);
set(l,'Interpreter','Latex');
xlabel('iterations','FontSize',20,'Interpreter','Latex');
ylabel('$\mathrm{dist}(\mu, X_t)$','FontSize',20,'Interpreter','Latex');
xlim([0 max_iter]);
ylim([y_lower y_upper]);




 
%% =========================== Plot time vs distance ===========================
figure(2)
h = semilogy(1,100);
set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
hold on;
set(gca,'FontSize',18);
legendInfo = {};


% =========================== Plot Bini ===========================
if(ismember('BINI', TestMethods))
    BINI = importdata([DataPath 'AveBINI.mat']);
    plot(BINI.dis_t, BINI.dis, 'linewidth', 3);
    legendInfo = union(legendInfo, 'RL Iteration', 'stable');
end


% =========================== Plot RSDQR ===========================
if(ismember('RSDQR', TestMethods))
    RSDQR = importdata([DataPath 'AveRSDQR.mat']);
    plot(RSDQR.dis_t, RSDQR.dis, 'linewidth', 3);
    legendInfo = union(legendInfo, 'RSD-QR', 'stable');
end


% =========================== Plot RSDBB ===========================
if(ismember('RSDBB', TestMethods))
    RSDBB = importdata([DataPath 'AveLRBFGSm0.mat']);
    plot(RSDBB.dis_t, RSDBB.dis, 'linewidth', 3);
    legendInfo = union(legendInfo, 'RBB', 'stable');
end

% =========================== Plot LRBFGS ===========================
if(ismember('LRBFGS', TestMethods))
    for m = 2 : length(MemorySize)
        LRBFGS = importdata([DataPath 'AveLRBFGSm' int2str(MemorySize{m}) '.mat']);
        stairs(LRBFGS.dis_t, LRBFGS.dis,'linewidth',3);
        legendInfo = union(legendInfo, ['LRBFGS: ' int2str(MemorySize{m})], 'stable');
    end
end


% =========================== Plot RBFGS ===========================
if(ismember('RBFGS', TestMethods))
    RBFGS = importdata([DataPath 'AveRBFGS.mat']);
    stairs(RBFGS.dis_t, RBFGS.dis, 'linewidth',3);
    legendInfo = union(legendInfo, 'RBFGS', 'stable');
end


l = legend(legendInfo);
set(l,'Interpreter','Latex');
xlabel('time (s)','FontSize',20,'Interpreter','Latex');
ylabel('$\mathrm{dist}(\mu, X_t)$','FontSize',20,'Interpreter','Latex');
xlim([0 max_time]);
ylim([y_lower y_upper]);





end