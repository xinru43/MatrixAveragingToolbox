function TestStepsizeVSEigV2

% investigate the relationship between the BB stepsize and eigenvalues of
% the Hessian

%% import data including initial point

% DataPath = sprintf('DataStein/try/');
% DataPath = sprintf('DataStein/k30n30ill/');
% Data = importdata([DataPath 'Kmean1.mat']);
% As = Data.As;
% n = size(As{1}, 1);
% k = length(As);


n = 30;
k = 10;
Type = 1;
As = DataGen(n, k, Type, 1);





Ls = zeros(n, n, k);
for i = 1 : k
    Ls(:, :, i) = As{i};
end


%% initial iterate
rng(37);
x1 = randn(n);
x1 = x1 * x1';

% x1 = diag([1 1e-3 1e-4]); % keep it, save plot!
%x1 = diag([1e-3 1 1e-5]);
%  [O, ~] = qr(randn(n));
%  x1 = O * x1 * O';


% a = ones(1, 5);
% x1 = diag([a, a * 1e-4, a, a*1e-4, a, a*5e-3]);
% [O, ~] = qr(n);
% x1 = O * x1 * O';
% 
% 
% indxeven = 2:2:n;
% indxodd = 1:2:n-1;
% elemeven = (1:15)*1e+2;
% elemodd = (1:15)*1e-5;
% D(indxeven) = elemeven;
% D(indxodd) = elemodd;
% [O, ~] = qr(n);
% x1 = O * diag(D) * O';

% x1 = computeAH(As, 0); % Type = 0: arithmetic mean; Type = 1: arithmetic-harmonic mean



%% compute Xtrue
Xtrue = ComputeXtrueSteinV2(As, 1e-16);



%% set parameters
SolverParams.IsCheckParams = 0;
SolverParams.Tolerance = 1e-17;
SolverParams.Accuracy = 1e-5;
%     SolverParams.Num_pre_funs = 5;
SolverParams.LS_ratio1 = 0.5;
SolverParams.LS_ratio2 = 0.5;
SolverParams.Minstepsize = eps;
SolverParams.Maxstepsize = 100;
SolverParams.Initstepsize = 1;
SolverParams.Finalstepsize = -1;
SolverParams.LineSearch_LS = 0; % 'ARMIJO'
SolverParams.DEBUG = 3;
SolverParams.HasHHR = 1;
SolverParams.isconvex = 1;
SolverParams.Max_Iteration = 60;
RBBiter = 60;
ABBiter = 50;
LRBFGSiter = 60;
FPiter = 60;



TestMethods = {'BB1', 'BB2', 'ABB'};
BB1 = {'m0WL', 'm0', 'm4', 'm16'};
BB2 = {'m0WL', 'm0', 'm4', 'm16'};
ABB = {'m0WL', 'm0', 'm4', 'm16'};


%% compute eigenvalues of Hessian
ischeck = 0;
if (ischeck == 1)
    SolverParams.DEBUG = 1;
    for i = 0 : 60
        SolverParams.method = 'RSD';
        SolverParams.Max_Iteration = i;
        SolverParams.BBratio = 0.8;
        SolverParams.Num_pre_BB = 5;
        SolverParams.InitSteptype = 1; % one
        Xopt = TestSPDMeanSteinV2(Ls, x1, 1, SolverParams, Xtrue);
        X = reshape(Xopt.main, n, n);
        [cond_Hess(i + 1), min_eig(i+1), max_eig(i+1)] = ComputeSteinCond(As, X);
%         [~, ~, ~, ~, Eig, ~] = ComputeSPDMeanSteinV2(As, x1, SolverParams, Xtrue);
%         min_eigWen(i+1) = Eig(3);
%         max_eigWen(i+1) = Eig(4);
    end
    
    % interpolation
    x = 1 : length(max_eig);
    ix = 1 : 0.2 : length(max_eig);
    imax_eig = interp1(x, max_eig, ix);
    imin_eig = interp1(x, min_eig, ix);
%     imax_eigWen = interp1(x, max_eigWen, ix,'spline');
%     imin_eigWen = interp1(x, min_eigWen, ix,'spline');
    
    save('./EigInfo', 'cond_Hess', 'max_eig', 'min_eig', 'imax_eig', 'imin_eig', 'x', 'ix');
end




% %% fixed point iteration
% if (ismember('FP', TestMethods))
%     AsFP = [];
%     for i = 1 : k
%         AsFP = [AsFP reshape(As{i}, 1, n * n)];
%     end
%     x1FP = reshape(x1, 1, n * n);
%     XtrueFP = reshape(Xtrue, 1, n * n);
%     Max_Iteration = FPiter;
%     Debug = 3;
%     [T0, ComTime0, iter, X, disFP, ~] = TestSPDSteinFP(AsFP, x1FP, XtrueFP, n, k, Max_Iteration, 1e-20, Debug);
% end


save('data');
 
%% BB1 family
if (ismember('BB1', TestMethods))
    
    % BB1 parameter
    SolverParams.BBratio = 0;
    SolverParams.Num_pre_BB = 0;
    
    % BB1 without locking condition
    if(ismember('m0WL', BB1))
        SolverParams.method = 'RSD';
        SolverParams.InitSteptype = 1; % BB stepsize
        [disBB1WL, ihBB1WL, hBB1WL, gradsBB1WL, ~] = ComputeSPDMeanSteinV2(As, x1, SolverParams, Xtrue);
        save('data', 'disBB1WL', 'ihBB1WL', 'hBB1WL', 'gradsBB1WL','-append');
    end
    
    
    if(ismember('m0', BB1))
        SolverParams.method = 'LRBFGS';
        SolverParams.LengthSY = 0;
        SolverParams.InitSteptype = 0; % one
        [disBB1m0, ihBB1m0, hBB1m0, gradsBB1m0, ~] = ComputeSPDMeanSteinV2(As, x1, SolverParams, Xtrue);
        save('data', 'disBB1m0', 'ihBB1m0', 'hBB1m0', 'gradsBB1m0','-append');
    end
    
    
    if(ismember('m4', BB1))
        SolverParams.method = 'LRBFGS';
        SolverParams.LengthSY = 4;
        SolverParams.InitSteptype = 0; % one
        [disBB1m4, ihBB1m4, hBB1m4, gradsBB1m4, ~] = ComputeSPDMeanSteinV2(As, x1, SolverParams, Xtrue);
        save('data', 'disBB1m4', 'ihBB1m4', 'hBB1m4', 'gradsBB1m4','-append');
    end
    
    
    if(ismember('m16', BB1))
        SolverParams.method = 'LRBFGS';
        SolverParams.LengthSY = 16;
        SolverParams.InitSteptype = 0; % one
        [disBB1m16, ihBB1m16, hBB1m16, gradsBB1m16, ~] = ComputeSPDMeanSteinV2(As, x1, SolverParams, Xtrue);
        save('data', 'disBB1m16', 'ihBB1m16', 'hBB1m16', 'gradsBB1m16','-append');
    end
    
end



%% BB2 family
if (ismember('BB2', TestMethods))
    
    % BB2 parameter
    SolverParams.BBratio = 1;
    SolverParams.Num_pre_BB = 0;
    
    % BB2 without locking condition
    if(ismember('m0WL', BB2))
        SolverParams.method = 'RSD';
        SolverParams.InitSteptype = 1; % BB stepsize
        [disBB2WL, ihBB2WL, hBB2WL, gradsBB2WL, ~] = ComputeSPDMeanSteinV2(As, x1, SolverParams, Xtrue);
        save('data', 'disBB2WL', 'ihBB2WL', 'hBB2WL', 'gradsBB2WL','-append');
    end
    
    
    if(ismember('m0', BB2))
        SolverParams.method = 'LRBFGS';
        SolverParams.LengthSY = 0;
        SolverParams.InitSteptype = 0; % one
        [disBB2m0, ihBB2m0, hBB2m0, gradsBB2m0, ~] = ComputeSPDMeanSteinV2(As, x1, SolverParams, Xtrue);
        save('data', 'disBB2m0', 'ihBB2m0', 'hBB2m0', 'gradsBB2m0','-append');
    end
    
    
    if(ismember('m4', BB2))
        SolverParams.method = 'LRBFGS';
        SolverParams.LengthSY = 4;
        SolverParams.InitSteptype = 0; % one
        [disBB2m4, ihBB2m4, hBB2m4, gradsBB2m4, ~] = ComputeSPDMeanSteinV2(As, x1, SolverParams, Xtrue);
        save('data', 'disBB2m4', 'ihBB2m4', 'hBB2m4', 'gradsBB2m4','-append');
    end
    
    
    if(ismember('m16', BB2))
        SolverParams.method = 'LRBFGS';
        SolverParams.LengthSY = 16;
        SolverParams.InitSteptype = 0; % one
        [disBB2m16, ihBB2m16, hBB2m16, gradsBB2m16, ~] = ComputeSPDMeanSteinV2(As, x1, SolverParams, Xtrue);
        save('data', 'disBB2m16', 'ihBB2m16', 'hBB2m16', 'gradsBB2m16','-append');
    end
    
end






%% ABB family
if (ismember('ABB', TestMethods))
    
    % BB2 parameter
    SolverParams.BBratio = 0.5;
    SolverParams.Num_pre_BB = 10;
    
    % BB2 without locking condition
    if(ismember('m0WL', ABB))
        SolverParams.method = 'RSD';
        SolverParams.InitSteptype = 1; % BB stepsize
        [disABBWL, ihABBWL, hABBWL, gradsABBWL, ~] = ComputeSPDMeanSteinV2(As, x1, SolverParams, Xtrue);
        save('data', 'disABBWL', 'ihABBWL', 'hABBWL', 'gradsABBWL','-append');
    end
    
    
    if(ismember('m0', ABB))
        SolverParams.method = 'LRBFGS';
        SolverParams.LengthSY = 0;
        SolverParams.InitSteptype = 0; % one
        [disABBm0, ihABBm0, hABBm0, gradsABBm0, ~] = ComputeSPDMeanSteinV2(As, x1, SolverParams, Xtrue);
        save('data', 'disABBm0', 'ihABBm0', 'hABBm0', 'gradsABBm0','-append');
    end
    
    
    if(ismember('m4', ABB))
        SolverParams.method = 'LRBFGS';
        SolverParams.LengthSY = 4;
        SolverParams.InitSteptype = 0; % one
        [disABBm4, ihABBm4, hABBm4, gradsABBm4, ~] = ComputeSPDMeanSteinV2(As, x1, SolverParams, Xtrue);
        save('data', 'disABBm4', 'ihABBm4', 'hABBm4', 'gradsABBm4','-append');
    end
    
    
    if(ismember('m16', ABB))
        SolverParams.method = 'LRBFGS';
        SolverParams.LengthSY = 16;
        SolverParams.InitSteptype = 0; % one
        [disABBm16, ihABBm16, hABBm16, gradsABBm16, ~] = ComputeSPDMeanSteinV2(As, x1, SolverParams, Xtrue);
        save('data', 'disABBm16', 'ihABBm16', 'hABBm16', 'gradsABBm16','-append');
    end
    
end



figure(1)
h = semilogy(1,1);
hold on;
set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
legendInfo = {};
if(ismember('BB1', TestMethods))
    if(ismember('m0WL', BB1))
        plot(disBB1WL,'-^', 'markersize',12,'linewidth',2);
        legendInfo = union(legendInfo, 'RBB1-w/o locking', 'stable'); 
    end
    if(ismember('m0', BB1))
        plot(disBB1m0,'-^', 'markersize',12,'linewidth',2);
        legendInfo = union(legendInfo, 'LRBFGS-BB1: 0', 'stable'); 
    end
    if(ismember('m4', BB1))
        plot(disBB1m4,'-^', 'markersize',12,'linewidth',2);
        legendInfo = union(legendInfo, 'LRBFGS-BB1: 4', 'stable');
    end
    if(ismember('m16', BB1))
        plot(disBB1m16,'-^', 'markersize',12,'linewidth',2);
        legendInfo = union(legendInfo, 'LRBFGS-BB1: 16', 'stable');
    end
end


if(ismember('BB2', TestMethods))
    if(ismember('m0WL', BB2))
        plot(disBB2WL,'-o', 'markersize',12,'linewidth',2);
        legendInfo = union(legendInfo, 'RBB2-w/o locking', 'stable'); 
    end
    if(ismember('m0', BB2))
        plot(disBB2m0,'-o', 'markersize',12,'linewidth',2);
        legendInfo = union(legendInfo, 'LRBFGS-BB2: 0', 'stable'); 
    end
    if(ismember('m4', BB2))
        plot(disBB2m4,'-o', 'markersize',12,'linewidth',2);
        legendInfo = union(legendInfo, 'LRBFGS-BB2: 4', 'stable');
    end
    if(ismember('m16', BB2))
        plot(disBB2m16,'-o', 'markersize',12,'linewidth',2);
        legendInfo = union(legendInfo, 'LRBFGS-BB2: 16', 'stable');
    end
end


if(ismember('ABB', TestMethods))
    if(ismember('m0WL', ABB))
        plot(disABBWL,'-*', 'markersize',12,'linewidth',2);
        legendInfo = union(legendInfo, 'RABB-w/o locking', 'stable'); 
    end
    if(ismember('m0', ABB))
        plot(disABBm0,'-*', 'markersize',12,'linewidth',2);
        legendInfo = union(legendInfo, 'LRBFGS-ABB: 0', 'stable'); 
    end
    if(ismember('m4', ABB))
        plot(disABBm4,'-*', 'markersize',12,'linewidth',2);
        legendInfo = union(legendInfo, 'LRBFGS-ABB: 4', 'stable');
    end
    if(ismember('m16', ABB))
        plot(disABBm16,'-*', 'markersize',12,'linewidth',2);
        legendInfo = union(legendInfo, 'LRBFGS-ABB: 16', 'stable');
    end
end
l1 = legend(legendInfo);
set(l1,'Interpreter','Latex');
set(gca,'FontSize',13)
xlabel('iterations','FontSize',16,'Interpreter','Latex');
ylabel('$dist(\mu, X_t)$','FontSize',16,'Interpreter','Latex');
axis([0 60 1e-16 1e+3]);




% legendInfo1 = {'RBB1', 'RBB1 Final', 'RBB2', 'RBB2 Final', '$\mathrm{RABB}_{\min}$', '$\mathrm{RABB}_{\min}$ Final'};
% 
% 
% eig_color = [255,105,180]./255;
% % eig_color = [30,144,255]./255;
% figure(1);
% EigInfo = importdata('EigInfo.mat');
% h_min = plot(EigInfo.ix, EigInfo.imin_eig, '.', 'Color', eig_color, 'markersize',4, 'linewidth', 2);
% hold on;
% h0_min = plot(EigInfo.x, EigInfo.min_eig,'o', 'Color', eig_color, 'markersize',3, 'linewidth', 2);
% h_max = plot(EigInfo.ix, EigInfo.imax_eig, '.','Color', eig_color, 'markersize',4, 'linewidth', 2);
% h0_max = plot(EigInfo.x, EigInfo.max_eig,'o', 'Color', eig_color, 'markersize',3, 'linewidth', 2);
% 
% % h_minWen = plot(EigInfo.ix, EigInfo.imin_eigWen, 'r.', 'markersize',4, 'linewidth', 2);
% % h0_minWen = plot(EigInfo.x, EigInfo.min_eigWen,'ro', 'markersize',3, 'linewidth', 2);
% % h_maxWen = plot(EigInfo.ix, EigInfo.imax_eigWen, 'r.', 'markersize',4, 'linewidth', 2);
% % h0_maxWen = plot(EigInfo.x, EigInfo.max_eigWen,'ro', 'markersize',3, 'linewidth', 2);
% % set(get(get(h_minWen,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% % set(get(get(h0_minWen,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% % set(get(get(h_maxWen,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% % set(get(get(h0_maxWen,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% 
% set(get(get(h_min,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% set(get(get(h0_min,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% set(get(get(h_max,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% set(get(get(h0_max,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% bb = lines(10);
% plot(1./ihBB1,'*', 'Color', bb(5, :), 'markersize',12,'linewidth',2);
% plot(1./hBB1,'o', 'Color', bb(1, :), 'markersize',12,'linewidth',2);
% plot(1./ihBB2,'x', 'markersize',12,'linewidth',2);
% plot(1./hBB2,'d', 'markersize',12,'linewidth',2);
% plot(1./ihABB(1:end-3),'+', 'Color', bb(3, :), 'markersize',12,'linewidth',2);
% plot(1./hABB(1:end-3),'s', 'Color', bb(2, :), 'markersize',12,'linewidth',2);
% axis([1 60 0 0.5])
% 
% 
% 
% disLRBFGSBB1((disLRBFGSBB1 == 0)) = 1.1e-16;
% disBB((disBB1 == 0)) = 1.1e-16;
% disBB((disBB2 == 0)) = 1.1e-16;
% disABB((disABB == 0)) = 1.1e-16;


% figure(1);
% l1 = legend(legendInfo1);
% set(l1,'Interpreter','Latex');
% set(gca,'FontSize',13)
% xlabel('iterations','FontSize',16,'Interpreter','Latex');
% ylabel('$1/\alpha_k$','FontSize',16,'Interpreter','Latex');

 

% figure(3)
% h_min = plot(EigInfo.ix, EigInfo.imin_eig, '.', 'Color', eig_color, 'markersize',4, 'linewidth', 2);
% hold on;
% h0_min = plot(EigInfo.x, EigInfo.min_eig,'o', 'Color', eig_color, 'markersize',3, 'linewidth', 2);
% h_max = plot(EigInfo.ix, EigInfo.imax_eig, '.','Color', eig_color, 'markersize',4, 'linewidth', 2);
% h0_max = plot(EigInfo.x, EigInfo.max_eig,'o', 'Color', eig_color, 'markersize',3, 'linewidth', 2);
% set(get(get(h_min,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% set(get(get(h0_min,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% set(get(get(h_max,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% set(get(get(h0_max,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% bb = lines(10);
% plot(1./ihBB1,'*', 'Color', bb(5, :), 'markersize',12,'linewidth',2);
% plot(1./hBB1,'o', 'Color', bb(1, :), 'markersize',12,'linewidth',2);
% axis([1 60 0 0.5]);
% 
% 
% 
% figure(4)
% h_min = plot(EigInfo.ix, EigInfo.imin_eig, '.', 'Color', eig_color, 'markersize',4, 'linewidth', 2);
% hold on;
% h0_min = plot(EigInfo.x, EigInfo.min_eig,'o', 'Color', eig_color, 'markersize',3, 'linewidth', 2);
% h_max = plot(EigInfo.ix, EigInfo.imax_eig, '.','Color', eig_color, 'markersize',4, 'linewidth', 2);
% h0_max = plot(EigInfo.x, EigInfo.max_eig,'o', 'Color', eig_color, 'markersize',3, 'linewidth', 2);
% set(get(get(h_min,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% set(get(get(h0_min,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% set(get(get(h_max,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% set(get(get(h0_max,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% bb = lines(10);
% plot(1./ihBB2,'*', 'Color', bb(5, :), 'markersize',12,'linewidth',2);
% plot(1./hBB2,'o', 'Color', bb(1, :), 'markersize',12,'linewidth',2);
% axis([1 60 0 0.5])
% 
% 
% figure(5)
% h_min = plot(EigInfo.ix, EigInfo.imin_eig, '.', 'Color', eig_color, 'markersize',4, 'linewidth', 2);
% hold on;
% h0_min = plot(EigInfo.x, EigInfo.min_eig,'o', 'Color', eig_color, 'markersize',3, 'linewidth', 2);
% h_max = plot(EigInfo.ix, EigInfo.imax_eig, '.','Color', eig_color, 'markersize',4, 'linewidth', 2);
% h0_max = plot(EigInfo.x, EigInfo.max_eig,'o', 'Color', eig_color, 'markersize',3, 'linewidth', 2);
% set(get(get(h_min,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% set(get(get(h0_min,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% set(get(get(h_max,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% set(get(get(h0_max,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% bb = lines(10);
% plot(1./ihABB,'*', 'Color', bb(5, :), 'markersize',12,'linewidth',2);
% plot(1./hABB,'o', 'Color', bb(1, :), 'markersize',12,'linewidth',2);
% axis([1 60 0 0.5])
% 



% figure(6)
% h = semilogy(1,1);
% hold on;
% set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% plot(gradsBB1,'-o',  'markersize',12,'linewidth',2);
% plot(gradsBB2,'-o', 'markersize',12,'linewidth',2);
% plot(gradsABB,'-s', 'markersize',12,'linewidth',2);
% plot(gradsLRBFGSBB1,'-*','markersize',12,'linewidth',2);
% plot(gradsLRBFGSBB2,'-*','markersize',12,'linewidth',2);
% plot(gradsLRBFGSABB,'-*','markersize',12,'linewidth',2);
% legendInfo3 = {'RBB1', 'RBB2', '$\mathrm{RABB}_{\min}$', 'LRBFGS-BB1', 'LRBFGS-BB2', 'LRBFGS-$\mathrm{ABB}_{\min}$'};
% l3 = legend(legendInfo3);
% set(l3,'Interpreter','Latex');
% set(gca,'FontSize',13)
% xlabel('iterations','FontSize',16,'Interpreter','Latex');
% ylabel('$dist(\mu, X_t)$','FontSize',16,'Interpreter','Latex');
% axis([0 60 1e-16 1e+3]);
end




