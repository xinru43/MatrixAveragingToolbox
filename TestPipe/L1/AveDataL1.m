function AveDataL1(inData, dt, SavePath, FileName)

Ndata = length(inData);

%ComTime = zeros(1, stabletime);

AveComTime = 0.0;
ComTime0 = 0.0;
nf = 0;
ng = 0;
nR = 0;
nV = 0;
nVp = 0;
nH = 0;
for times = 1 : Ndata
    AveComTime = AveComTime + inData{times}.AveComTime;
    AveComTimeStat(times) = inData{times}.AveComTime;
    ComTime0 = ComTime0 + inData{times}.ComTime0;
    nf = nf + inData{times}.nf;
    ng = ng + inData{times}.ng;
    nR = nR + inData{times}.nR;
    nV = nV + inData{times}.nV;
    nVp = nVp + inData{times}.nVp;
    nH = nH + inData{times}.nH;
    Eig{times} = inData{times}.Eig;
end
AveComTime = AveComTime/Ndata;
nf = nf/Ndata;
ng = ng/Ndata;
nR = nR/Ndata;
nV = nV/Ndata;
nVp = nVp/Ndata;
nH = nH/Ndata;
[dis, dis_t] = ave_distance(inData, Ndata, dt); % average distance with time
dis_iter = ave_dis_iter(inData, Ndata); % average distance with iterations

[funs, funs_t] = ave_f(inData, Ndata, dt); % average distance with time
funs_iter = ave_f_iter(inData, Ndata);

grad = ave_grad(inData, Ndata); % average |grad|
T_iter = ave_T_iter(inData, Ndata);
save([SavePath FileName '.mat'], ...
    'dis', 'grad', 'dis_t', 'dis_iter', 'AveComTime', 'ComTime0', 'AveComTimeStat', ...
    'nf', 'ng', 'nR', 'nV', 'nVp', 'nH', 'Eig', 'T_iter', 'funs_iter', 'funs', 'funs_t');
end



%% compute the average of gradient norm
function X_grad = ave_grad(X,times)
Xtemp = [X{:}];
X_max_iter = max([Xtemp.iter])+1;
X_min_iter = min([Xtemp.iter])+1;
X_grad = zeros(1,X_max_iter);
for i = 1: times
    Grad = X{i}.G(end);
    X_grad = X_grad + [X{i}.G Grad*ones(1,length(X_grad) - length(X{i}.G))];
end
X_grad = X_grad/times;
end



%% average the distance along with iterations
function X_dis = ave_dis_iter(X,times)
Xtemp = [X{:}];
X_max_iter = max([Xtemp.iter])+1;
X_dis = zeros(1,X_max_iter);
for i = 1: times
    dis = X{i}.dis(end);
    X_dis = X_dis + [X{i}.dis dis*ones(1,length(X_dis) - length(X{i}.dis))];
end
X_dis = X_dis/times;
end




%% average the distance along with iterations
function X_f = ave_f_iter(X, times)
Xtemp = [X{:}];
X_max_iter = max([Xtemp.iter])+1;
X_f = zeros(1, X_max_iter);
for i = 1: times
    f = X{i}.F(end);
    X_f = X_f + [X{i}.F f * ones(1,length(X_f) - length(X{i}.F))];
end
X_f = X_f/times;
end





%% average the stepsize along with iterations
function X_truestepsize = ave_stepsize_iter(X,times)
Xtemp = [X{:}];
X_max_iter = max([Xtemp.iter])+1;
X_truestepsize = zeros(1,X_max_iter);
for i = 1:times
    truestepsize = X{i}.truestepsize(end);
    X_truestepsize = X_truestepsize + [X{i}.truestepsize truestepsize*ones(1,length(X_truestepsize) - length(X{i}.truestepsize))];
end
X_truestepsize = X_truestepsize/times;
end



%% average df/f along with iterations
function X_relf = ave_relf_iter(X,times)
Xtemp = [X{:}];
X_max_iter = max([Xtemp.iter])+1;
X_relf = zeros(1,X_max_iter);
for i = 1:times
    relf = X{i}.relf(end);
    X_relf = X_relf + [X{i}.relf relf*ones(1,length(X_relf) - length(X{i}.relf))];
end
X_relf = X_relf/times;
end



%% compute the average of distance
function [AveDis,t] = ave_distance(X, times, dt)
Xtemp = [X{:}];
X_max_time = max([Xtemp.T]);

ngrid = floor(X_max_time/dt);

AveDis = zeros(1,ngrid);
t = 0 : dt : (ngrid-1) * dt;
for i = 1:times
    Dis = X{i}.dis(end);%------------
    Grid{i} = Dis * ones(1,ngrid);
    for k = 1 : length(X{i}.T)
        ind(k) = floor(X{i}.T(k)/dt);
    end
    Grid{i}(1:ind(1)) = X{i}.dis(1);
    for j = 1 : length(X{i}.dis)-1
        Grid{i}(ind(j)+1 : ind(j+1)) = X{i}.dis(j);
    end
    AveDis = AveDis + Grid{i};
end
AveDis = AveDis/times;
end




%% compute the average of function
function [Avef, t] = ave_f(X, times, dt)
Xtemp = [X{:}];
X_max_time = max([Xtemp.T]);

ngrid = floor(X_max_time/dt);
Avef = zeros(1, ngrid);
t = 0 : dt : (ngrid-1) * dt;
for i = 1 : times
    f = X{i}.F(end);
    Grid{i} = f * ones(1,ngrid);
    for k = 1 : length(X{i}.T)
        ind(k) = floor(X{i}.T(k)/dt);
    end
    Grid{i}(1:ind(1)) = X{i}.F(1);
    for j = 1 : length(X{i}.F)-1
        Grid{i}(ind(j)+1 : ind(j+1)) = X{i}.F(j);
    end
    Avef = Avef + Grid{i};
end
Avef = Avef/times;
end




function X_T = ave_T_iter(X,times)
Xtemp = [X{:}];
X_max_iter = max([Xtemp.iter])+1;
X_T = zeros(1,X_max_iter);
for i = 1: times
    T = X{i}.T(end);
    X_T = X_T + [X{i}.T T*ones(1,length(X_T) - length(X{i}.T))];
end
X_T = X_T/times;
end