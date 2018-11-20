function AveDataFP(inData, dt, stabletime, SavePath, FileName)

Ndata = length(inData);

ComTime = zeros(1, stabletime);
ComTime0 = 0.0;

for times = 1 : Ndata
    ComTime = ComTime + inData{times}.ComTime;
    ComTime0 = ComTime0 + inData{times}.ComTime0;
end
ComTime = ComTime/Ndata;

[dis, dis_t] = ave_distance(inData, Ndata, dt); % average distance with time
dis_iter = ave_dis_iter(inData, Ndata); % average distance with iterations
grad = ave_grad(inData, Ndata); % average |grad|
save([SavePath FileName '.mat'], ...
    'dis', 'grad', 'dis_t', 'dis_iter', 'ComTime', 'ComTime0');
end




%%
%%
%%
%%


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
