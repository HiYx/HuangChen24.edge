clc, clear
%% ================= Simulation of LODCO Algorithm =================

%% 基本参数设置
k = 1e-28;                        % 有效开关电容
tau = 0.002;                      % 时间片长度(s)
tau_d = 0.002;                    % 计算任务执行时间的deadline(s)
phi = 0.002;                      % 任务丢弃的惩罚项权重(s)
omega = 1e6;                      % 服务器带宽(Hz)
sigma = 1e-13;                    % 接收端的噪声功率(W)
p_tx_max = 1;                     % 移动设备最大传输功率(W)
f_max = 1.5e9;                    % 移动设备最大CPU时钟周期频率(Hz)
E_max = 0.002;                    % 电池允许的最大放电量(J)
L = 1000;                         % 一项计算任务的大小(bit)
X = 737.5;                        % 移动设备执行1bit任务所需的时钟周期个数
W = L*X;                          % 移动设备本地执行一项计算任务所需的时钟周期个数
E_H_max = 48e-6;                  % 收集的能量服从的均匀分布上限(J)
g0 = power(10,-4);                % 路径损失常数(dB转化之后的数值比)
d = 50;                           % 服务器和移动设备之间的相对距离(m)

%% 变量控制
T = 200;                         % 时间片个数
E_min = 0.02e-3;                  % 能量使用下界(J)
V = 1e-5;                         % LODCO中penalty项的权重(J^2/second)
rho = 0.6;                        % 计算任务抵达的概率

% 实际能耗上限
E_max_hat = min(max(k*W*(f_max^2), p_tx_max*tau), E_max);
theta = E_max_hat + V*phi/E_min;        % 扰动参数


%% 中间变量存储
B = zeros(T, 1);                        % 实际电量
B_hat = zeros(T, 1);                    % 虚拟电量
e = zeros(T, 1);                        % 能量收集
indicator = zeros(T, 1);                % 1/2/3/4分别对应local,remote,dropped,无任务产生
f = zeros(T, 1);                        % 移动设备本地执行的频率
p = zeros(T, 1);                        % 移动设备卸载执行的传输功率
cost = zeros(T, 3);                     % 每一列分别对应local,remote,最后选择的模式下的execution cost
E = zeros(T, 3);                        % 每一列分别对应local,remote,最后选择的模式下的能耗

% 关闭fslove的输出
opt = optimset('Display', 'off');

t = 1;
while t <= T
    %% 阶段初始化
    % 以伯努利分布产生计算任务
    zeta = binornd(1, rho);
    % 产生虚拟电量队列值
    B_hat(t) = B(t) - theta;

    %% 求解optimal energy harvesting e* (不管是否有任务产生，能量收集都不能停！)
    % 产生E_H_t
    E_H_t = unifrnd(0, E_H_max);
    if B_hat(t) <= 0
        e(t) = E_H_t;
    end
    
    if zeta == 0
        % 没有计算任务产生
        indicator(t) = 4;
        % f(t) = 0; p(t) = 0;  默认即为0
    else
        % 产生信道功率增益
        h = exprnd(g0/power(d,4));

        %% 求解P_ME
        f_L = max(sqrt(E_min/(k*W)), W/tau_d);
        f_U = min(sqrt(E_max/(k*W)), f_max);
        if f_L <= f_U
            % P_ME有解
            f0 = power(V/(-2*B_hat(t)*k), 1/3);
            if f0 > f_U
                f(t) = f_U;
            elseif f0 >= f_L && f0 <= f_U && B_hat(t) < 0
                f(t) = f0;
            elseif f0 < f_L
                f(t) = f_L;
            end
            % 计算local的execution delay
            cost(t, 1) = W / f(t);
            % 计算此时的能耗
            E(t, 1) = k * W * (f(t)^2);
            if E(t, 1) >= B(t)
                disp(['本地执行电量不足![此时t为', num2str(t), ']']);
                % 当电量不足或者本地执行子问题无解的时候，将该子问题的目标函数值设为inf
                J_m = inf;
            else
                % 计算此时的J_m，即本地执行子问题的目标函数值
                J_m = -B_hat(t)*E(t,1) + V*cost(t,1);
            end

        else
            disp(['P_ME无解![此时t为', num2str(t), ']']);
            J_m = inf;
        end

        %% 求解P_SE
        E_tmp = sigma*L*log(2) / (omega*h);
        p_L_taud = (power(2, L/(omega*tau_d)) - 1) * (sigma/h);
        if E_tmp >= E_min
            p_L = p_L_taud;
        else
            % 计算p_Emin
            y = @(x)x*L-omega*log2(1+h*x/sigma)*E_min;
            %p_Emin = double(vpa(solve(y, 1)));
            tmp = fsolve(y, [0.001, 1], opt);
            p_Emin = max(tmp);
            p_L = max(p_L_taud, p_Emin);
        end
        if E_tmp >= E_max
            p_U = 0;
        else
            % 计算p_Emax
            y = @(x)x*L-omega*log2(1+h*x/sigma)*E_max;
            p_Emax = max(fsolve(y, [0.001, 100]));
            p_U = min(p_tx_max, p_Emax);
        end
        if p_L <= p_U
            % P_SE有解
            % 计算p0
            tmp = B_hat(t);
            syms x
            y = tmp*log2(1+h*x/sigma) + h*(V-tmp*x)/(log(2)*(sigma+h*x));
            p0 = double(vpasolve(y));
            if p_U < p0
                p(t) = p_U;
            elseif p_L > p0 && B_hat(t) < 0
                p(t) = p_L;
            elseif p_L <= p0 && p_U >= p0 && B_hat(t) < 0
                p(t) = p0;
            end
             % 计算achievable rate
            r = calAchieveRate(h, p(t), omega, sigma);
            % 计算此时的execution delay
            cost(t, 2) = L / r;
            % 计算此时的能耗
            E(t, 2) = p(t) * L / r;
            if E(t, 2) >= B(t)
                disp(['卸载执行电量不足![此时t为', num2str(t), ']']);
                J_s = inf;
            else
                % 计算此时的J_s卸载执行子问题的目标函数值
                J_s = -B_hat(t)*E(t,2) + V*cost(t,2);
            end
        else
            disp(['P_SE无解![此时t为', num2str(t), ']']);
            J_s = inf;
        end

        %% 选取最佳模式
        [~, mode] = min([J_m, J_s, phi]);
        indicator(t) = mode;
    end

    % 计算实际选择的模式下的execution cost和能耗
    if indicator(t) == 1
        cost(t, 3) = cost(t, 1);
        E(t, 3) = E(t, 1);
    elseif indicator(t) == 2
        cost(t, 3) = cost(t, 2);
        E(t, 3) = E(t, 2);
    % 剩下两种情况能耗为0，不变
    elseif indicator(t) == 3
        cost(t, 3) = phi;
    else
        % 即没有任务产生的情况
        cost(t, 3) = 0;
    end
    % 电量迭代
    B(t+1) = B(t) - E(t, 3) + e(t);
    % 时间片迭代
    t = t + 1;
end

%% results
% 产生任务的时间片个数
num = T - sum(indicator == 4);
disp(['当前任务到达的频率为：', num2str(num/T)]);
mode1 = sum(indicator==1) / (num);
disp(['任务本地执行的比率：', num2str(mode1)]);
mode2 = sum(indicator==2) / (num);
disp(['任务卸载执行的比率:', num2str(mode2)]);
mode3 = sum(indicator==3) / (num);
disp(['任务被抛弃的比率：', num2str(mode3)]);


figure
plot(1:T, B(1:T, :))
hold on
plot(1:T, repmat(theta + E_H_max, [T, 1]), '-')
title('Fig. 2. Battery energy level of each mobile device vs. time.')
xlabel('time slot')
ylabel('battery energy level $B_t$ of each mobile device', 'Interpreter','latex')

N=1;
accumulated = 0;
average_cost = zeros(T, N);

chosen_mode= indicator;
figure
% for i = 1: N
%     % draw for each mobile device
%     request_num = 0;
%     for t = 1: T
%         accumulated = accumulated + cost(t, i);
%         if chosen_mode(t, i) ~= 4
%             % there exists task request
%             request_num = request_num + 1;
%         end
%         average_cost(t, i) = accumulated / request_num;
%     end
%     %plot(1:T, average_cost(:, i));
%     
%     hold on
% end
% plot(1:T, mean(average_cost(:, i),2));
plot(1:T, mean(cost,2),'--b')
title('Fig. 3. Average QoE-Cost of all mobile device vs. time')
xlabel('time slot')
ylabel('average execution cost $\frac{1}{T} \sum_{t=0}^{T-1} cost^t$ of each mobile device', 'Interpreter','latex')



% Fig. 4. Average energy level of each mobile device.
figure
all_ave = []
for i =1:N
    all_ave(i) = mean( B(:,i));
end
bar(1:N, all_ave)
hold on
%plot(1:N, repmat(theta + E_H_max, [N, 1]), '-') %  实际能耗上限
title('Average energy level of each mobile device. ')
xlabel('mobile N=10')
ylabel('battery energy level', 'Interpreter','latex')



%Fig. 5. The ratio of each chosen modes vs. time.
% 移动端选择的三种模式的比例
average_ratio = zeros(T, 3);
mobile_exe = 0; server_exe = 0; drop = 0;
request_num = 0;
i = 1;          % we simply choose the first mobile device
for t = 1: T
    if cost(t, i) == 0
        if request_num < 0.9
            average_ratio(t, :) = [0, 0, 0];
        else
            average_ratio(t, :) = [mobile_exe, server_exe, drop] / request_num;
        end
        continue
    else
        request_num = request_num + 1;
        if chosen_mode(t, i) == 1
            mobile_exe = mobile_exe + 1;
        elseif chosen_mode(t, i) == 2
            server_exe = server_exe + 1;
        else
            drop = drop + 1;
        end
    end
    average_ratio(t, :) = [mobile_exe, server_exe, drop] / request_num;
end
figure
plot(1:T, average_ratio(:, 1));
hold on
plot(1:T, average_ratio(:, 2));
hold on
plot(1:T, average_ratio(:, 3));
legend('mobile execution', 'MEC server execution', 'drop')
title('Fig. 5. The ratio of each chosen modes vs. time.')
xlabel('time slot')
%ylabel('average  ratio of chosen modes $\frac{1}{T} \sum_{t=0}^{T-1} \{I_m^t, I_s^t, I_d^t\}$ of the i-th mobile device', 'Interpreter','latex')
ylabel('average  ratio of chosen modes  of the i-th mobile device', 'Interpreter','latex')





% Fig. 6. Average ratio of offloading tasks by different algorithms.
% 在图 6 中，在第二个优化目标（最大卸载计算任务的数量）上，
% 比较了基于 LODCO 的贪心策略遗传算法与基准算法（基于 LODCO 的贪心算法）的性能。
% 可以看到，基于 LODCO 的贪心策略遗传算法得到的卸载任务的平均比率（95.0698%）大于基于 LODCO 的贪心算法得到的卸载任务比率（92.8549%5），
% 而 这两个比率都大于 LODCO 算法得到的卸载任务的平均比率（即 84.4401%）。
% 结果证明，我们的算法可以获得比基准算法更好的小区容量。 
figure
plot(1:T, average_ratio(:, 2));
hold on
title('load')
xlabel('time slot')
ylabel('server_exe load ratio', 'Interpreter','latex')