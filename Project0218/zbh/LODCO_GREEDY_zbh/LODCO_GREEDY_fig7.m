% the LODCO-based Greedy algorithm.贪心算法的LODCO卸载决策

clc, clear






%% 距离的影响  如图7所示
%  描述了随机移动设备和随机 MEC 服务器之间的最大距离对小区容量的影响。
%  可以看出，随着最大距离的增加，卸载任务的平均比例逐渐减小。
%  当距离任意远时（在最大距离约束下），卸载计算任务的数量会非常少。
%  原因是信道功率增益随着每个移动设备和每个MEC服务器之间的距离而增长，
%  这将导致更大的能量消耗和更长的执行延迟。 
%  因此，更多的移动设备选择在本地执行计算任务。 

%Fig. 7. Average ratio of offloading tasks vs. maximum distance.
max_distance= [40,50,60,70,80,90,100];
ratio=[];
for i =1:7
ratio(i,:) = func_lodco(max_distance(i),0.002,0.6,200);
end
figure
%plot(max_distance, mean(ratio,2));
plot(max_distance,ratio(:,200));
hold on
title('load')
xlabel('time slot')
ylabel('Average ratio of offloading tasks1 ', 'Interpreter','latex')



%% 任务量的影响，如图8所示
% 随着γ的增加，属于[0,1]，卸载任务的平均比例随着减速率逐渐增加，最终收敛到特定值（即97.5731%）。 
% 原因是基于 LODCO 的贪心策略遗传算法是基于 γ-贪婪策略的，即更大的 γ 将带来更大的选择卸载模式的概率。 
% Fig. 8. Average ratio of offloading tasks vs. phi  and rho, respectively

phi= [0.001,0.002,0.003,0.004,0.005,0.006,0.007];
ratio=[];
for i =1:7
ratio(i,:) = func_lodco(50,phi(i),0.6,200);
end
figure
%plot(phi, mean(ratio,2));
plot(phi,ratio(:,200));
hold on
title('load')
xlabel('time slot')
ylabel('Average ratio of offloading tasks2 ', 'Interpreter','latex')


rho= [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1];
ratio=[];
for i =1:10
ratio(i,:) = func_lodco(50,0.002,rho(i),200);
end
figure
plot(rho, mean(ratio,2));
%plot(rho,ratio(:,200));
hold on
title('load')
xlabel('time slot')
ylabel('Average ratio of offloading tasks3 ', 'Interpreter','latex')


%输入：最远距离 rho，任务丢弃的惩罚项权重 phi，计算任务抵达的概率 rho，运行的时间片长度T
% 默认 max_distance=50（取值范围40-110）  ,phi=0.002（取值范围0.001-0.007）,rho（取值范围0.1-1）,
% 默认 T =500
function ratio = func_lodco(max_distance,phi,rho,T)
opt = optimset('Display', 'none');
%% basic parameter settings (had better not change those paras) 基本参数设置（最好不要修改）
k = 1e-28;                % effective switched capacitance (a constant decided by the chip architecture)有效开关电容（设备固有属性）
tau = 0.002;              % the length of time slot (in second)时间片长度(0.002s)
%phi = 0.002;              % the cost of task dropping (in second) 任务丢弃的惩罚项权重(0.002s)
omega = 1e6;              % the bandwidth of MEC server (in Hz)服务器带宽(Hz)
sigma = 1e-13;            % the noise power of the receiver (in W)接收端的噪声功率(W)
p_tx_max = 1;             % the maximum transmit power of mobile device (in W)移动设备最大传输功率(W)
f_max = 1.5e9;            % the maximum CPU-cycle frequency of mobile device (in Hz)移动设备最大CPU时钟周期频率(Hz)
E_max = 0.002;            % the maximum amout of battery output energy (in J)电池允许的最大放电量(J)
L = 1000;                 % the input size of the computation task (in bit)一项计算任务的大小(bit)
X = 737.5;                % the number of CPU cycles needed on processing one bit of task移动设备执行1bit任务所需的时钟周期个数
W = L * X;                % the number of CPU cycles needed on processing one task移动设备本地执行一项计算任务所需的时钟周期个数737500个
E_H_max = 48e-6;          % the upper bound of the energy arrive at the mobile device (in J)收集的能量服从的均匀分布上限(J)
p_H = E_H_max / (2*tau);  % the average Energy Harvesting (EH) power (in W) 平均能量收集(EH)功率0.0120W
g0 = power(10, -4);       % the path-loss constant 路径损失常数(dB转化之后的数值比)
d0 = 1;                   % the relative distance between each mobile device and each MEC server 服务器和移动设备之间的相对距离(m)

%% parameter control
N = 10;                   % the number of mobile devices  移动设备数目
M = 8;                    % the number of MEC servers MEC服务器个数
%T = 500;                 % the number of time slot (a.k.a. the size of the time horizon) 时间片个数
tau_d = 0.002;            % execution deadline (in second)计算任务执行时间的deadline(s)
d = 50;                   % the distance between the mobile device and the MEC server (in meter)移动设备与MEC服务器之间的距离(米)
E_min = 0.02e-3;          % the minimum amout of battery output energy (in J) 能量使用下界(J)2%的Emax能量
V = 1e-5;                 % the weight of penalty (the control parameter introduced by Lyapunov Optimization)  LODCO中penalty项的权重(J^2/second)
%rho = 0.6;                % the probability that the computation task is requested 计算任务抵达的概率
max_connects = 4;         % the maximum number of processible mobile devices for each MEC server ($ \frac{f_s^{max} \tau}{L X} $) 每台MEC服务器的最大可处理移动设备数
min_distance = 1;        % the minimum distance from mobile device to MEC server移动设备与MEC服务器之间的距离(米)
%max_distance = 50;        % the maximum distance from mobile device to MEC server移动设备与MEC服务器之间的距离(米)


% the lower bound of perturbation parameter   实际能耗上限
E_max_hat = min(max(k * W * (f_max)^2, p_tx_max * tau), E_max);
theta = E_max_hat + V * phi / E_min;  % 扰动参数

%% allocate storage for valuable results
B = zeros(T, N);                   % the battery energy level (in J) N个移动设备的实际电量
B_hat = zeros(T, N);               % the virtual battery energy level ($B_hat = B - theta$) N个移动设备的虚拟电量
e = zeros(T, N);                   % the amout of the harvested and stored energy (in J) N个移动设备的能量收集
chosen_mode = zeros(T, N);         % {1: local, 2: remote, 3: drop, 4: no task request} 用1,2,3,4分别表示本地执行、卸载执行、drop(前三个均意味着有任务到达)以及没有任务到达
chosen_server = zeros(T, N);       % record the index of chosen server for each mobile device if its choice is MEC server execution如果选择的是MEC服务器执行，则记录每个移动设备选择的服务器的索引
f = zeros(T, N);                   % the CPU-cycle frequency of local execution (in Hz) N个移动设备本地执行的频率
p = zeros(T, N);                   % the transmit power of computation offloading (in W)  N个移动设备卸载执行的传输功率


mobile_exe_cost = zeros(T, N);     % the mobile execution cost (delay) (in second) N个移动设备的local execution delay
server_exe_cost = zeros(T, N);     % the MEC server execution cost (delay) (in second) N个移动设备的remote execution delay
final_chosen_cost = zeros(T, N);   % the final execution delay under currently chosen modes (in second) N个移动设备最后选择的模式下的execution cost

mobile_exe_E = zeros(T, N);        % the energy consumption for mobile execution (in J)  N个移动设备本地执行的能耗
server_exe_E = zeros(T, N);        % the energy consumption for MEC server execution (in J) N个移动设备卸载执行的能耗
final_chosen_E = zeros(T, N);      % the energy consumption of the final chosen modes (in J)  N个移动设备在最后选择的模式下的实际能耗

mode_num = zeros(T,3);                  % 每一列分别表示每轮中本地执行、远程执行及任务丢弃的比率(分母为N减去没有任务到达的)

%% simulation begin
t = 1;
while t <= T
    disp(['===> Time slot #', num2str(t), ' <==='])
    
    %% allocate storage for mode-chosen
    device_server_pairs = [];      % each column represents i, j, J_s^{\star}(i, j), respectively % 用一个“不定长的3列”的矩阵描述【"移动设备i，MEC服务器j，i到j的最小卸载执行子问题目标函数值】"
    remained_connects = max_connects * ones(M, 1);    % the available connections of MEC servers 存储每一个MEC服务器连接到的移动设备个数
    J_s = zeros(N, M);             % the matrices for J_m and J_s values 分别保存当前移动设备到各MEC服务器的J_s值，J_s是卸载执行子问题的目标函数值，而非延迟！！
    J_m = zeros(N, 1);              % 分别存储每一个移动设备的本地执行的J_m，供二次决策时使用

    p_mat = zeros(N, M);                              % the matrix for transmit power 当前移动设备到各MEC服务器的传送的能量
    server_cost_mat = zeros(N, M);                    % the matrix for MEC server execution cost MEC服务器执行cost矩阵
    server_E_mat = zeros(N, M);                       % the matrix for energy consumption of MEC server MEC服务器能耗矩阵
    
    %% initialization 初始化
    % generate the virtual battery energy level 产生虚拟电量队列值
    B_hat(t, :) = B(t, :) - theta;
    % generate the channel power gain (from each mobile device to each MEC sever)产生信道功率增益(从每个移动设备到每个MEC服务器)
    distances = unifrnd(min_distance, max_distance, N, M);% 随机产生服务器和移动设备的距离(限定在10 ~ 50之内)
    gamma = exprnd(1, N, M);% 服从lambda=1的指数分布的小尺度衰落信道功率收益
    h_mat = g0 * gamma .* power(d0 ./ distances, 4); % 从任意移动设备到任意服务器的信道功率增益
    
    %% step 1: for each mobile device, choose the initial mode 对每一个移动设备，阶段初始化
    for i = 1: N
        %disp(['Mobile device #', num2str(i)])
        
        %% step 1.1: get the optimal energy harvesting no matter whether task is requested
        %% 求解optimal energy harvesting e* (不管是否有任务产生，能量收集都不能停！)
         % 产生E_H_t，unifrnd生成被0和E_H_max指定上下端点[0,E_H_max]的连续均匀分布的随机数组E_H_t。
        E_H_t = unifrnd(0, E_H_max);% 虚拟电量队列值
        if B_hat(t, i) <= 0
            e(t, i) = E_H_t;
        end
        
        %% step 1.2: get the (initial) optimal computation offloading strategy (I_m(i), I_s(i, :), I_d(i), f(t, i), p(t, i))
        % generate the task request
        % 以伯努利分布产生计算任务
        zeta = binornd(1, rho);
        if zeta == 0
            % chosen mode has to be 4
             % 没有计算任务产生
            %disp('no task request generated!')
            chosen_mode(t, i) = 4;
            % f(t, i) = 0; p(t, i) = 0;    默认即为0
        else
            % chosen_mode is chosen from {1, 2, 3}
            % 计算任务产生
            %disp('task request generated!')
            
            %% step 1.2.1: solve the optimization problem $\mathcal{P}_{ME}$ (f(t, i) > 0)
            %% 求解P_ME
            % calculate f_L and f_U
            f_L = max(sqrt(E_min / (k * W)), W / tau_d);
            f_U = min(sqrt(E_max / (k * W)), f_max);
            if f_L <= f_U
                % the sub-problem is feasible
                %% disp
%                disp('mobile execution ($\mathcal{P}_{ME}$) is feasible!')
                % P_ME有解
                if B_hat(t, i) < 0
                    f_0 = (V / (-2 * B_hat(t, i) * k))^(1/3);
                else
                    % complex number may exist, which may lead to error
                    % 存在复数的话可能导致错误
                    f_0 = -(V / (2 * B_hat(t, i) * k))^(1/3);
                end
                
                if (f_0 > f_U && B_hat(t, i) < 0) || (B_hat(t, i) >= 0)
                    f(t, i) = f_U;
                elseif f_0 >= f_L && f_0 <= f_U && B_hat(t, i) < 0
                    f(t, i) = f_0;
                elseif f_0 < f_L && B_hat(t, i) < 0
                    f(t, i) = f_L;
                end
                % check whether f(t, i) is zero
                if f(t, i) == 0
                    disp('Something wrong! f is 0!')
                end
                
                % calculate the delay of mobile execution
                % 计算此时的execution delay
                mobile_exe_cost(t, i) = W / f(t, i);
                % calculate the energy consumption of mobile execution
                % 计算此时的能耗
                mobile_exe_E(t, i) = k * W * (f(t, i)^2);
                % calculate the value of optimization goal
                % 计算优化目标值
                %  J_m代表的是子问题的目标函数值！是能耗和延迟统一变换后的结果，而非仅仅是延迟！
                J_m(i) = -B_hat(t, i) * k * W * (f(t, i)^2 + V * W / f(t, i));
            else
                % the sub-problem is not fasible because (i) the limited 
                % computation capacity or (ii) time cosumed out of deadline 
                % or (iii) the energy consumed out of battery energy level
                % If it is not feasible, it just means that we cannot choose 
                % 'I_m(i)=1'. It dosen't mean that the task has to be dropped.
                % P_ME无解
                %% disp
%                disp('mobile execution ($\mathcal{P}_{ME}$) is not feasible!')
                f(t, i) = 0;
                mobile_exe_cost(t, i) = 0;
                mobile_exe_E(t, i) = 0;
                % 'I_m(i)=1' can never be chosen if mobile execution goal is inf
                J_m(i) = inf;
            end
            
            %% step 1.2.2: solve the optimization problem $\mathcal{P}_{SE}$ (p(t, i) > 0)
            %% 求解P_SE
            % calculate J_s(i, j) from mobile device i to each MEC server j
            % 计算从移动设备i到每个MEC服务器j的J_s(i，j)
            for j = 1: M
                %disp(['MEC server #', num2str(j)])
                h = h_mat(i, j);
                
                E_tmp = sigma * L * log(2) / (omega * h);
                p_L_taud = (power(2, L / (omega * tau_d)) - 1) * sigma / h;
                % calculate p_L
                if E_tmp >= E_min
                    p_L = p_L_taud;
                else
                    % calculate p_E_min (use inline function and fsolve)
                    % 计算p_Emin
                    y = @(x) x * L - omega * log2(1 + h*x/sigma) * E_min;%% 函数句柄  y(x)  代替  “x * L - omega * log2(1 + h*x/sigma) * E_min”
                    % accroding to the function figure, p_L_taud is a positive 
                    % number around 0.2
                    p_E_min = fsolve(y, 0.2, opt);
                    p_L = max(p_L_taud, p_E_min);
                end
                % calculate p_U
                if E_tmp >= E_max
                    p_U = 0;
                else
                    % caculate p_E_max (use inline function and fsolve)
                    % 计算p_Emax
                    y = @(x) x * L - omega * log2(1 + h*x/sigma) * E_max;
                    % accroding to the function figure, p_E_max is a large positive
                    % number around 20
                    p_E_max = fsolve(y, 100, opt);
                    p_U = min(p_tx_max, p_E_max);
                end
                
                if p_L <= p_U
                    % the sub-problem is feasible
                    % P_SE有解
                    %% disp
%                    disp('MEC server execution ($\mathcal{P}_{SE}$) is feasible!')
                    % calculate p_0
                    % 计算p0
                    virtual_battery = B_hat(t, i);
                    y = @(x) virtual_battery * log2(1 + h*x/sigma) + ...
                        h * (V - virtual_battery*x) / log(2) / (sigma + h*x);
                    p_0 = fsolve(y, 0.5, opt);
                    
                    if (p_U < p_0 && B_hat(t, i) < 0) || B_hat(t, i) >= 0
                        p_mat(i, j) = p_U;
                    elseif p_0 < p_L && B_hat(t, i) < 0
                        p_mat(i, j) = p_L;
                    elseif p_0 >= p_L && p_0 <= p_U && B_hat(t, i) < 0
                        p_mat(i, j) = p_0;
                    end
                    % check whether p_mat(i, j) is zero
                    
                    if p_mat(i, j) == 0
                        disp('Something wrong! p is 0!')
                    end
                    
                    % calculate the delay of MEC server execution  计算此时的延迟
                    server_cost_mat(i, j) = L / (omega * log2(1 + h*p_mat(i, j)/sigma));
                    % calculate the energy consumption of MEC server execution计算此时的能耗
                    server_E_mat(i, j) = p_mat(i, j) * server_cost_mat(i, j);
                    % calculate the value of optimization goal 计算此时的J_s， J_s是卸载执行子问题的目标函数值，而非延迟！！
                    J_s(i, j) = (-B_hat(t, i) * p_mat(i, j) + V) * server_cost_mat(i, j);
                    
                    % (we can not set server_exe_cost(t, i) and server_exe_E(t, i) for now)
                else
                    % the sub-problem is not feasible because (i) the limited transmit 
                    % power or (ii) time cosumed out of deadline or (iii) the energy 
                    % consumed out of battery energy level
                    % If it is not feasible, it just means that we cannot choose 
                    % 'I_s(i,j)=1'. It dosen't mean that the task has to be dropped.
                    
                    %% disp
%                    disp('MEC server execution ($\mathcal{P}_{SE}$) is not feasible!')
                    p_mat(i, j) = 0;
                    server_cost_mat(i, j) = 0;
                    server_E_mat(i, j) = 0;
                    % 'I_s(i,j)=1' can never be chosen if MEC server execution goal is inf
                    % 因此indicator(t, 2) = 0不变
                    J_s(i, j) = inf;
                    
                    % (we can not set server_exe_cost(t, i) and server_exe_E(t, i) for now)
                % Similarly, we do not check whether the energy cunsumed is larger than
                % battery energy level because the problem $\mathcal{J}_{CO}$ does
                % not have constraint (8).
                end
            end
            
            %% step 1.2.3: choose the (initial) optimal execution mode
            
            J_d = V * phi;
%             disp(['J_m(i):', num2str(J_m(i))])
%             disp(['J_s(i,:):', num2str(J_s(i, :))])
            %% 为第i个移动设备选取最佳模式
            [~, mode] = min([J_m(i), J_s(i, :), J_d]);
            
            
            if mode == 1
                % mobile execution
                chosen_mode(t, i) = 1;
                final_chosen_cost(t, i) = mobile_exe_cost(t, i);
                final_chosen_E(t, i) = mobile_exe_E(t, i);
            elseif mode == (M+2)
                % drop
                chosen_mode(t, i) = 3;
                final_chosen_cost(t, i) = phi;
                final_chosen_E(t, i) = 0;
            else
                % MEC servre execution
                chosen_mode(t, i) = 2;
                % add i, the chosen j, and their J_s(i, j) value into the pairs
                % 形成最优组合矩阵
                device_server_pairs = [device_server_pairs; [i, mode-1, J_s(i, mode-1)]];
            end
        end
    end
    
    %% step 2: allocate connection right for each mobile device who chooses MEC server execution
    %% 为那些选择卸载执行的任务分配服务器(或另选模式)
    while ~isempty(device_server_pairs)
        for j = 1: M
            % find every pair who chooses the MEC server j  查找选择MEC服务器j的移动端
            device_j_pairs = device_server_pairs(device_server_pairs(:, 2) == j, :);
            % find those devices is
            is = device_j_pairs(:, 1);% 移动端列表
            if isempty(is) % 该服务端没有被任何移动设备选择
                %disp(['For current MEC server #', num2str(j), ', no mobile device choose it!'])
                % go to handle next MEC server
                %  执行下一个服务器的检查
                continue;
            end
            if remained_connects(j) >= length(is)% 大于等于一个移动设备选择时，同时一个服务端不能链接超过四个设备
                % every mobile device who chooses j can be connected,
                % set their final modes as 2 and record the final cost and energy consumption for them
                % 设置当前模式为2，卸载执行
                chosen_mode(t, is) = 2;     % this is not necessary % 当前最优模式仍为卸载执行
                p(t, is) = transp(p_mat(is, j));
                server_exe_cost(t, is) = transp(server_cost_mat(is, j));% 计算了当前卸载执行的cost
                server_exe_E(t, is) = transp(server_E_mat(is, j));%计算当前卸载执行的能量消耗
                chosen_server(t, is) = repmat(j, 1, length(is)); %服务器被占用的数量
                final_chosen_cost(t, is) = server_exe_cost(t, is);% 总cost
                final_chosen_E(t, is) = server_exe_E(t, is);% 总能量
                
                % update remained_connects for j
                % %更新当前的链接到服务端j的数量，减去已经链接的数量，最大不超过4个
                remained_connects(j) = remained_connects(j) - length(is);
                % finally, remove them from global pairs, they are done
                % %清空当前的服务端j的链接列表
                device_server_pairs(device_server_pairs(:, 2) == j, :) = [];
            else%服务端链接数量超过了允许范围
                if remained_connects(j) == 0 % 如果服务端剩余可以被接入数量为0
                    % no mobile device can be connected, remove all mobile devices who chooses j from global pairs
                    % 无法接入任何移动设备，清空当前配对
                    device_server_pairs(device_server_pairs(:, 2) == j, :) = [];
                    % set those mobile devices' J_s(i, j) as inf
                    % (which means no matter what happens, the MEC server j can never be chosen)
                    % J_s 设置为无限大，服务器不可以被选择
                    J_s(is, j) = inf;
                    % choose mode again for those i and insert new potential pairs into global pairs
                    % 再次为这些I选择模式，并将新的潜在对插入到全局对中，让这些移动端选择其他的服务器
                    for idx = 1: numel(is)
                        i = is(idx);
                        [~, mode] = min([J_m(i), J_s(i, :), J_d]);
                        if mode == 1
                            chosen_mode(t, i) = 1;
                            final_chosen_cost(t, i) = mobile_exe_cost(t, i);
                            final_chosen_E(t, i) = mobile_exe_E(t, i);
                        elseif mode == (M+2)
                            chosen_mode(t, i) = 3;
                            final_chosen_cost(t, i) = phi;
                            final_chosen_E(t, i) = 0;
                        else
                            chosen_mode(t, i) = 2;
                            device_server_pairs = [device_server_pairs; [i, mode-1, J_s(i, mode-1)]];
                        end
                    end
                else% 服务端还可以被链接，但是移动端的接入数量大于服务端剩余承受能力
                    % some mobile devices can be connected, set their final modes as 2 and remove them from global pairs
                    % besides, record the final cost and energy consumption for them
                    % for the left i', remove them from global pairs and set J_s(i', j) as inf, 
                    % choose mode again for those i' and insert new potential pairs into global pairs
                    % 无法链接的移动端i，重新插入全局列表，分配其他服务端， 此时这个设备只能重新选择次最优的卸载对象
                    
                    % sort the J_s(is, j) and return those lucky idxs
                    % 可以被接入的移动端：
                    [~, idxs] = sort(device_j_pairs(:, 3));
                    for idx = 1: remained_connects(j)
                        i = idxs(idx);
                        chosen_mode(t, i) = 2;
                        p(t, i) = p_mat(i, j);
                        server_exe_cost(t, i) = server_cost_mat(i, j);
                        server_exe_E(t, i) = server_E_mat(i, j);
                        chosen_server(t, i) = j;
                        final_chosen_cost(t, i) = server_exe_cost(t, i);
                        final_chosen_E(t, i) = server_exe_E(t, i);
                        
                        % remove i from global pairs
                        % 从 device_server_pairs 中删除该键值对并同步一系列共同维护的变量，因为 device_server_pairs 只存放欲图卸载执行的设备
                        device_server_pairs(device_server_pairs(:, 1) == i, :) = [];
                    end
                    
                    % update remained_connects for j
                    % 更新服务端剩余可接入的数量为0
                    remained_connects(j) = 0;
                    
                    % for those unlucky mobile devices i', set J_s(i', j) as inf (i' are in currently device_server_pairs(: , j) now)
                    % 那些不可以被接入的移动端（因为服务端已经满了）
                    residual_is = device_server_pairs(device_server_pairs(:, 2) == j, 1);
                    J_s(residual_is, j) = inf;
                    for idx = 1: numel(residual_is)
                        residual_i = residual_is(idx);
                        %  没有服务器可以选了，只能在另外两种mode中选取
                        [~, mode] = min([J_m(residual_i), J_s(residual_i, :), J_d]);
                        if mode == 1
                            chosen_mode(t, residual_i) = 1;
                            final_chosen_cost(t, residual_i) = mobile_exe_cost(t, residual_i);
                            final_chosen_E(t, residual_i) = mobile_exe_E(t, residual_i);
                        elseif mode == (M+2)
                            chosen_mode(t, residual_i) = 3;
                            final_chosen_cost(t, residual_i) = phi;
                            final_chosen_E(t, residual_i) = 0;
                        else
                            chosen_mode(t, residual_i) = 2;
                            device_server_pairs = [device_server_pairs; [residual_i, mode-1, J_s(residual_i, mode-1)]];
                        end
                    end
                end
            end
        end
    end
    
    %% 这里设置保存中间变量
     % 计算每一个移动设备在当前模式下的execution cost
        mobile_exe_cost=mobile_exe_cost; % N个移动设备的local execution delay
        server_exe_cost=server_exe_cost;     %N个移动设备的remote execution delay
        final_chosen_cost= final_chosen_cost;    % N个移动设备最后选择的模式下的execution cost


    
     % 计算每一个移动设备在当前模式下的能耗
        mobile_exe_E =mobile_exe_E;  %  N个移动设备本地执行的能耗
        server_exe_E = server_exe_E;     %  N个移动设备卸载执行的能耗
        final_chosen_E = final_chosen_E;   % N个移动设备在最后选择的模式下的实际能耗
       
    %% step 3: update the battery energy level and go to the next time slot
    %% 更新电池电量并转到下一个时隙
    B(t + 1, :) = B(t, :) - final_chosen_E(t, :) + e(t, :);
    t = t + 1;
    
end

%% step 4: evaluate the simulation results
average_ratio = zeros(T, 3);
mobile_exe = 0; server_exe = 0; drop = 0;
request_num = 0;
i = 1;          % we simply choose the first mobile device
for t = 1: T
    if final_chosen_cost(t, i) == 0
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
% figure
% plot(1:T, average_ratio(:, 1));
% hold on
% plot(1:T, average_ratio(:, 2));
% hold on
% plot(1:T, average_ratio(:, 3));
% legend('mobile execution', 'MEC server execution', 'drop')
% title('Envolution of average ratio of chosen modes')
% xlabel('time slot')
% ylabel('average  ratio of chosen modes $\frac{1}{T} \sum_{t=0}^{T-1} \{I_m^t, I_s^t, I_d^t\}$ of the i-th mobile device', 'Interpreter','latex')

%ratio =  mean(average_ratio(:,2));
ratio =  average_ratio(:,2);
end





