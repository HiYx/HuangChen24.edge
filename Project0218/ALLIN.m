%% 计算卸载技术主要包括卸载决策、资源分配和卸载系统实现这三方面。
%%  本项目为：卸载决策方案
% 基于 Lyapunov 优化的动态卸载算法，属于全部卸载思路，以降低时延为目标 

% 为了降低时延，优化目标包含了执行时延和执行故障优化两部分。
% 对这两个方面的优化，不仅能使任务时延最小化，还能保证故障率最低，降低了卸载失败的风险。
% 作者提出采用动态电压频率调整和功率控制的技术分别优化计算执行过程和计算卸载的数据传送。
% 一种基于 Lyapunov 优化的动态卸载（LODCO, low-complexity Lyapunov optimization based dynamic computation offloading）算法。
% LODCO 算法会在每个时隙中进行卸载决定，然后为 UE 分配 CPU 周期（在本地执行）或分配传输功率（卸载到 MEC），
% 结果表明能将运行时间缩短64%。

%  LODCO-Based eps-Greedy Algorithm  对应  LODCO_based_e_Greedy.m （贪心算法） epsilon-Greedy算法算是MBA(Multiarmed Bandit Algorithms)算法中最简单的一种
%  LODCO-Based Genetic Algorithm with Greedy Policy  对应  LODCO_Based_Genetic_Algorithm.m（遗传算法）
% 论文中提出的算法



