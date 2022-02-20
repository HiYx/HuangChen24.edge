function r = calAchieveRate(h, p, omega, sigma)
%% calculate achievable rate 计算可达到的率 
r = omega*log2(1+h*p/sigma);

end