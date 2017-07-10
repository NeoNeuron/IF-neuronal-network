% Estimate the proper intensity for Poisson feedforward driving
% Parameters
g_m = 5e-2;
f_e = 5e-3;
e_e = 14/3;
tau = 2;
n = 35;

g = @(t) n * f_e * exp(- t / tau);

tmax = 100;
dt = 1e-3;

v = zeros(tmax / dt + 1, 1);

for i = 2:length(v)
	v(i) = v(i - 1) + dt*(-(g_m + g(i * dt))*v(i-1) + e_e*g(i*dt));
end

plot(0:dt:tmax, v)

find(v(:) == max(v))