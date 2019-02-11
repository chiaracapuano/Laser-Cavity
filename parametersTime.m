function [ dt,t,w] = parametersTime( t_max,nn )
dt =2* t_max/(nn); % time spacing in s
t = linspace(-t_max,t_max,nn);
w = linspace(0,1,nn)/dt;

end

