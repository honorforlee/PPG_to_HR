function [t_ s_] = time_div(t,s,dt,m)
range = (1 : (10/dt)) * dt;

%   - Divide timeline -
for k = 0 : length(s) / length(range) - 1
    
    t_div (k+1,:) =  t(  k*(length(range)) + 1 : (k+1)*(length(range)) ) ;
    s_div(k+1,:)= s(  k*(length(range)) + 1 : (k+1)*(length(range)) ) ;
    
end

t_ = t_div(m,:);    s_ = s_div(m,:);    