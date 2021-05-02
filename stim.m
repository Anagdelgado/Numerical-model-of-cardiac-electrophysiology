function [I_stim] = stim(n,dt)
if n*dt <= 5
    I_stim = 8;
else
    I_stim=0;
end