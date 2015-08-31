function [idx_signal, dur_signal] = cum_time_for_segment(signal_vec, time_vec)
% find the location and duration of a signal satifying certain constaints
% signal_vec - 1D

% n = length(signal_vec);
% 
% h = 1;
idx_signal = [];
dur_signal = [];

idx_start = (signal_vec(1:end-1) == 0) & (signal_vec(2:end) == 1);
idx_end = (signal_vec(1:end-1) == 1) & (signal_vec(2:end) == 0);

if sum(idx_start) == sum(idx_end)
    idx_signal = idx_start;
    dur_signal = time_vec(idx_end) - time_vec(idx_start);
elseif sum(idx_start) < sum(idx_end)
    idx_signal = idx_start;
    dur_signal = time_vec(idx_end(2:end)) - time_vec(idx_start);
elseif sum(idx_start) > sum(idx_end)
    idx_signal = idx_start(1:end-1);
    dur_signal = time_vec(idx_end) - time_vec(idx_start(1:end-1));
end
