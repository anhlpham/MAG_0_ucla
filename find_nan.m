function [arr_out] = find_nan(arr_in)
inds_nan = find(isnan(arr_in)==1);
arr_out = arr_in;
arr_out(inds_nan)=0;
