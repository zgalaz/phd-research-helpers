function [split_data, intervals] = ...
    split_data(data, num_intervals, split_function)

% [split_data, intervals] = split_data(data, num_intervals, split_function)
%
% This function splits an input data into selected number of intervals. It
% returns a matrix of logical values (rows: input values of 'data' vector;
% cols: intervals) denoting the presence of data values in each interval.
% 
% data              - input data (column vector)
% num_intervals     - number of intervals to split the data into
% split_function    - function to split the data (e.g. @linspace)
%
% intervals         - intervals determined by the chosen function
% split_data        - matrix of logical values (rows: input values of 
%                     data vector; cols: intervals) denoting the presence
%                     of data values in each interval
%
%
%
% --
% ing. Zoltán Galáž
% xgalaz00@stud.feec.vutbr.cz       
% 
% Department of Telecommunications
% Faculty of Electrical Engineering and Communication
% Brno University of Technology

%% Paths and variables
if ((nargin < 3) || (isempty(split_function)))
    split_function = @linspace;
end
if ((nargin < 2) || (isempty(num_intervals)))
    if (length(data) >= 2)
        num_intervals = 2;
    else
        num_intervals = 1;
    end
end

%% Split the data into num_interval using chosen split technique
%
%  Split technique: @linspace, @logspace
%  num_intervals:    whole number, e.g. 4 intervals for data in range 
%                    0 - 100 => (<0, 25) <25, 50) <50, 75) <75, 100))
%                    splitted linearly using @linspace (by default)
split_data = false(length(data), num_intervals);
intervals  = split_function(min(data), max(data), num_intervals + 1);

for int = 1:num_intervals
    for dat = 1:length(data)
        if ((data(dat) >= intervals(int)) && ...
            (data(dat) <  intervals(int + 1)))
            split_data(dat, int) = true;
        end
    end
end