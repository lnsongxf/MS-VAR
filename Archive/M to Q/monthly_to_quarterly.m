function [Q,M] = monthly_to_quarterly(data,method)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert monthly to quarterly data
%
% Written by Danilo Cascaldi-Garcia
%
% Fisrt observation of data is the first month of the quarter
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Method: 1, if mean
%         2, if last observation

control =0;
while control ==0
    T = size(data,1);
    test_size = mod(T,3);
    if test_size ==0
        control =1;
    else
        data = [data;NaN(1,size(data,2))];
        continue        
    end
end
Q = NaN(T/3,size(data,2));
for jj=1:T/3
    temp = data((jj-1)*3+1:(jj-1)*3+3,:);
    if method ==1
        Q(jj,:) = mean(temp(~isnan(temp(:,1)),:));
    elseif method ==2
        Q(jj,:) = temp(end,:);
    end
end
M = data;