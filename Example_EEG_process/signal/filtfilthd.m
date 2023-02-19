function filtered_data=filtfilthd(Hd,data)
%Takes the output of the filter design, and applies it to data in
%forward and backward direction to obtain a zero phase filter. 
%Always applies filter to the first dimension of the data matrix.
%filter design is generated using makefilter
%Input: Hd - output of makefilter
%       data - data matrix with time as the first dimension
%Output: filtered_data = data after filter is applied.    
filtered_data = filtfilt(Hd.sosMatrix,Hd.ScaleValues,data);