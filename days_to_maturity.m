function dtm = days_to_maturity(nowDate, matureDate, tradingDays)

% Function to calculate days to maturity of an option
% Input:
% nowDate: current date (T*1 array)
% matureDate: mature date (T*1 array)
% tradingDays: an array of trading days (be sure to cover the earliest 
%              nowDate and the latest matureDate)

dtm = NaN(size(nowDate,1),1);

for i=1:size(nowDate,1)
    dtm(i) = length(tradingDays(tradingDays<=matureDate(i) & tradingDays>nowDate(i)));
end

end