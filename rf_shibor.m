function r = rf_shibor(yearToMaturity,rf_1D,rf_1W,rf_2W,rf_1M,rf_3M,rf_6M,rf_9M,rf_1Y)

% Function to calculate risk free rate based on SHIBOR
% Method: linear interpolation

r = NaN(size(yearToMaturity,1),1);

c_1 = 1 / 365;
c_2 = 7 / 365;
c_3 = 14 / 365;
c_4 = 1 / 12;
c_5 = 3 / 12;
c_6 = 6 / 12;
c_7 = 9 / 12;
c_8 = 1;

for i=1:size(yearToMaturity,1)

    if (yearToMaturity(i) < c_1)
        r(i) = NaN;
    elseif (c_1 <= yearToMaturity(i) && yearToMaturity(i) < c_2)
        r = rf_1D + (yearToMaturity(i) - c_1) / (c_2 - c_1) * (rf_1W(i) - rf_1D(i));
    elseif (c_2 <= yearToMaturity(i) && yearToMaturity(i) < c_3)
        r = rf_1W + (yearToMaturity(i) - c_2) / (c_3 - c_2) * (rf_2W(i) - rf_1W(i));
    elseif (c_3 <= yearToMaturity(i) && yearToMaturity(i) < c_4)
        r = rf_2W + (yearToMaturity(i) - c_3) / (c_4 - c_3) * (rf_1M(i) - rf_2W(i));
    elseif (c_4 <= yearToMaturity(i) && yearToMaturity(i) < c_5)
        r = rf_1M + (yearToMaturity(i) - c_4) / (c_5 - c_4) * (rf_3M(i) - rf_1M(i));
    elseif (c_5 <= yearToMaturity(i) && yearToMaturity(i) < c_6)
        r = rf_3M + (yearToMaturity(i) - c_5) / (c_6 - c_5) * (rf_6M(i) - rf_3M(i));
    elseif (c_6 <= yearToMaturity(i) && yearToMaturity(i) < c_7)
        r = rf_6M + (yearToMaturity(i) - c_6) / (c_7 - c_6) * (rf_9M(i) - rf_6M(i));
    elseif (c_7 <= yearToMaturity(i) && yearToMaturity(i) <= c_8)
        r = rf_9M + (yearToMaturity(i) - c_7) / (c_8 - c_7) * (rf_1Y(i) - rf_9M(i));
    elseif (c_8 < yearToMaturity(i))
        r(i) = NaN;
        disp('Invalid yearToMaturity');
    
    end
end

end
