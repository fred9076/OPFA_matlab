function dataout = db20(data)

dataout = 20*log10(abs(data)/(eps+max(abs(data(:))))+eps);
end