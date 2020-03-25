function dataout = db20(data)
% dataout = db20(data)
% й╣ож20db

dataout = 20*log10(abs(data)/(eps+max(abs(data(:))))+eps);
end