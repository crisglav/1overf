function parsave(fname,data)
% Saves variables to a file. To be called inside a parfor loop.
s.(inputname(2)) = data;
save(fname,'-struct','s');
end

