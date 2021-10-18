%%
%
%
%
%

for num=100:num<235
	archivo = ['/home/eliasg/Documents/MATLAB/record/',num2str(num)];
	if exist(cat(2,archivo,'.dat'))
		load([num2str(num),'_QRS_detection.mat']);
		csvwrite('csvlist.dat',series_performance.conf_mat{:,:,1})
		type csvlist.dat


	end
end
