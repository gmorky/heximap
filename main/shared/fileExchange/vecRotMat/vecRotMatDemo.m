function vecRotMatDemo()
vec1 = rand(100000,3);
vec1 = bsxfun(@rdivide,vec1,sqrt(sum(vec1.^2,2)));
vec2 = rand(100000,3);
vec2 = bsxfun(@rdivide,vec2,sqrt(sum(vec2.^2,2)));
R = vecRotMat(vec1,vec2);
%Test to ensure that R*vec1' = vec2'
sol = [squeeze(sum(bsxfun(@times,R(1,:,:),permute(vec1,[3 2 1])),2)), ...
       squeeze(sum(bsxfun(@times,R(2,:,:),permute(vec1,[3 2 1])),2)), ...
       squeeze(sum(bsxfun(@times,R(3,:,:),permute(vec1,[3 2 1])),2))];
%Compute and Report Maximum Error in Rotation Matrices   
RMS = sqrt(sum((sol-vec2).^2,2)); 
disp(['MAX RMS ERROR:' num2str(max(RMS))]);
if isequal(single(sol),single(vec2))
   disp('TEST PASSED');
else
   error('TEST FAILED');
end

end