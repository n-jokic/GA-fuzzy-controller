function retArr = MergeMatrix(arr1, arr2)

retArr = [];

[d1 d2] = size (arr1);
[d3 d4] = size (arr2);

if d2 == d4
   retArr = zeros(d1+d3, d2);
   
   for i = 1 : d1
      for j = 1 : d2
          retArr(i,j) = arr1(i,j);
      end
   end
   
   for i =  1 : d3
      for j = 1 : d4
          retArr(d1 + i,j) = arr2(i,j);
      end
   end 
else
   disp('Error');
   %stop;
end

end