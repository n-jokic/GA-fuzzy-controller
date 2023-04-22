function retArr = MergeMatrix2(arr1, arr2)

retArr = [];

[d1 d2] = size (arr1);
[d3 d4] = size (arr2);

if d1 == d3
   retArr = zeros(d1, d2 + d4);
   
   for i = 1 : d1
      for j = 1 : d2
          retArr(i,j) = arr1(i,j);
      end
   end
   
   for i =  1 : d3
      for j = 1 : d4
          retArr(i,j + d2) = arr2(i,j);
      end
   end 
else
   disp('Error');
   %stop;
end

end