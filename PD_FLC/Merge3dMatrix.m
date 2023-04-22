function retArr = Merge3dMatrix(arr1, arr2)

retArr = [];

[d1 d2 d5] = size (arr1);
[d3 d4 d6] = size (arr2);

if d5 == d6 && d2 == d4
   retArr = zeros(d1+d3, d2 , d5);
   
   for i = 1 : d1
      for j = 1 : d2
          for k = 1 : d5
             retArr(i,j,k) = arr1(i,j,k);
          end
      end
   end
   
   for i =  1 : d3
      for j = 1 : d4
          for k = 1 : d5
             retArr(d1 + i,j,k) = arr2(i,j,k);
          end
      end
   end 
else
   disp('Error');
   stop;
end

end