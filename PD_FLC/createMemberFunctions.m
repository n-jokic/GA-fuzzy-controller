function  Member_functions= createMemberFunctions(n, numBitsOnePin)

%Now we will create a population of defininitions for our member functions
%for both of the inputs and our desired output

All_mf = rand(3*n*numBitsOnePin,1);
for i=1:length(All_mf)
   if All_mf(i)<0.5
       All_mf(i) = 0;
   else
      All_mf(i) = 1; 
   end
end

Member_functions = zeros(3*numBitsOnePin,n);

for j=1:3*numBitsOnePin
   for i=1:n
       Member_functions(j,i)=All_mf((j-1)*n+i);
   end
end

end