function modified = decodeMembershipFunctions(Member_functions, numBitsOnePin)

   epsilon = 0.015;
   input1 = zeros(7,3);
   input2 = zeros(7,3);
   output1 = zeros(7,3);
   %now we must include our membership functions
   
   
   for j=1:7
      input1((j-1) + 1,2) = sign((Member_functions( (j-1)*7 + 3)-0.5))*(1/(1+Member_functions((j-1)*7 +4)*2+Member_functions((j-1)*7 +5))) + j*epsilon;
      input1((j-1) + 1,1) = input1((j-1) + 1,2) - 1/(1+Member_functions((j-1)*7 +1)*2+Member_functions((j-1)*7 +2)) + j*epsilon;
      input1((j-1) + 1,3) = input1((j-1) + 1,2) + 1/(1+Member_functions((j-1)*7 +6)*2+Member_functions((j-1)*7 +7)) - j*epsilon;
   end
   for j=1:7
      input2((j-1) + 1,2) = sign((Member_functions(numBitsOnePin + (j-1)*7 + 3)-0.5))*(1/(2+Member_functions(49 + (j-1)*7 +4)*2+Member_functions(49 + (j-1)*7 +5))) + j*epsilon;
      input2((j-1) + 1,1) = input2((j-1) + 1,2) - 1/(2+Member_functions(49 + (j-1)*7 +1)*2+Member_functions(49 + (j-1)*7 +2)) + j*epsilon;
      input2((j-1) + 1,3) = input2((j-1) + 1,2) + 1/(2+Member_functions(49 + (j-1)*7 +6)*2+Member_functions(49 + (j-1)*7 +7)) - j*epsilon;
   end
   for j=1:7
      output1((j-1) + 1,2) = sign((Member_functions(2*numBitsOnePin + (j-1)*7 + 3)-0.5))*(1/(2+Member_functions(2*49 + (j-1)*7 +4)*2+Member_functions(2*49 + (j-1)*7 +5))) + j*epsilon;
      output1((j-1) + 1,1) = output1((j-1) + 1,2) - 1/(2+Member_functions(2*49 + (j-1)*7 +1)*2+Member_functions(2*49 + (j-1)*7 +2)) + j*epsilon;
      output1((j-1) + 1,3) = output1((j-1) + 1,2) + 1/(2+Member_functions(2*49 + (j-1)*7 +6)*2+Member_functions(2*49 + (j-1)*7 +7)) - j*epsilon;
   end
   %now we must sort our central points and prolongations towards both
   %sides by the positive of central points 
   
   swapped = 0;
   n = 7;
   flag=1;
   
   pos = n;
   while flag==1
       %swapped = 0;
       if pos == 0
           break;
       end
       bound = pos;
       pos = 0;
       for j=1:bound-1
           if input1(j,2)>input1(j+1,2)
               temp2 = input1(j,2);
               input1(j,2) = input1(j+1,2);
               input1(j+1,2) = temp2;
               temp3 = input1(j,3);
               input1(j,3) = input1(j + 1,3);
               input1(j + 1,3) = temp3;
               temp1 = input1(j,1);
               input1(j,1) = input1(j + 1,1);
               input1(j + 1,1) = temp1;
               %swapped = 1;
               pos = j;
           end
       end
     
       %n = n - 1;
       %if swapped == 0
       %    break;
       %end
   end
   
     if input1(1,1)>-1 && input1(1,2)>-1
           input1(1,1)=-1;
       end
        if input1(7,3)<1 && input1(7,2)<1
           input1(7,3)=1;
       end
   
         pos = n;
   while flag==1
       %swapped = 0;
       if pos == 0
           break;
       end
       bound = pos;
       pos = 0;
       
       for j=1:bound-1
           if input2(j,2)>input2(j+1,2)
               temp2 = input2(j,2);
               input2(j,2) = input2(j+1,2);
               input2(j+1,2) = temp2;
               temp3 = input2(j,3);
               input2(j,3) = input2(j + 1,3);
               input2(j + 1,3) = temp3;
               temp1 = input2(j,1);
               input2(j,1) = input2(j + 1,1);
               input2(j + 1,1) = temp1;
               %swapped = 1;
               pos = j;
           end
       end
       %n = n - 1;
       %if swapped == 0
       %    break;
       %end
   end
   
     if input2(1,1)>-1 && input2(1,2)>-1
           input2(1,1)=-1;
       end
        if input2(7,3)<1 && input2(7,2)<1
           input2(7,3)=1;
       end   
   
       pos = n;
   while flag==1
       %swapped = 0;
       if pos == 0
           break;
       end
       bound = pos;
       pos = 0;
       
       for j=1:bound-1
           if output1(j,2)>output1(j+1,2)
               temp2 = output1(j,2);
               output1(j,2) = output1(j+1,2);
               output1(j+1,2) = temp2;
               temp3 = output1(j,3);
               output1(j,3) = output1(j + 1,3);
               output1(j + 1,3) = temp3;
               temp1 = output1(j,1);
               output1(j,1) = output1(j + 1,1);
               output1(j + 1,1) = temp1;
               %swapped = 1;
               pos = j;
           end
       end
       %n = n - 1;
       %if swapped == 0
       %    break;
       %end
   end
   
     if output1(1,1)<-1 &&  output1(1,2)<-1 && output1(1,3) < -1
         output1(1,3) = -0.95;
     end    
     if output1(7,3) > 1 && output1(7,2) > 1 && output(7,1) > 1
          output1(7,1)= 0.95;
     end    
     
     
     if output1(1, 1) > - 1
          output1(1,1) = -1;
      end
      if output1(7,3) < 1
          output1(7,3) = 1;
          
      end
      
      modified = [{input1}, {input2}, {output1}];
      
      
      
      
      disp('input1:');      
      input1
      disp('input2:');
      input2
      disp('output1:');
      output1



end