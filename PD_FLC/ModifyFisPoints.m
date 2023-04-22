function modified = ModifyFisPoints(Member_functions, fis, i, numBitsOnePin)

   epsilon = 0.015;
   %now we must include our membership functions
   for j=1:7
      fis.input(1).mf((j-1) + 1).params(2) = sign((Member_functions( (j-1)*7 + 3,i)-0.5))*(1/(1+Member_functions((j-1)*7 +4,i)*2+Member_functions((j-1)*7 +5,i))) + j*epsilon;
      fis.input(1).mf((j-1) + 1).params(1) = fis.input(1).mf((j-1) +1).params(2) - 1/(1+Member_functions((j-1)*7 +1,i)*2+Member_functions((j-1)*7 +2,i)) + j*epsilon;
      fis.input(1).mf((j-1) + 1).params(3) = fis.input(1).mf((j-1) +1).params(2) + 1/(1+Member_functions((j-1)*7 +6,i)*2+Member_functions((j-1)*7 +7,i)) - j*epsilon;
   end
   for j=1:7
      fis.input(2).mf((j-1) + 1).params(2) = sign((Member_functions(numBitsOnePin + (j-1)*7 + 3,i)-0.5))*(1/(2+Member_functions(numBitsOnePin + (j-1)*7 +4,i)*2+Member_functions(numBitsOnePin + (j-1)*7 +5,i))) + j*epsilon;
      fis.input(2).mf((j-1) + 1).params(1) = fis.input(2).mf((j-1) +1).params(2) - 1/(2+Member_functions(numBitsOnePin + (j-1)*7 +1,i)*2+Member_functions(numBitsOnePin + (j-1)*7 +2,i)) + j*epsilon;
      fis.input(2).mf((j-1) + 1).params(3) = fis.input(2).mf((j-1) +1).params(2) + 1/(2+Member_functions(numBitsOnePin + (j-1)*7 +6,i)*2+Member_functions(numBitsOnePin + (j-1)*7 +7,i)) - j*epsilon;
   end
   for j=1:7
      fis.output(1).mf((j-1) + 1).params(2) = sign((Member_functions(2*numBitsOnePin + (j-1)*7 + 3,i)-0.5))*(1/(2+Member_functions(2*numBitsOnePin + (j-1)*7 +4,i)*2+Member_functions(2*numBitsOnePin + (j-1)*7 +5,i))) + j*epsilon;
      fis.output(1).mf((j-1) + 1).params(1) = fis.output(1).mf((j-1) +1).params(2) - 1/(2+Member_functions(2*numBitsOnePin + (j-1)*7 +1,i)*2+Member_functions(2*numBitsOnePin + (j-1)*7 +2,i)) + j*epsilon;
      fis.output(1).mf((j-1) + 1).params(3) = fis.output(1).mf((j-1) +1).params(2) + 1/(2+Member_functions(2*numBitsOnePin + (j-1)*7 +6,i)*2+Member_functions(2*numBitsOnePin + (j-1)*7 +7,i)) - j*epsilon;
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
           if fis.input(1).mf(j).params(2)>fis.input(1).mf(j+1).params(2)
               temp2 = fis.input(1).mf(j).params(2);
               fis.input(1).mf(j).params(2) = fis.input(1).mf(j+1).params(2);
               fis.input(1).mf(j+1).params(2) = temp2;
               temp3 = fis.input(1).mf(j).params(3);
               fis.input(1).mf(j).params(3) = fis.input(1).mf(j+1).params(3);
               fis.input(1).mf(j+1).params(3) = temp3;
               temp1 = fis.input(1).mf(j).params(1);
               fis.input(1).mf(j).params(1) = fis.input(1).mf(j+1).params(1);
               fis.input(1).mf(j+1).params(1) = temp1;
               %swapped = 1;
               pos = j;
           end
       end
     
       %n = n - 1;
       %if swapped == 0
       %    break;
       %end
   end
   
     if fis.input(1).mf(1).params(1)>-1 && fis.input(1).mf(1).params(2)>-1
           fis.input(1).mf(1).params(1)=-1;
       end
        if fis.input(1).mf(7).params(3)<1 && fis.input(1).mf(7).params(2)<1
           fis.input(1).mf(7).params(3)=1;
       end
   
         pos = n;
   while flag==1
       %swapped = 0;
       if pos == 0
           break;
       end
       bound = pos;
       pos = 0;
       
       for j=1:n-1
           if fis.input(2).mf(j).params(2)>fis.input(2).mf(j+1).params(2)
               temp2 = fis.input(2).mf(j).params(2);
               fis.input(2).mf(j).params(2) = fis.input(2).mf(j+1).params(2);
               fis.input(2).mf(j+1).params(2) = temp2;
               temp3 = fis.input(2).mf(j).params(3);
               fis.input(2).mf(j).params(3) = fis.input(2).mf(j+1).params(3);
               fis.input(2).mf(j+1).params(3) = temp3;
               temp1 = fis.input(2).mf(j).params(1);
               fis.input(2).mf(j).params(1) = fis.input(2).mf(j+1).params(1);
               fis.input(2).mf(j+1).params(1) = temp1;
               %swapped = 1;
               pos = j;
           end
       end
       %n = n - 1;
       %if swapped == 0
       %    break;
       %end
   end
   
     if fis.input(2).mf(1).params(1)>-1 && fis.input(2).mf(1).params(2)>-1 
           fis.input(2).mf(1).params(1)=-1;
       end
        if fis.input(2).mf(7).params(3)<1 && fis.input(2).mf(7).params(2)<1
           fis.input(2).mf(7).params(3)=1;
       end   
   
       pos = n;
   while flag==1
       %swapped = 0;
       if pos == 0
           break;
       end
       bound = pos;
       pos = 0;
       
       for j=1:n-1
           if fis.output(1).mf(j).params(2)>fis.output(1).mf(j+1).params(2)
               temp2 = fis.output(1).mf(j).params(2);
               fis.output(1).mf(j).params(2) = fis.output(1).mf(j+1).params(2);
               fis.output(1).mf(j+1).params(2) = temp2;
               temp3 = fis.output(1).mf(j).params(3);
               fis.output(1).mf(j).params(3) = fis.output(1).mf(j+1).params(3);
               fis.output(1).mf(j+1).params(3) = temp3;
               temp1 = fis.output(1).mf(j).params(1);
               fis.output(1).mf(j).params(1) = fis.output(1).mf(j+1).params(1);
               fis.output(1).mf(j+1).params(1) = temp1;
               %swapped = 1;
               pos = j;
           end
       end
       %n = n - 1;
       %if swapped == 0
       %    break;
       %end
   end
    
      if fis.output(1).mf(1).params(1)<-1 && fis.output(1).mf(1).params(2)<-1 && fis.output(1).mf(1).params(3)<-1
          
          fis.output(1).mf(1).params(3) = -0.95;
          
          
      end
        
      if fis.output(1).mf(7).params(3)>1 &&  fis.output(1).mf(7).params(2)>1 && fis.output(1).mf(1).params(1)>1
          
           fis.output(1).mf(1).params(1) = 1;
           
      end
      
      
     if fis.output(1).mf(1).params(1) > - 1
         fis.output(1).mf(1).params(1) = -1;
     end
     if fis.output(1).mf(7).params(3) <  1
        fis.output(1).mf(7).params(3) = 1;
     end
         
     modified = fis;
      
end