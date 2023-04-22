function reasoning = createReasoning(n,numRules,numBitsRules,fis)

%Although we have added fis - we won't be using it - it was added so that
%we could only make an easier attempt on generating bits for ant

Reasoning_rules = zeros(numRules,numBitsRules);
All_reasoning_rules = zeros(n,numRules,numBitsRules);
antc1 = [];
antc2 = [];
cons = [];

for i = 1 : numRules
   antc1 = [antc1 fis.rule(i).antecedent(1)]; 
end

for i = 1 : numRules
   antc2 = [antc2 fis.rule(i).antecedent(2)]; 
end

for i = 1 : numRules
   cons = [cons fis.rule(i).consequent]; 
end

%cons = pravilaodvojeno.VarName3;
for i=1:49
        if antc1(i) == 1
         Reasoning_rules(i,1) = 0 ; 
         Reasoning_rules(i,2) = 0 ;
         Reasoning_rules(i,3) = 1 ;
        elseif antc1(i) == 2
         Reasoning_rules(i,1) = 0 ; 
         Reasoning_rules(i,2) = 1 ;
         Reasoning_rules(i,3) = 0 ;
        elseif antc1(i) == 3
         Reasoning_rules(i,1) = 0 ; 
         Reasoning_rules(i,2) = 1 ;
         Reasoning_rules(i,3) = 1 ;
        elseif antc1(i) == 4
         Reasoning_rules(i,1) = 1 ; 
         Reasoning_rules(i,2) = 0 ;
         Reasoning_rules(i,3) = 0 ;
        elseif antc1(i) == 5
         Reasoning_rules(i,1) = 1 ; 
         Reasoning_rules(i,2) = 0 ;
         Reasoning_rules(i,3) = 1 ;
        elseif antc1(i) == 6
         Reasoning_rules(i,1) = 1 ; 
         Reasoning_rules(i,2) = 1 ;
         Reasoning_rules(i,3) = 0 ;
        elseif antc1(i) == 7
         Reasoning_rules(i,1) = 1 ; 
         Reasoning_rules(i,2) = 1 ;
         Reasoning_rules(i,3) = 1 ;
        else
        end
        if antc2(i) == 1
         Reasoning_rules(i,4) = 0 ; 
         Reasoning_rules(i,5) = 0 ;
         Reasoning_rules(i,6) = 1 ;
        elseif antc2(i) == 2
         Reasoning_rules(i,4) = 0 ; 
         Reasoning_rules(i,5) = 1 ;
         Reasoning_rules(i,6) = 0 ;
        elseif antc2(i) == 3
         Reasoning_rules(i,4) = 0 ; 
         Reasoning_rules(i,5) = 1 ;
         Reasoning_rules(i,6) = 1 ;
        elseif antc2(i) == 4
         Reasoning_rules(i,4) = 1 ; 
         Reasoning_rules(i,5) = 0 ;
         Reasoning_rules(i,6) = 0 ;
        elseif antc2(i) == 5
         Reasoning_rules(i,4) = 1 ; 
         Reasoning_rules(i,5) = 0 ;
         Reasoning_rules(i,6) = 1 ;
        elseif antc2(i) == 6
         Reasoning_rules(i,4) = 1 ; 
         Reasoning_rules(i,5) = 1 ;
         Reasoning_rules(i,6) = 0 ;
        elseif antc2(i) == 7
         Reasoning_rules(i,4) = 1 ; 
         Reasoning_rules(i,5) = 1 ;
         Reasoning_rules(i,6) = 1 ;
        else
        end
        if cons(i) == 1
         Reasoning_rules(i,7) = 0 ; 
         Reasoning_rules(i,8) = 0 ;
         Reasoning_rules(i,9) = 1 ;
        elseif cons(i) == 2
         Reasoning_rules(i,7) = 0 ; 
         Reasoning_rules(i,8) = 1 ;
         Reasoning_rules(i,9) = 0 ;
        elseif cons(i) == 3
         Reasoning_rules(i,7) = 0 ; 
         Reasoning_rules(i,8) = 1 ;
         Reasoning_rules(i,9) = 1 ;
        elseif cons(i) == 4
         Reasoning_rules(i,7) = 1 ; 
         Reasoning_rules(i,8) = 0 ;
         Reasoning_rules(i,9) = 0 ;
        elseif cons(i) == 5
         Reasoning_rules(i,7) = 1 ; 
         Reasoning_rules(i,8) = 0 ;
         Reasoning_rules(i,9) = 1 ;
        elseif cons(i) == 6
         Reasoning_rules(i,7) = 1 ; 
         Reasoning_rules(i,8) = 1 ;
         Reasoning_rules(i,9) = 0 ;
        elseif cons(i) == 7
         Reasoning_rules(i,7) = 1 ; 
         Reasoning_rules(i,8) = 1 ;
         Reasoning_rules(i,9) = 1 ;
        else
        end
end

    for i=1:n
        for j=1:numRules
            for k=1:numBitsRules
                All_reasoning_rules(i,j,k) = Reasoning_rules(j,k);
            end
        end
    end
    
    %now we will randomly generate our consequences - so that this process
    %can be started from a random point in our search space
    
%     generateCons = rand(n*numRules*numBitsRules/3,1);
%     for i = 1 : length(generateCons)
%         if generateCons(i)<0.5
%             generateCons(i) = 0;
%         else 
%             generateCons(i) = 1;
%         end
%     end
%     
%     for i=1:n
%         for j=1:numRules
%             for k = 2/3*numBitsRules + 1 : numBitsRules
%                 All_reasoning_rules(i,j,k) = generateCons((i-1)*numRules*numBitsRules/3  + (j-1)*numBitsRules/3   + k - 2/3*numBitsRules);
%             end
%         end
%     end
%     for i=1:n
%         for j=1:numRules
%             test = 0;
%             for k = 2/3*numBitsRules + 1 : numBitsRules
%                 test = test + All_reasoning_rules(i,j,k);
%             end
%             if test == 0 
%                 All_reasoning_rules(i,j,numBitsRules) = 1;
%                 All_reasoning_rules(i,j,numBitsRules-1) = 1;
%                 
%             end
%         end
%     end

    reasoning = All_reasoning_rules;
end