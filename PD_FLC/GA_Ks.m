clear all
close all
clc

%This is a piece of code which we'll use to implement our final
%optimization - the one that includes optimization of Ks, Reasoning rules,
%and definitions of membership functions

nP = 30; %population size
rC = 0.7; %crossover probability
rM = 0.08; % mutation probability
maxIter = 50; %criterion for convergence
maxKp = 50; %the following three will be the constants which we'll use to scale Ks to their respective ranges
maxK = 255;
maxKd = 1.5;
numRules = 49; 
numBitsRules = 9;
numBitsOnePin = 7*7;


K = 10
Kp = 10
Kd = 10

%Now we list parameters for our simulation
t_end = 10;
t_start = 0;
Num = [2];
Den = [1 12 24];
fis_file_name = 'pravila';
fis = readfis(fis_file_name);
sim_file_name = 'PD_FLC_simulink';
open_system(sim_file_name);



%we'll use back up fis file just in case we need the original one for samo
%modification

fis = readfis(fis_file_name);
writefis(fis, 'rezervni fis');
fis = readfis ('rezervni fis');




K_correct = 1 + 0.5 + 0.5*0.5 + 0.5*0.5*0.5 + 0.5*0.5*0.5*0.5 + 0.5*0.5*0.5*0.5*0.5 + 0.5*0.5*0.5*0.5*0.5*0.5;
%Coefficient which we'll use to normalize our Ks in range from 0 to 255

numBitsCoef = 7;

newChildren = round(rC * nP);

%Now we will generate a starting number of possible nP + newChildren population members from
%which we'll extract nP members which will represent our base population
initialSizeKpKdK = rand(3*(newChildren+nP)*numBitsCoef,1);
KpInitial = zeros(nP + newChildren, numBitsCoef);
KInitial = zeros(nP + newChildren, numBitsCoef);
KdInitial = zeros(nP + newChildren, numBitsCoef);

for i=1:length(initialSizeKpKdK)
    if initialSizeKpKdK(i)<0.5
        initialSizeKpKdK(i)=0;
    else
        initialSizeKpKdK(i)=1;
    end
end

for i=1:nP + newChildren
   for j=1:numBitsCoef
       KpInitial(i,j)=initialSizeKpKdK((i-1)*numBitsCoef+j);
   end
end

for i=1:nP + newChildren
   for j=1:numBitsCoef
       KdInitial(i,j)=initialSizeKpKdK((nP + newChildren)*numBitsCoef + (i-1)*numBitsCoef + j);
   end
end

for i=1:nP + newChildren
   for j=1:numBitsCoef
       KInitial(i,j)=initialSizeKpKdK(2*(nP + newChildren)*numBitsCoef+(i-1)*numBitsCoef+j);
   end
end

%Now we will create a population of the reasoning rules for our fuzzy
%controler - we will have a hundred sets of rules for them
All_reasoning_rules_initial = createReasoning(nP+newChildren,numRules,numBitsRules,fis);

%Now we will create a population of parameters which describe our
%membership functions

All_mf_initial = createMemberFunctions(nP + newChildren,numBitsOnePin); 

%Now we'll do the extraction part
KpPopulation = zeros(nP, numBitsCoef);
KPopulation = zeros(nP, numBitsCoef);
KdPopulation = zeros(nP, numBitsCoef);
FitnessPopulation = zeros(nP,1);
ReasoningRulesPopulation = zeros(nP, numRules, numBitsRules);
MembershipFunctionsPopulation = zeros(3*numBitsOnePin, nP);

fitness_vector = zeros(nP + newChildren, 1);

for i=1: nP + newChildren
 
    Kp = 0; Kd = 0; K = 0;
    for j=1:numBitsCoef
        Kp = Kp + KpInitial(i,j)*(0.5^(j-1));
        K = K + KInitial(i,j)*(0.5^(j-1));
        Kd = Kd + KdInitial(i,j)*(0.5^(j-1));
    end
    Kp = Kp*maxKp/K_correct;
    K = K*maxK/K_correct;
    Kd = Kd*maxKd/K_correct;
    
    %we must not forget to include our consequences into consideration
    %in case a consequence is set to zero, we'll handle it by automatically
    %trying to set it to a value between 1 and 7, and if, in the process,
    %we accidentally get 8, we'll set it to 7 by default
    
       for j = 1 : numRules
           
       fis.rule(j).consequent = 0;    
       for k = 2*numBitsRules/3 + 1 : numBitsRules
             fis.rule(j).consequent = fis.rule(j).consequent + All_reasoning_rules_initial(i,j,k)*(2^(numBitsRules - k));
       end
       
       if fis.rule(j).consequent == 0
           fis.rule(j).consequent = 1;
           All_reasoning_rules_initial(i,j,numBitsRules) = 1;
       end
       if fis.rule(j).consequent == 8
           fis.rule(j).consequent = 7;
           for k = 2*numBitsRules/3 + 1 : numBitsRules
               All_reasoning_rules_initial(i,j,k) = 1;
           end
       end
    end
    
    fis = ModifyFisPoints(All_mf_initial, fis, i, numBitsOnePin);

    fitness_vector(i) = returnFitness(Kp,K,Kd,t_end, t_start, Num, Den, fis,sim_file_name);
end

%Now we will sort Ks, reasoning rules and function points going by their fitness value (note that for two
%fitness values the better one will be the one that's smaller - IAE
%criterion

len = length(fitness_vector);
tempRules = zeros(numRules, numBitsRules);
tempFunctionPoints = zeros(3*numBitsOnePin);
for i = 1 : len-1
   swapped = 0;
   
   for j = 1 : len - 1
       % now we will swap both the fitness values and their corresponding Ks 
        if fitness_vector(j)>fitness_vector(j+1)
            tempFV = fitness_vector(j);
            fitness_vector(j) = fitness_vector(j + 1);
            fitness_vector(j + 1) = tempFV;
            
            for k = 1 : numBitsCoef
                tempKpkthBit = KpInitial(j,k);
                KpInitial(j,k) = KpInitial(j+1,k);
                KpInitial(j+1,k) = tempKpkthBit;
                
                tempKkthBit = KInitial(j,k);
                KInitial(j,k) = KInitial(j+1,k);
                KInitial(j+1,k) = tempKkthBit;
                
                tempKdkthBit = KdInitial(j,k);
                KdInitial(j,k) = KdInitial(j+1,k);
                KdInitial(j+1,k) = tempKdkthBit;
            end
            
            %We must also sort our reasoning rules in the appropriate order 
            
            for m = 1 : numRules
               for k = 1 : numBitsRules
                   tempRules(m,k) = All_reasoning_rules_initial(j,m,k);
               end
            end
            
             for m = 1 : numRules
               for k = 1 : numBitsRules
                   All_reasoning_rules_initial(j,m,k) = All_reasoning_rules_initial(j + 1,m,k);
               end
             end
             
             for m = 1 : numRules
               for k = 1 : numBitsRules
                   All_reasoning_rules_initial(j + 1,m,k) = tempRules(m,k);
               end
             end
             
             %and now we must also sort our membership function points
             for m = 1 : 3*numBitsOnePin
                 tempFunctionPoints(m) = All_mf_initial(m,j);
             end
             
             for m = 1 : 3*numBitsOnePin
                 All_mf_initial(m,j) =  All_mf_initial(m,j + 1);
             end
             
             for m = 1 : 3*numBitsOnePin
                 All_mf_initial(m,j + 1) =  tempFunctionPoints(m);
             end
            
            swapped = 1;
        end
   end
   if swapped == 0
       break;
   end
end

%now it's time we selected our initial base population as the best members
%of this raw set

KpPopulation = KpInitial(1:nP,:);
KPopulation = KInitial(1:nP,:);
KdPopulation = KdInitial(1:nP,:);
FitnessPopulation = fitness_vector(1:nP);
ReasoningRulesPopulation = All_reasoning_rules_initial(1:nP, :, :);
MembershipFunctionsPopulation = All_mf_initial(:,1:nP);

%now we'll create maxIter generations by creating newChildren number of
%childen - for the selection part we will be using roulette selection
%principle; we'll be using crossover principles which bases on breaking
%point; for mutation we'll simulate mutations via new vector of mutations
%and we'randomly mutate a bit in our child

%first, we'll need supporting vectors for our children, their respective
%fintess values and vector of mutations
KpChildren = zeros(newChildren,numBitsCoef);
KChildren = zeros(newChildren,numBitsCoef);
KdChildren = zeros(newChildren,numBitsCoef);
FitnessChildren = zeros(newChildren,1);
ReasoningRulesChildren = zeros(newChildren, numRules, numBitsRules);
MembershipFunctionsChildren = zeros(3*numBitsOnePin, newChildren);

for i = 1 : newChildren
    for j = 1 : numRules
        for k = 1 : 2/3 * numBitsRules
            ReasoningRulesChildren(i,j,k)=All_reasoning_rules_initial(i,j,k);
        end
    end
end

%now let's start iterating through our generations
iter = 0;
num_bits_cons = numBitsRules/3;
num_of_bits_functions = 3 * numBitsOnePin;

childFunctions = zeros(num_of_bits_functions,1);

while iter < maxIter 
    array_of_mutations = rand(newChildren,1);
    for i = 1 : length(array_of_mutations)
        if array_of_mutations(i)<rM
            array_of_mutations(i)=1;
        else
            array_of_mutations(i)=0;
        end
    end
    
    for i = 1 : newChildren
       %first, we must select our parents
        parent1 = floor(rand(1,1)*nP) + 1;
        parent2 = floor(rand(1,1)*nP) + 1;
   
        %now we will use crossover
        childKp = zeros(numBitsCoef,1); 
        childKd = zeros(numBitsCoef,1); 
        childK = zeros(numBitsCoef,1); 

        childsetOfCons = zeros(numRules,numBitsRules/3);
        
        boundary_bit = floor(rand(1,1)*numBitsCoef) + 1;
        boundary_bit_cons = floor(rand(1,1)*num_bits_cons) + 1;
        boundary_bit_functions = floor(rand(1,1)*num_of_bits_functions) + 1;
        for j = 1 : numBitsCoef
            if j < boundary_bit
              childKp(j) = KpPopulation(parent1,j); 
              childKd(j) = KdPopulation(parent1,j);
              childK(j) = KPopulation(parent1,j); 
            else
              childKp(j) = KpPopulation(parent2,j); 
              childKd(j) = KdPopulation(parent2,j);
              childK(j) = KPopulation(parent2,j); 
            end
        end
   
   %Now we will apply mutation operator
        if array_of_mutations(i) == 1
            mutated_bit = floor(rand(1,1)*numBitsCoef) + 1;
            childK(mutated_bit) = mod (childK(mutated_bit)+1,2);
            childKp(mutated_bit) = mod (childKp(mutated_bit)+1,2);
            childKd(mutated_bit) = mod (childKd(mutated_bit)+1,2);
        end 
        
        Kp=0; K=0;Kd=0;
        for j = 1 : numBitsCoef
            KpChildren(i,j) = childKp(j);   
            KChildren(i,j) = childK(j);
            KdChildren(i,j) = childKd(j); 
            
            Kp = Kp + childKp(j)*(0.5^(j-1));
            K = K + childK(j)*(0.5^(j-1));
            Kd = Kd + childKd(j)*(0.5^(j-1));
        end
        
        Kp = Kp*maxKp/K_correct;
        K = K*maxK/K_correct;
        Kd = Kd*maxKd/K_correct;
        
        for k = 1 : numRules
           for j = 1 : num_bits_cons 
              if j<boundary_bit_cons
                childsetOfCons(k,j)=ReasoningRulesPopulation(parent1,k,j+2*numBitsRules/3);
              else
               childsetOfCons(k,j)=ReasoningRulesPopulation(parent2,k,j+2*numBitsRules/3);
              end
           end
        end
        if array_of_mutations(i) == 1
            mutated_bit = floor(rand(1,1)*numBitsRules/3) + 1;
            for k = 1 : numRules
                childsetOfCons(k,mutated_bit) = mod(childsetOfCons(k,mutated_bit)+1,2);
            end
        end
        
        %Now we must secure ourselves that our childsetOfCons is different
        %from zzero
        
        for j=1:numRules
            test = 0;
            for k = 2/3*numBitsRules + 1 : numBitsRules
                test = test + childsetOfCons(j,k - 2/3*numBitsRules);
            end
            if test == 0 
                spareCons = rand(2*numBitsRules/3,1);
                    for k = 1 : length(spareCons)
                        if spareCons(k) < 0.5
                            spareCons(k) = 0;
                        else
                            spareCons(k) = 1;
                        end
                        test = test + spareCons(k);
                    end
                
                        if test == 0 
                            childsetOfCons(j,1) = 1;
                        else
                        for k = 1 : length(spareCons)
                            childsetOfCons(j,k) = spareCons(k);
                        end
                    end
                end
        end

        
        for k = 1 : numRules
           for j = 1 : num_bits_cons 
                ReasoningRulesChildren(i,k,2/3*numBitsRules + j) = childsetOfCons(k,j);
           end
        end
        
        %And now we must apply this new consequences to our rules
        
        for j = 1 : numRules
           
            fis.rule(j).consequent = 0;    
            for k = 2*numBitsRules/3 + 1 : numBitsRules
                fis.rule(j).consequent = fis.rule(j).consequent + ReasoningRulesChildren(i,j,k)*(2^(numBitsRules - k));
            end
       
            if fis.rule(j).consequent == 0
                fis.rule(j).consequent = 2;
                ReasoningRulesChildren(i,j,numBitsRules-1)=1;
            end
            if fis.rule(j).consequent == 8
                fis.rule(j).consequent = 7;
                for k = 2*numBitsRules/3 + 1 : numBitsRules
                    ReasoningRulesChildren(i,j,k) = 1;
                end
            end
        end
        
        %now we must permorfm these operators on our membership functions'
        %points and apply them to our fis file
        
        %first we'll apply crossover
        for m = 1 : num_of_bits_functions
            if i<boundary_bit_functions
                childFunctions(m)=MembershipFunctionsPopulation(m,parent1);
            else
                childFunctions(m)=MembershipFunctionsPopulation(m,parent2);
            end
        end
        
        %and now its time to apply possible mutation
        
       if array_of_mutations(i) == 0
            mutated_bit_functions = floor(rand(1,1)*num_of_bits_functions) + 1;
            childFunctions(mutated_bit_functions) = mod (childFunctions(mutated_bit_functions)+1,2);
       end
       
       %now it's time to insert childFunctions bits into our
       %MembershipFunctionsChildren and to apply it to our fis
       
       for m = 1 : num_of_bits_functions
           MembershipFunctionsChildren(m,i) = childFunctions(m);
       end
       
       fis = ModifyFisPoints(MembershipFunctionsChildren, fis, i, numBitsOnePin);
       
       
       
        %And now it's time we estimated our new child
        FitnessChildren(i) = returnFitness(Kp,K,Kd,t_end, t_start, Num, Den,fis,sim_file_name);
        
    end
    
    %Now, we'll compare both our parents and pur children and we'll choose
    %nP fittest 
     
    KpInitial = MergeMatrix(KpPopulation, KpChildren);
    KInitial = MergeMatrix(KPopulation, KChildren);
    KdInitial = MergeMatrix(KPopulation, KChildren);
    All_reasoning_rules_initial = Merge3dMatrix(ReasoningRulesPopulation, ReasoningRulesChildren);
    All_mf_initial = MergeMatrix2(MembershipFunctionsPopulation, MembershipFunctionsChildren);
    fitness_vector = MergeMatrix(FitnessPopulation, FitnessChildren);
    
    %Now, it's choosing time
    
    %First comes the sorting
    len = length(fitness_vector);
    for i = 1 : len-1
        swapped = 0;
   
        for j = 1 : len - 1
            % now we will swap both the fitness values and their corresponding Ks 
            if fitness_vector(j)>fitness_vector(j+1)
                tempFV = fitness_vector(j);
                fitness_vector(j) = fitness_vector(j + 1);
                fitness_vector(j + 1) = tempFV;
            
                for k = 1 : numBitsCoef
                    tempKpkthBit = KpInitial(j,k);
                    KpInitial(j,k) = KpInitial(j+1,k);
                    KpInitial(j+1,k) = tempKpkthBit;
                
                    tempKkthBit = KInitial(j,k);
                    KInitial(j,k) = KInitial(j+1,k);
                    KInitial(j+1,k) = tempKkthBit;
                
                    tempKdkthBit = KdInitial(j,k);
                    KdInitial(j,k) = KdInitial(j+1,k);
                    KdInitial(j+1,k) = tempKdkthBit;
                end
                
                %We must also sort our reasoning rules in the appropriate order 
            
                for m = 1 : numRules
                    for k = 1 : numBitsRules
                        tempRules(m,k) = All_reasoning_rules_initial(j,m,k);
                    end
                end
            
                for m = 1 : numRules
                    for k = 1 : numBitsRules
                        All_reasoning_rules_initial(j,m,k) = All_reasoning_rules_initial(j + 1,m,k);
                    end
                end
             
                for m = 1 : numRules
                    for k = 1 : numBitsRules
                        All_reasoning_rules_initial(j + 1,m,k) = tempRules(m,k);
                    end
                end
                
                %and now we must also sort our membership function points
                for m = 1 : 3*numBitsOnePin
                    tempFunctionPoints(m) = All_mf_initial(m,j);
                end
             
                for m = 1 : 3*numBitsOnePin
                     All_mf_initial(m,j) =  All_mf_initial(m,j + 1);
                end
             
                for m = 1 : 3*numBitsOnePin
                     All_mf_initial(m,j + 1) =  tempFunctionPoints(m);
                end
            
                swapped = 1;
            end
        end
        if swapped == 0
            break;
        end
    end
    
    %and now the choosing
    KpPopulation = KpInitial(1:nP,:);
    KPopulation = KInitial(1:nP,:);
    KdPopulation = KdInitial(1:nP,:);
    FitnessPopulation = fitness_vector(1:nP);
    ReasoningRulesPopulation = All_reasoning_rules_initial(1:nP, :, :);
    MembershipFunctionsPopulation = All_mf_initial(:,1:nP);
    
    iter = iter + 1;
end

%Now, when we have achieved our convergence criterion, we'll choose the
%best solution from our final population


bestSolutionKp = KpPopulation(FitnessPopulation == min(FitnessPopulation),:); 
bestSolutionK = KPopulation(FitnessPopulation == min(FitnessPopulation),:);
bestSolutionKd = KdPopulation(FitnessPopulation == min(FitnessPopulation),:); 
bestSolutionReasoningRules = ReasoningRulesPopulation((FitnessPopulation == min(FitnessPopulation)),:,:);
bestSolutionMembershipFunctions = MembershipFunctionsPopulation(:,FitnessPopulation == min(FitnessPopulation));

bestKp = 0;
bestK = 0;
bestKd = 0;
for j = 1 : numBitsCoef
        bestKp = bestKp + bestSolutionKp(j)*(0.5^(j-1));
        bestK = bestK + bestSolutionK(j)*(0.5^(j-1));
        bestKd = bestKd + bestSolutionKd(j)*(0.5^(j-1));
end
bestKp = bestKp / K_correct * maxKp
bestK = bestK / K_correct * maxK
bestKd = bestKd / K_correct * maxKd

%now, the only thing left for us is to find interpretation of bestSolutionReasoningRules

bestSolutionAntc1 = zeros(numRules,1);
bestSolutionAntc2 = zeros(numRules,1);
bestSolutionCons = zeros(numRules,1);

for i = 1 : numRules
    for k = 1 : 1/3*numBitsRules
        bestSolutionAntc1(i) = bestSolutionAntc1(i) + bestSolutionReasoningRules(1,i,k)*(2^(1/3*numBitsRules-k));
    end
end

for i = 1 : numRules
    for k = 1/3*numBitsRules + 1 : 2/3*numBitsRules
        bestSolutionAntc2(i) = bestSolutionAntc2(i) + bestSolutionReasoningRules(1,i,k)*(2^(2/3*numBitsRules-k));
    end
end

for i = 1 : numRules
    for k = 2/3*numBitsRules + 1 : numBitsRules
        bestSolutionCons(i) = bestSolutionCons(i) + bestSolutionReasoningRules(1,i,k)*(2^(numBitsRules-k));
    end
    if bestSolutionCons(i)==0
        bestSolutionCons(i) = 4; %exception value
    end;
end

bestReasoning = zeros(numRules,3);

for i = 1 : numRules
    bestReasoning(i,1) = bestSolutionAntc1(i);
end
for i = 1 : numRules
    bestReasoning(i,2) = bestSolutionAntc2(i);
end
for i = 1 : numRules
    bestReasoning(i,3) = bestSolutionCons(i);
end
disp('Best reasoning is:');
bestReasoning
disp('Best membership function points are:');

bestSolutionMembershipFunctions = decodeMembershipFunctions(bestSolutionMembershipFunctions,numBitsOnePin); %Povratna vrednost treba da bude matrica [...]

K = bestK;
Kd = bestKd;
Kp = bestKp;

for j = 1 : numRules
          
       
       fis.rule(j).consequent = bestReasoning(2*49+j);  
       fis.rule(j).antecedent(1) = bestReasoning(0*49+j);  
       fis.rule(j).antecedent(2) = bestReasoning(1*49+j);  

end

for j = 1 : 7
    
    fis.input(1).mf(j).params(1) = bestSolutionMembershipFunctions{1}(j, 1);
    fis.input(1).mf(j).params(2) = bestSolutionMembershipFunctions{1}(j, 2);
    fis.input(1).mf(j).params(3) = bestSolutionMembershipFunctions{1}(j, 3);
    
    
    fis.input(2).mf(j).params(1) = bestSolutionMembershipFunctions{2}(j, 1);
    fis.input(2).mf(j).params(2) = bestSolutionMembershipFunctions{2}(j, 2);
    fis.input(2).mf(j).params(3) = bestSolutionMembershipFunctions{2}(j, 3);
    
    fis.output(1).mf(j).params(1) = bestSolutionMembershipFunctions{3}(j, 1);
    fis.output(1).mf(j).params(2) = bestSolutionMembershipFunctions{3}(j, 2);
    fis.output(1).mf(j).params(3) = bestSolutionMembershipFunctions{3}(j, 3);
    
end

sim(sim_file_name);

f = figure();
plot(time, y); xlabel('t[s]'); ylabel('y(t)'); title('Odskocni odziv PD-FLC'); 
grid on;

figure(f);
print -deps PD_Fuzzy;

savefig('PD_FLC.fig');

steady_state = y(end);

k = find(y >= 0.9*steady_state);
Tr = time(k(1));
k = find(y >= 0.98*steady_state);
for i = 0 : length(k)-1
    if(i == length(k)-1) 
        Ts = time(k(end));
        break
    end
    if k(length(k) - i) ~= k(length(k) -i -1)
        Ts = time(k(i+1));
        break
    end
end
P = (max(y) - steady_state)/steady_state;
Peak = max(y);



disp('Ref. Time-Span No.Osc. Amplitude     Tr         P[%]         Ts        Peak');

disp(sprintf('                              %3.4f     %10.4f      %10.4f      %10.4f',Tr, P*100, Ts, Peak));
writefis(fis, 'final fis');

close_system(sim_file_name);

    
       