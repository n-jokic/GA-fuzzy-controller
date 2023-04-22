clear
close all
clc

nP = 30; %population size
rC = 0.7; %crossover probability
rM = 0.08; % mutation probability
maxIter = 50; %criterion for convergence
maxKp = 100;
maxK = 100;
maxKd = 50;

%Now we list parameters for our simulation
t_end = 10;
t_start = 0;
Num = [2];
Den = [1 12 24];

sim_file_name = 'PID_GA_simulink';

open_system(sim_file_name);




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

%Now we'll do the extraction part
KpPopulation = zeros(nP, numBitsCoef);
KPopulation = zeros(nP, numBitsCoef);
KdPopulation = zeros(nP, numBitsCoef);
FitnessPopulation = zeros(nP,1);

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
    






    sim(sim_file_name);


    display(IAE(end));

    fitness_vector(i) = IAE(end); 
end

%Now we will sort Ks going by their fitness value (note that for two
%fitness values the better one will be the one that's smaller - IAE
%criterion

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

%now let's start iterating through our generations
iter = 0;
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

        
        boundary_bit = floor(rand(1,1)*numBitsCoef) + 1;
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
        
        sim(sim_file_name);


        display(IAE(end));

        FitnessChildren(i) = IAE(end); 
        
    end
    
    %Now, we'll compare both our parents and pur children and we'll choose
    %nP fittest 
     
    KpInitial = MergeMatrix(KpPopulation, KpChildren);
    KInitial = MergeMatrix(KPopulation, KChildren);
    KdInitial = MergeMatrix(KPopulation, KChildren);
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
    
    iter = iter + 1;
end

%Now, when we have achieved our convergence criterion, we'll choose the
%best solution from our final population

bestSolutionKp = KpPopulation(1,:); 
bestSolutionK = KPopulation(1,:);
bestSolutionKd = KdPopulation(1,:); 

 
bestKp = 0; bestKd = 0; bestK = 0;
for j = 1 : numBitsCoef
        bestKp = bestKp + bestSolutionKp(j)*(0.5^(j-1));
        bestK = bestK + bestSolutionK(j)*(0.5^(j-1));
        bestKd = bestKd + bestSolutionKd(j)*(0.5^(j-1));
end
bestKp = bestKp*maxKp / K_correct
bestK = bestK*maxK / K_correct
bestKd = bestKd*maxKd / K_correct


K = bestK;
Kd = bestKd;
Kp = bestKp;



sim(sim_file_name);

f = figure();
plot(time, y); xlabel('t[s]'); ylabel('y(t)'); title('Odskocni odziv PID'); 
grid on;

figure(f);
print -deps PID;

savefig('PID.fig');

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


close_system(sim_file_name);