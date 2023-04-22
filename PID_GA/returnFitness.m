function foundIAE = returnFitness(sim_file_name)

IAE = 0;


open_system(sim_file_name);



sim(sim_file_name);



V=IAE(end); 
foundIAE=V;
end