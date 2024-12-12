mkdir Intervention6
cd Intervention6
rm -r *_*
cd ../
cat run_numbers | parallel --jobs 12 'echo {}; mkdir {}_1; cd {}_1 ; /Users/christopher.illingworth/Documents/Coronavirus/Simulation/UDCA_Project/WriteupCI/Code/run_sim_process --sim {} --verb 1 --intervention 6 > Run{}.out ; cd ../'
cat run_numbers | parallel --jobs 12 'echo {}; mkdir {}_2; cd {}_2 ; /Users/christopher.illingworth/Documents/Coronavirus/Simulation/UDCA_Project/WriteupCI/Code/run_sim_process --sim {} --verb 1 --intervention 6 > Run{}.out ; cd ../'
cat run_numbers | parallel --jobs 12 'echo {}; mkdir {}_3; cd {}_3 ; /Users/christopher.illingworth/Documents/Coronavirus/Simulation/UDCA_Project/WriteupCI/Code/run_sim_process --sim {} --verb 1 --intervention 6 > Run{}.out ; cd ../'
cat run_numbers | parallel --jobs 12 'echo {}; mkdir {}_4; cd {}_4 ; /Users/christopher.illingworth/Documents/Coronavirus/Simulation/UDCA_Project/WriteupCI/Code/run_sim_process --sim {} --verb 1 --intervention 6 > Run{}.out ; cd ../'
cat run_numbers | parallel --jobs 12 'echo {}; mkdir {}_5; cd {}_5 ; /Users/christopher.illingworth/Documents/Coronavirus/Simulation/UDCA_Project/WriteupCI/Code/run_sim_process --sim {} --verb 1 --intervention 6 > Run{}.out ; cd ../'
mv ?_? ??_? ???_? Intervention6
