rm final_data.dat
ls Brahe-3141-f/* > file_names.dat
cc execute.c 
./a.out 
chmod u+x run_files.sh
cc Least_square.c
./run_files.sh
rm run_files.sh
rm file_names.dat