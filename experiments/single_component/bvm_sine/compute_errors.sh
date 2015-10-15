g++ compute_errors.cpp -lboost_program_options -lboost_system -o compute_errors

./compute_errors --n 10 --prior 2
#./compute_errors --n 20 --prior 2
#./compute_errors --n 30 --prior 2
#./compute_errors --n 40 --prior 2
./compute_errors --n 50 --prior 2
./compute_errors --n 100 --prior 2

./compute_errors --n 10 --prior 3
#./compute_errors --n 20 --prior 3
#./compute_errors --n 30 --prior 3
#./compute_errors --n 40 --prior 3
./compute_errors --n 50 --prior 3
./compute_errors --n 100 --prior 3
