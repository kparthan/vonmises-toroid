g++ compute_kldivs.cpp -lboost_program_options -lboost_system -o compute_kldivs

./compute_kldivs --n 10 --prior 2
#./compute_kldivs --n 20 --prior 2
#./compute_kldivs --n 30 --prior 2
#./compute_kldivs --n 40 --prior 2
./compute_kldivs --n 50 --prior 2
./compute_kldivs --n 100 --prior 2

./compute_kldivs --n 10 --prior 3
#./compute_kldivs --n 20 --prior 3
#./compute_kldivs --n 30 --prior 3
#./compute_kldivs --n 40 --prior 3
./compute_kldivs --n 50 --prior 3
./compute_kldivs --n 100 --prior 3
