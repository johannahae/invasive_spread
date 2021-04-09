
g++ -std=c++11 make_webs.cpp -I/usr/local/include/ -L/usr/local/lib/ -o webs -lm -lgsl -lgslcblas

g++ invasive_spread.cpp -I/usr/local/include/ -L/usr/local/lib/ -o simulation -lm -lgsl -lgslcblas -lsundials_cvode -lsundials_nvecserial
