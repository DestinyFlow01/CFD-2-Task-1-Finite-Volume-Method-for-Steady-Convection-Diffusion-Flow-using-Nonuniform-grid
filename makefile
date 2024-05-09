file_obj = FVM.o main.o
file_cpp = FVM.cpp main.cpp

run :
	make allclean
	g++ -c $(file_cpp)
	g++ $(file_obj) -o main.exe
	./main.exe > "outputFVM.txt"
allclean : 
	rm -f *.csv *.o *.exe

