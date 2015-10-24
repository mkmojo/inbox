###
#  Project: TreeWorks
#
#  File: Makefile
#  Created: Mar 15, 2009
#
#  Author: Abhinav Sarje <abhinav.sarje@gmail.com>
###

CXX = mpicxx
CXXFLAGS = -O3 -std=c++0x #-mpe=mpilog
CXXLIBFLAGS = 
CXXIFLAGS =

test: user_test.o user_read_input_data.o
	$(CXX) -o test $(CXXFLAGS) $(CXXIFLAGS) user_test.o user_read_input_data.o $(CXXLIBFLAGS)

user_read_input_data.o: user_read_input_data.cpp user_read_input_data.hpp input_data.hpp mpi_data.hpp sample_sort.hpp user_head.hpp
	$(CXX) -c -o user_read_input_data.o $(CXXFLAGS) $(CXXIFLAGS) user_read_input_data.cpp $(CXXLIBFLAGS)

user_test.o: user_test.cpp user_combine_function.hpp user_generate_function.hpp base_tree.hpp tree_construct_octree.hpp sample_sort.hpp utility.hpp adjacent_nodes_generate_function.hpp interaction_list_generate_function.hpp tree_compute.hpp specialized_generate.hpp knn_generate.hpp knn_combine.hpp user_head.hpp user_read_input_data.o 
	$(CXX) -c -o user_test.o $(CXXFLAGS) $(CXXIFLAGS) user_test.cpp $(CXXLIBFLAGS)

clean:
	rm -f *.o test
