/***
 *  Project: TreeWorks
 *
 *  File: user_read_input_data.hpp
 *  Created: Mar 17, 2009
 *
 *  Author: Abhinav Sarje <abhinav.sarje@gmail.com>
 */

#ifndef USER_READ_INPUT_DATA_HPP
#define USER_READ_INPUT_DATA_HPP

#include "mpi_data.hpp"
#include "user_input.hpp"
#include "const_data.hpp"

bool read_point_data(const char* input_file, const MPI_data& mpi_data, InputData& data);

#endif // USER_READ_INPUT_DATA_HPP
