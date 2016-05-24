  # Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project
  # Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
  #   Max-Planck-Institute for Polymer Research, Theory Group
  
  # This file is part of ESPResSo.
  
  # ESPResSo is free software: you can redistribute it and/or modify
  # it under the terms of the GNU General Public License as published by
  # the Free Software Foundation, either version 3 of the License, or
  # (at your option) any later version.
  
  # ESPResSo is distributed in the hope that it will be useful,
  # but WITHOUT ANY WARRANTY; without even the implied warranty of
  # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  # GNU General Public License for more details.
  
  # You should have received a copy of the GNU General Public License
  # along with this program.  If not, see <http://www.gnu.org/licenses/>. 

# Target for the unit tests
add_custom_target(check_unit_tests COMMAND ${CMAKE_CTEST_COMMAND})

# Run unit tests on check
add_dependencies(check check_unit_tests)


# Vector
add_executable(Vector_test Vector_test.cpp)
target_link_libraries(Vector_test ${LIBRARIES})
add_test(vector Vector_test)
add_dependencies(check_unit_tests Vector_test)

# RuntimeError
add_executable(RuntimeError_test RuntimeError_test.cpp)
target_link_libraries(RuntimeError_test ${LIBRARIES})
add_test(runtimeerror RuntimeError_test)
add_dependencies(check_unit_tests RuntimeError_test)

# RuntimeErrorCollector
add_executable(RuntimeErrorCollector_test RuntimeErrorCollector_test.cpp ../RuntimeErrorCollector.cpp ../RuntimeError.cpp)
target_link_libraries(RuntimeErrorCollector_test ${LIBRARIES})
add_test(runtimeerrorcollector_test RuntimeErrorCollector_test)
add_dependencies(check_unit_tests RuntimeErrorCollector_test)

# RunningAverage
add_executable(RunningAverage_test RunningAverage_test.cpp)
target_link_libraries(RunningAverage_test ${LIBRARIES})
add_test(runing_average RunningAverage_test)

