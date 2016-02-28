add_library(reconstructionfit SHARED reconstruction_fit.cpp)
target_link_libraries(reconstructionfit FitGen ${ROOT_LIBRARIES})