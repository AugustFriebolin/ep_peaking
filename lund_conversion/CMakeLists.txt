# Add executable for converting generator output to lund format.
add_executable(lund_convert lund_convert.cpp )
target_link_libraries(lund_convert ${ROOT_LIBRARIES} )
