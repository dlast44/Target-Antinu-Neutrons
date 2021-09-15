#add_library(event CVUniverse.cpp)
add_library(event NeutCands.cpp)
target_link_libraries(event ${ROOT_LIBRARIES})
install(TARGETS event DESTINATION lib)
