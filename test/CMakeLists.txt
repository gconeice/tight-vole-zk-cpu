#Testing macro
macro (add_test_executable_with_lib folder _name)
	add_executable(test_${folder}_${_name} "${folder}/${_name}.cpp")
	target_link_libraries(test_${folder}_${_name} ${EMP-OT_LIBRARIES} emp-zk)
endmacro()

macro (add_test_case folder _name)
	add_test_executable_with_lib(${folder} ${_name})
  	add_test(NAME ${folder}_${_name} COMMAND "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_${folder}_${_name}" WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/")
endmacro()

macro (add_test_case_with_run folder _name)
	add_test_executable_with_lib(${folder} ${_name})
	add_test(NAME ${folder}_${_name} COMMAND "./run" "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_${folder}_${_name}" WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/")

endmacro()

# Tetst Path ZK
add_test_case_with_run(arith stack_batched_matmul_v1)

add_test_case_with_run(pathzk test)
add_test_case_with_run(pathzk comp_with_pub_cir)
add_test_case_with_run(pathzk comp_with_pub_cir_opt)
add_test_case_with_run(pathzk pub_cir)

add_test_case_with_run(pathzk comp_batchman_balance)
add_test_case_with_run(pathzk comp_batchman_unbalance)
add_test_case_with_run(pathzk comp_batchman_uniform)
add_test_case_with_run(pathzk comp_batchman_balance_opt)

add_test_case_with_run(pathzk fine_grain)