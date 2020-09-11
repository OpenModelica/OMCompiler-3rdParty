cmake_minimum_required(VERSION 3.14)
project(lpsolve)


find_package(BISON REQUIRED)
find_package(FLEX REQUIRED)


set(LP_SOURCES lp_MDO.c shared/commonlib.c shared/mmio.c
	            shared/myblas.c	ini.c fortify.c	colamd/colamd.c
	            lp_crash.c bfp/bfp_LUSOL/lp_LUSOL.c
				bfp/bfp_LUSOL/LUSOL/lusol.c lp_Hash.c
				lp_lib.c lp_wlp.c lp_matrix.c
				lp_mipbb.c lp_MPS.c lp_params.c
				lp_presolve.c lp_price.c lp_pricePSE.c
				lp_report.c lp_scale.c lp_simplex.c
				lp_SOS.c lp_utils.c yacc_read.c
	)




# Put the generated files in a specific directory to simplify debugging.
set (GENERATED_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/generated)
file(MAKE_DIRECTORY ${GENERATED_DIRECTORY})

############################ process lp_rlp.y ###################################
# Add a bison compilation of lp_rlp.y. The generated C file will be modified afterwards
# to change yy to lp_yy in teh code. See next custom command
BISON_TARGET(lp_rlp_y # name
				${CMAKE_CURRENT_SOURCE_DIR}/lp_rlp.y # input
				${GENERATED_DIRECTORY}/y.tab.c # output a.k.a BISON_lp_rlp_y_OUTPUT_SOURCE
				COMPILE_FLAGS "--no-lines -y -d"
			 )
# change yy to lp_yy and put the result in lp_rlp.c
add_custom_command(OUTPUT ${GENERATED_DIRECTORY}/lp_rlp.c
				   COMMAND sed -e "s/yy/lp_yy/g" ${BISON_lp_rlp_y_OUTPUT_SOURCE} > ${GENERATED_DIRECTORY}/lp_rlp.c
				   DEPENDS ${BISON_lp_rlp_y_OUTPUT_SOURCE}
				   COMMENT "Generating ${GENERATED_DIRECTORY}/lp_rlp.c from ${BISON_lp_rlp_y_OUTPUT_SOURCE}  (change yy to lp_yy)")

# set_source_files_properties(${GENERATED_DIRECTORY}/lp_rlp.c GENERATED)
set(LP_SOURCES ${LP_SOURCES} ${GENERATED_DIRECTORY}/lp_rlp.c)
############################ process lp_rlp.y ###################################


############################ process lp_rlp.l ###################################
# Add a flex compilation of lp_rlp.l. The generated C file will be modified afterwards
# to change yy to lp_yy in the code. See next custom command
FLEX_TARGET(lp_rlp_l # name
				${CMAKE_CURRENT_SOURCE_DIR}/lp_rlp.l # input
				${GENERATED_DIRECTORY}/lex.yy.c # output , there is no FLEX_lp_rlp_l_OUTPUT_SOURCE unfortunately.
				COMPILE_FLAGS "-L"
			 )
# change yy to lp_yy and put the result in lp_rlp.h
add_custom_command(OUTPUT ${GENERATED_DIRECTORY}/lp_rlp.h
				   COMMAND sed -e "s/yy/lp_yy/g" ${GENERATED_DIRECTORY}/lex.yy.c > ${GENERATED_DIRECTORY}/lp_rlp.h
				   DEPENDS ${GENERATED_DIRECTORY}/lex.yy.c
				   COMMENT "Generating ${GENERATED_DIRECTORY}/lp_rlp.h from ${GENERATED_DIRECTORY}/lex.yy.c  (change yy to lp_yy)")

# set_source_files_properties(${GENERATED_DIRECTORY}/lp_rlp.c GENERATED)
# This is important! One it causes the custom command above to run.
# It also means we do not need to add the 'generated' directory as include directory
# just for this file.
set(LP_SOURCES ${LP_SOURCES} ${GENERATED_DIRECTORY}/lp_rlp.h)
############################ process lp_rlp.y ###################################



add_library(lpsolve55 SHARED ${LP_SOURCES})
# This is a bad way to include things. But we will not modify the lpsolve code yet.
# For other targets that depend on lpsolve we use the custom target below to copy
# tre needed headers to one place and we provide that as an interface include directory.
target_include_directories(lpsolve55 PRIVATE ./ ./shared ./bfp ./bfp/bfp_LUSOL ./bfp/bfp_LUSOL/LUSOL ./colamd)

target_compile_definitions(lpsolve55 PRIVATE YY_NEVER_INTERACTIVE PARSER_LP INVERSE_ACTIVE=INVERSE_LUSOL RoleIsExternalInvEngine)



# Copy the "public" headers to one place and provide that for targets
# that depend on lpsolve.
file(MAKE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/lpsolve/)
add_custom_command(TARGET lpsolve55
                   POST_BUILD
				   COMMAND ${CMAKE_COMMAND} -E copy
						${CMAKE_CURRENT_SOURCE_DIR}/lp_lib.h
						${CMAKE_CURRENT_SOURCE_DIR}/lp_types.h
						${CMAKE_CURRENT_SOURCE_DIR}/lp_utils.h
						${CMAKE_CURRENT_SOURCE_DIR}/lp_Hash.h
						${CMAKE_CURRENT_SOURCE_DIR}/lp_matrix.h
						${CMAKE_CURRENT_SOURCE_DIR}/lp_mipbb.h
						${CMAKE_CURRENT_SOURCE_DIR}/lp_SOS.h
						${CMAKE_CURRENT_SOURCE_DIR}/include/lpsolve/
				   COMMENT "Copying lpsolve header files to ${CMAKE_CURRENT_SOURCE_DIR}/include/"
				   )
target_include_directories(lpsolve55 INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include/)

