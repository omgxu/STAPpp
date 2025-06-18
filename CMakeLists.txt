cmake_minimum_required(VERSION 3.15)  # 推荐使用较新版本
project(stap++
  VERSION 1.0
  LANGUAGES CXX
  DESCRIPTION "Structural Analysis Program in C++"
)

# 1. 编译器标准设置（更严格的设置）
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)  # 禁用编译器扩展

# 2. 调试选项（更规范的现代写法）
option(STAP++_DEBUG "Enable debug output and symbols" OFF)
if(STAP++_DEBUG)
  add_compile_definitions(_DEBUG_)
  set(CMAKE_BUILD_TYPE Debug)  # 自动设置调试标志
endif()

# 3. 项目结构定义
set(SOURCE_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(CPP_DIR ${SOURCE_ROOT}/cpp)
set(INCLUDE_DIR ${SOURCE_ROOT}/h)

# 4. 源文件收集（显式列出更可靠）
set(SRC_FILES
  ${CPP_DIR}/Clock.cpp
  ${CPP_DIR}/Domain.cpp
  ${CPP_DIR}/ElementGroup.cpp
  ${CPP_DIR}/Elements/Bar.cpp
  ${CPP_DIR}/Elements/H8.cpp
  ${CPP_DIR}/Elements/MindlinPlate.cpp
  ${CPP_DIR}/Elements/T3.cpp
  ${CPP_DIR}/LoadCaseData.cpp
  ${CPP_DIR}/Material.cpp
  ${CPP_DIR}/Node.cpp
  ${CPP_DIR}/Outputter.cpp
  ${CPP_DIR}/Solver.cpp
  ${CPP_DIR}/main.cpp
)

# 5. 可执行目标配置
add_executable(${PROJECT_NAME} ${SRC_FILES})

# 6. 包含目录（现代写法）
target_include_directories(${PROJECT_NAME} PRIVATE
  ${INCLUDE_DIR}
  ${INCLUDE_DIR}/Elements
)

# 7. 编译特性（可选但推荐）
target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_11)

# 8. 安装规则（可选）
install(TARGETS ${PROJECT_NAME}
  RUNTIME DESTINATION bin
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
)

# 9. 单元测试支持（可选）
# enable_testing()
# add_test(NAME ${PROJECT_NAME}_test COMMAND ${PROJECT_NAME})

# 10. 包依赖查找（示例）
# find_package(Boost 1.70 REQUIRED COMPONENTS system filesystem)
# target_link_libraries(${PROJECT_NAME} PRIVATE Boost::system Boost::filesystem)

# 11. 调试信息输出
message(STATUS "Source files: ${SRC_FILES}")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Install prefix: ${CMAKE_INSTALL_PREFIX}")