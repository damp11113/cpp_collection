# TJC Studio: GUI player/encoder for Windows and Linux (X11).
#
#   cmake -S . -B build && cmake --build build --config Release
#
# GLFW, Dear ImGui and miniaudio are downloaded at configure time. For an
# offline build, point FETCHCONTENT_SOURCE_DIR_GLFW / _IMGUI / _MINIAUDIO at
# local copies.
cmake_minimum_required(VERSION 3.16)
project(tjc_studio C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release)
endif()
if(POLICY CMP0135)
    cmake_policy(SET CMP0135 NEW)
endif()

include(FetchContent)
set(GLFW_BUILD_DOCS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(GLFW_INSTALL OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_WAYLAND OFF CACHE BOOL "" FORCE)  # X11 (also runs under XWayland)
FetchContent_Declare(glfw
    URL https://github.com/glfw/glfw/releases/download/3.4/glfw-3.4.zip
    URL_HASH SHA256=b5ec004b2712fd08e8861dc271428f048775200a2df719ccf575143ba749a3e9)
# ImGui and miniaudio are used as plain sources; SOURCE_SUBDIR points at a
# directory without a CMakeLists.txt so their own build files are not run.
FetchContent_Declare(imgui
    URL https://github.com/ocornut/imgui/archive/refs/tags/v1.91.9b.tar.gz
    SOURCE_SUBDIR no-cmake)
FetchContent_Declare(miniaudio
    URL https://github.com/mackron/miniaudio/archive/refs/tags/0.11.22.tar.gz
    SOURCE_SUBDIR no-cmake)
FetchContent_MakeAvailable(glfw imgui miniaudio)

find_package(OpenGL REQUIRED)
find_package(Threads REQUIRED)

add_executable(tjc_studio WIN32
    main.cpp
    player.cpp
    encode_job.cpp
    audio_out.cpp
    platform.cpp
    file_browser.cpp
    tjc_impl.cpp
    ${imgui_SOURCE_DIR}/imgui.cpp
    ${imgui_SOURCE_DIR}/imgui_draw.cpp
    ${imgui_SOURCE_DIR}/imgui_tables.cpp
    ${imgui_SOURCE_DIR}/imgui_widgets.cpp
    ${imgui_SOURCE_DIR}/backends/imgui_impl_glfw.cpp
    ${imgui_SOURCE_DIR}/backends/imgui_impl_opengl2.cpp)
target_include_directories(tjc_studio PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/..
    ${imgui_SOURCE_DIR}
    ${imgui_SOURCE_DIR}/backends
    ${miniaudio_SOURCE_DIR})
target_link_libraries(tjc_studio PRIVATE glfw OpenGL::GL Threads::Threads ${CMAKE_DL_LIBS})
target_compile_definitions(tjc_studio PRIVATE _FILE_OFFSET_BITS=64)

if(WIN32)
    target_link_libraries(tjc_studio PRIVATE comdlg32 shell32)
    target_compile_definitions(tjc_studio PRIVATE UNICODE _UNICODE _WIN32_WINNT=0x0601 _CRT_SECURE_NO_WARNINGS)
    if(MSVC)
        # WIN32 executable with a regular main().
        target_link_options(tjc_studio PRIVATE /ENTRY:mainCRTStartup)
        target_compile_options(tjc_studio PRIVATE /utf-8)
    elseif(MINGW)
        option(TJC_STUDIO_STATIC "Link the MinGW runtime statically (no extra DLLs needed)" ON)
        if(TJC_STUDIO_STATIC)
            target_link_options(tjc_studio PRIVATE -static)
        endif()
    endif()
else()
    target_link_libraries(tjc_studio PRIVATE m)
endif()
