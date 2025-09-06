#!/bin/bash

# --- 이 스크립트는 rxnmech 프로젝트를 빌드하고 실행합니다. ---
# 사용법:
#    ./build.sh                  (기본값: 릴리즈(-O3) 모드로 빌드 및 실행)
#    ./build.sh release          (릴리즈(-O3) 모드로 빌드 및 실행)
#    ./build.sh debug            (디버그 모드로 빌드 및 실행)
#    ./build.sh debug 50 40 30   (디버그 모드로 빌드 후, 인자를 전달하여 실행)

# --- 변수 설정 ---
BUILD_DIR="build"
# macOS에서는 기본적으로 "Unix Makefiles" 또는 "Ninja" 등을 사용합니다.
# "Ninja"가 일반적으로 더 빠르므로, Ninja가 설치되어 있다면 이를 사용하는 것이 좋습니다.
# CMAKE_GENERATOR="Unix Makefiles"
CMAKE_GENERATOR="Ninja" # Ninja가 설치되어 있다면 이 줄을 활성화하세요.

# 기본 빌드 타입을 Release로 설정
BUILD_TYPE="Release"

# --- 첫 번째 인자로 빌드 타입 결정 ---
if [ "$1" == "debug" ]; then
    BUILD_TYPE="Debug"
    shift # 첫 번째 인자 제거 (나머지 인자를 앱에 전달하기 위함)
elif [ "$1" == "release" ]; then
    # 'release'는 기본값이므로, 이 인자는 그냥 무시하고 넘어감
    shift # 첫 번째 인자 제거
fi

echo "================================================="
echo "                  BUILD MODE: ${BUILD_TYPE}"
echo "================================================="

# --- 빌드 디렉토리 확인 및 생성 ---
if [ ! -d "${BUILD_DIR}" ]; then
    echo "Creating build directory: ${BUILD_DIR}"
    mkdir -p "${BUILD_DIR}"
fi

# --- CMake 실행 ---
echo ""
echo "Running CMake to configure the project..."
# 현재 디렉토리를 기준으로 빌드 디렉토리로 이동
(cd "${BUILD_DIR}" && cmake .. -G "${CMAKE_GENERATOR}" -D CMAKE_BUILD_TYPE="${BUILD_TYPE}")

# CMake 실행 결과 확인
if [ $? -ne 0 ]; then
    echo ""
    echo "CMake configuration failed. Aborting."
    exit 1
fi

# --- 컴파일 실행 (make 또는 ninja) ---
echo ""
echo "Compiling the project with ${CMAKE_GENERATOR}..."
# `make`는 `Unix Makefiles`를 사용할 때, `ninja`는 `Ninja`를 사용할 때
# macOS에서는 기본적으로 `make`가 설치되어 있으나, `ninja`가 더 빠를 수 있습니다.
# -j 옵션은 코어 수에 따라 병렬 빌드를 수행합니다.
# `nproc` 또는 `sysctl -n hw.ncpu`로 코어 수를 가져올 수 있습니다.
NUM_PROCESSORS=$(sysctl -n hw.ncpu)

if [ "${CMAKE_GENERATOR}" == "Ninja" ]; then
    (cd "${BUILD_DIR}" && ninja -j"${NUM_PROCESSORS}")
else
    (cd "${BUILD_DIR}" && make -j"${NUM_PROCESSORS}")
fi


# 컴파일 결과 확인
if [ $? -ne 0 ]; then
    echo ""
    echo "Compilation failed. Aborting."
    exit 1
fi

# --- 실행 ---
echo ""
echo "Build successful."

# build.sh에 전달된 나머지 모든 인자를 $@ 로 받아서 실행 파일에 그대로 전달
# 실행 파일 경로는 CMakeLists.txt에 정의된 프로젝트 이름에 따라 달라질 수 있습니다.
# 여기서는 예시로 rxnmech_app을 사용합니다.
APP_NAME="rxnmech_app" # CMakeLists.txt의 add_executable 이름과 일치해야 합니다.
APP_PATH="${BUILD_DIR}/${APP_NAME}"

echo "Running the application: ${APP_PATH} $@"
echo "--------------------------------------------------"

# 애플리케이션 실행
"${APP_PATH}" "$@"

echo "--------------------------------------------------"
echo ""
echo "Application finished."
# Windows의 'pause'와 달리, macOS에서는 명시적으로 기다릴 필요가 없는 경우가 많습니다.
# 필요하다면 다음 줄의 주석을 해제하세요.
# read -p "Press any key to continue..." -n1 -s