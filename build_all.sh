#!/bin/bash

# Pastikan path ANDROID_HOME sudah benar sesuai instalasi Arch Linux Anda
export ANDROID_HOME=${ANDROID_HOME:-/opt/android-sdk}
NDK_VERSION="25.2.9519653"
TOOLCHAIN="$ANDROID_HOME/ndk/$NDK_VERSION/build/cmake/android.toolchain.cmake"

# Daftar arsitektur target Android yang Anda minta
ABIS=("armeabi-v7a" "arm64-v8a" "x86" "x86_64")

for ABI in "${ABIS[@]}"; do
    BUILD_DIR="build_$ABI"
    echo "========================================="
    echo "Memulai build untuk arsitektur: $ABI"
    echo "========================================="
    
    mkdir -p "$BUILD_DIR" && cd "$BUILD_DIR" || exit 1
    
    # Pengaturan khusus tambahan untuk perangkat 32-bit (v7) agar perhitungan matematika lebih cepat
    EXTRA_FLAGS=""
    if [ "$ABI" == "armeabi-v7a" ]; then
        EXTRA_FLAGS="-DANDROID_ARM_NEON=ON"
    fi

    # Konfigurasi CMake menggunakan toolchain NDK
    cmake .. \
      -DCMAKE_TOOLCHAIN_FILE="$TOOLCHAIN" \
      -DANDROID_ABI="$ABI" \
      -DANDROID_PLATFORM=android-21 \
      $EXTRA_FLAGS

    # Proses kompilasi menggunakan semua core CPU yang tersedia
    make -j$(nproc)
    
    # Kembali ke folder root proyek
    cd ..
    echo "Selesai build untuk arsitektur: $ABI"
    echo ""
done

echo "Seluruh arsitektur berhasil dikompilasi!"

