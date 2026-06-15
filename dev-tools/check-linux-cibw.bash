#!/bin/bash

# Run cibuildwheel locally for Linux; useful as a pre-push check

CIBW_ENVIRONMENT="CMAKE_ARGS='-DSLYCOT_BUNDLE_OPENBLAS=ON'" CIBW_BUILD=cp313-manylinux_x86_64 cibuildwheel
