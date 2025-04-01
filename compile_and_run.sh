#!/bin/bash

g++ -std=c++11 -o cube_raytracer cube_raytracer.cpp

if [ $? -eq 0 ]; then
    echo "Compilation successful. Running the program..."
    ./cube_raytracer
else
    echo "Compilation failed. Please check for errors."
    exit 1
fi