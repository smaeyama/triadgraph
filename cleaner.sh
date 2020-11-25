#!/bin/sh

find ./png/ -name *.png  | xargs rm -f
find ./dot/ -name *.dot  | xargs rm -f
find . -type d -name .ipynb_checkpoints | xargs rm -rf
find . -type d -name __pycache__ | xargs rm -rf

jupyter nbconvert --ClearOutputPreprocessor.enabled=True --inplace *.ipynb
jupyter nbconvert --to python *.ipynb


