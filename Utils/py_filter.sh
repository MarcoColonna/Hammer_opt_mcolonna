#!/bin/bash
if [ -x "$(command -v doxypypy)" ]; then
   doxypypy -a -c $1
fi
