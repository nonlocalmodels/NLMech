#!/bin/bash

##############################################################################
#  Copyright (c) 2019 Prashant K. Jha
#  Copyright (c) 2019 Patrick Diehl
#
#  Distributed under the Boost Software License, Version 1.0. (See accompanying
#  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
##############################################################################

branch=$(git rev-parse --abbrev-ref HEAD)
echo "Checking branch ${branch}"

list="$(git diff --name-only ${branch})"

if [ "${#list}" -eq "0"  ]; then
  echo "No files to check"
  exit
fi

files=$(echo $list | tr " " "\n")

for file in $files
do
  echo "Checking file ${file}"
  aspell --mode=ccpp check -l En_us $file
done


