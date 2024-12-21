#!/bin/bash
find . -type f \( -name "*.cpp" -o -name "*.cxx" -o -name "*.C" -o -name "*.h" \) \
    ! -path "./macros" ! -path "./sofia" ! -path "./frs" ! -path "./asyeos" ! -path "./glad-tpc" -print0 |
    xargs -0 -L 1 clang-format-15 -i

echo "Use git add -A ; git commit -m \"clang-format all files\" --author=\"white space <whitespace@example.com>\" to commit changes."
