#! /usr/bin/env sh
#
sed 's/^[Cc]/\/\//g' | sed 's/\(^.*CPP_OPTIONS.*$\)/\/\/\1/'