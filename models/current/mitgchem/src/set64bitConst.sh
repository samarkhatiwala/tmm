#! /usr/bin/env sh
#
# 64-bit single precision
sym64bitConst=E
#
# 32-bit single precision
sym64bitConst=D
sed s'/ * _d  */'${sym64bitConst}'/g'
