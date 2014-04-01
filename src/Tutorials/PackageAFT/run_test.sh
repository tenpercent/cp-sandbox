#!/bin/bash

[[ "x$1" == "x" ]] && FILE=../data/sculpture.frt || FILE=$1

echo "Testing $FILE"
echo "  cf nf_aft F         q_avg         q_min   R   C   message"
for cf in 1.05 1.10 1.15 1.20 1.25 1.30 1.35 1.40 1.45 1.50 1.55 1.60
do
    printf "%s:" $cf
    MSG=$(nice ./FrontToTet_aft.exe $FILE $cf 2>&1 | egrep "^(AFT returned|AFTINFO|MESH|QS:|ln |libaft\(3d\) warning: Check:)")
    NF=$(echo "$MSG" | grep "AFTINFO" | egrep -o "[0-9]+")
    echo "$MSG" | grep -q "MESH: ok" && OK=" " || OK="x"
    printf " %5d " $NF
    echo -n "$OK"
    QSFINAL=$(echo "$MSG" | sed -ne "/ln /{s/ln  q_min: \+\([0-9e.+-]\+\) \+#.*/\1/;p};/QS:/{s/QS: q_avg: \+\([0-9e.+-]\+\) \+#.*/\1 /;p}" | tail -2)
    printf "  %12.6le %13.6le  " $QSFINAL
    echo "$MSG" | grep -q "AFT returned: ok" && AFT=" . " || AFT="[X]"
    echo "$MSG" | grep -q "MESH: ok" && MESH=" . " || MESH="[X]"
    RES="unknown"
    echo "$MSG" | grep -q "tetras with negative volume" && RES="FAIL (negative volume)"
    echo "$MSG" | grep -q "AFT returned: FAIL" && RES="FAIL (internal stop)"
    echo "$MSG" | grep -q "Tetra does not have neighbour" && RES="FAIL (wrong topology)"
    echo "$MSG" | grep -q "Duplicate tetras" && RES="FAIL (wrong topology)"
    echo "$MSG" | grep -q "MESH: ok" && RES="pass"
    echo "$AFT $MESH  $RES"
done
echo
echo
