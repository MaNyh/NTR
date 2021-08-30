currentOutFile=$(find . -type f -name '*.out' -printf '%T@ %p\n' | sort -n | tail -1 | cut -f2- -d" ")
#tail -f $currentOutFile
pimpleIterations=$(cut -f 1 $currentOutFile | grep -i "PIMPLE: iteration 1" | wc -l)
echo "iteration number of last simulation:"
echo $pimpleIterations

