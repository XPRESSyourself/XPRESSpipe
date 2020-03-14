OS=$(uname -s)

echo ${OS}

if [ ${OS} == "Darwin" ]
then
    make -f Makefile_macOS
else
    make -f Makefile_Linux
fi
