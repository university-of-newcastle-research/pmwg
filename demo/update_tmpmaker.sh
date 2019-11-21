#!/bin/bash
if [[ $# -ne 3 ]]
then
	echo "Usage: ./update_tmpmaker.sh (-l(local)/-g(github)) <destination> <test_source_dir>"
	exit 1
fi

loc=0
if [[ $1 -eq '-l' ]]
then
	echo "Grabbing from local"
	loc=1
fi
dest_dir=$(readlink -f $2)
test_scripts=$(readlink -f $3)

if [ ! -d "$dest_dir" ]
then
	echo "Directory $dest_dir does not already exist, please create with tmpmaker.sh script"
	exit 1
fi

cd $dest_dir
echo "Installing psamplers"
if [[ $loc -eq 1 ]]
then
	Rscript -e "devtools::install_local('$test_scripts/..', force=TRUE)"
else
	Rscript -e "devtools::install_github('gjcooper/samplers')"
fi
cp -a $test_scripts/* .
echo "Try running Rscript PMwG.R from $dest_dir directory now..."
