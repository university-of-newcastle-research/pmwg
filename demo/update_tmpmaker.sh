#!/bin/bash
if [[ $# -ne 3 ]]
then
	echo "Usage: ./update_tmpmaker.sh <destination> <test_source_dir>"
	exit 1
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
Rscript -e "devtools::install_local('$test_scripts/..', force=TRUE)"
cp -a $test_scripts/* .
echo "Try running Rscript <master|single|parallel>PMwG.R from $dest_dir directory now..."
