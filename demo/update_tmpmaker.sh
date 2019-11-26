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
cd $test_scripts
echo "Copying parallel and single core demo files"
git show parallelDemo:demo/parallelPMwG.R > $dest_dir/parallelPMwG.R
git show singleDemo:demo/singlePMwG.R > $dest_dir/singlePMwG.R
cd $dest_dir
echo "Try running Rscript <master|single|parallel>PMwG.R from $dest_dir directory now..."
