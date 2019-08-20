#!/bin/bash
if [[ $# -ne 2 ]]
then
	echo "Usage: ./update_tmpmaker.sh <destination> <test_source_dir>"
	exit 1
fi

dest_dir=$(readlink -f $1)
test_scripts=$(readlink -f $2)

if [ ! -d "$dest_dir" ]
then
	echo "Directory $dest_dir does not already exist, please create with tmpmaker.sh script"
	exit 1
fi

cd $dest_dir
echo "Installing psamplers"
Rscript -e "devtools::install_github('gjcooper/samplers', auth_token='c71394bda6a42db165ee8b114dbfa1688ed5882f')"
cp -a $test_scripts/* .
echo "Try running Rscript PMwG.R from $dest_dir directory now..."
