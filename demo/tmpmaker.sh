#!/bin/bash
if [[ $# -ne 2 ]]
then
	echo "Usage: ./tmpmaker.sh <destination> <test_source_dir>"
	exit 1
fi

dest_dir=$(readlink -f $1)
test_scripts=$(readlink -f $2)

if [ -d "$dest_dir" ]
then
	echo "Directory $dest_dir already exists, please delete or choose a new directory"
	exit 1
fi

echo "Creating test directory $dest_dir"
mkdir $dest_dir
cd $dest_dir
echo "Initialising Environment"
Rscript -e "packrat::init()"
sed -i 's/external.packages: /external.packages: devtools/g' packrat/packrat.opts
echo "Installing psamplers"
Rscript -e "devtools::install_github('gjcooper/samplers', auth_token='c71394bda6a42db165ee8b114dbfa1688ed5882f')"
cp -a $test_scripts/* .
echo "Try running Rscript PMwG.R from $dest_dir directory now..."
