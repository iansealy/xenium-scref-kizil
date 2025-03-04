#!/usr/bin/env bash

for version in `cut -f1 data/versions.txt`; do

# Skip if version already done
if [ -e "data/$version/ensembl.txt" ]; then
  continue
fi

# Make directory for genes
mkdir -p data/$version

# Clone BioPerl repository if haven't already
if [ ! -e "bioperl-live" ]; then
  git clone -q https://github.com/bioperl/bioperl-live.git
  cd bioperl-live
  git checkout -q bioperl-release-1-6-9
  cd ..
fi

# Clone Ensembl repository if haven't already
if [ ! -e "ensembl" ]; then
  git clone -q https://github.com/Ensembl/ensembl.git
fi

# Change Ensembl version
cd ensembl
git pull -q
git checkout -q release/$version
cd ..

# Dump Ensembl genes
perl -I ensembl/modules -I bioperl-live script/dump_ensembl_genes.pl \
> data/$version/ensembl.txt.tmp \
2> data/$version/ensembl.stdout.tmp
mv data/$version/ensembl.txt.tmp data/$version/ensembl.txt
mv data/$version/ensembl.stdout.tmp data/$version/ensembl.stdout

done

# Remove Ensembl and BioPerl repositories
rm -rf bioperl-live ensembl

# Check for errors
errors=""
for version in `cut -f1 data/versions.txt`; do
  count=`cat "data/$version/ensembl.stdout" | grep -cv '^#'`
  if [ $count -ne 0 ]; then
    errors="$errors e$version"
  fi
done
if [ ! -z "$errors" ]; then
  echo "Versions with errors:"
  echo $errors
  exit 1
fi
