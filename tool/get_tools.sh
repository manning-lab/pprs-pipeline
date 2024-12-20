#!/bin/bash

pwd | grep -q "tool$" || { echo 'Please run this script from the tool/ folder.'; exit 1; }

install_bgenix() {
  echo 'Getting bgenix... (this takes a long time to compile)'
   if [ -d 'bgen.tgz' ]; then
     rm -r bgen.tgz # In case user got impatient and leaves a half-untarred directory
   fi
   curl -LO http://code.enkre.net/bgen/tarball/release/bgen.tgz
   tar zxf bgen.tgz
   mv bgen.tgz bgen # Extracted directory retains .tgz in its name for some reason
   cd bgen
   ./waf configure
   ./waf
   cd -
}

# Function to install bcftools
install_bcftools() {
  echo 'Getting bcftools...'
  git clone --depth 1 --recurse-submodules https://github.com/samtools/htslib.git
  git clone --depth 1 https://github.com/samtools/bcftools.git
  cd bcftools
  make -j
  cd -
}

# Function to install plink
install_plink() {
  echo 'Getting plink2...'
  curl -LO https://s3.amazonaws.com/plink2-assets/alpha6/plink2_linux_avx2_20241111.zip
  unzip -o plink2_linux_avx2_20241111.zip plink2
}


if [ $# -eq 0 ]; then
  install_bgenix
  install_bcftools
  install_plink
else
  case "$1" in
    "plink")    install_plink ;;
    "bcftools") install_bcftools ;;
    "bgenix")   install_bgenix ;;
    *) echo 'Invalid argument. Please provide "bgenix", "bcftools", or "plink".'
       exit 1 ;;
  esac
fi

