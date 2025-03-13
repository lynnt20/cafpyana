# Setup grid submission

outDir=$1
echo "@@ outDir : ${outDir}"
DFPREFIX=$2

source /cvmfs/larsoft.opensciencegrid.org/spack-packages/setup-env.sh

echo "@@ pwd"
pwd
echo "@@ ls -alh"
ls -alh
echo "@@ git clone cafpyana"
git clone https://github.com/sungbinoh/cafpyana.git
echo "@@ cd to cafpyana dir"
cd cafpyana
echo "@@ ls -alh"
ls -alh
echo "@@ check if there is cmake"
spack find cmake
spack load cmake@3.27.7
which cmake
echo "@@ check if other spack packages"
spack find hdf5
spack find xrootd
echo "@@ check if there is util-linux"
rpm -q util-linux 
echo "@@ run init.sh"
source init.sh
echo "@@ ls -alh"
ls -alh
echo "@@ mkdir output"
mkdir output
echo "@@ Done!"
thisOutputCreationDir=`pwd`
filesFromSender=${CONDOR_DIR_INPUT}/bin_dir/

export IFDH_CP_MAXRETRIES=2

echo "@@ outDir : "${outDir}
echo "@@ ifdh  mkdir_p "${outDir}
ifdh  mkdir_p ${outDir}

echo "@@ source ${filesFromSender}/run_"${nProcess}".sh "
source ${filesFromSender}/run_${nProcess}.sh &> log_${nProcess}.log
echo "@@ Check output : ${DFPREFIX}_${nProcess}.df"
ls -alh ${DFPREFIX}_${nProcess}.df

outFILE=${thisOutputCreationDir}/${DFPREFIX}_${nProcess}.df
if [ -f "$outFILE" ]; then
  echo "ifdh cp ${thisOutputCreationDir}/${DFPREFIX}_${nProcess}.df ${outDir}/${DFPREFIX}_${nProcess}.df"
  ifdh cp ${thisOutputCreationDir}/${DFPREFIX}_${nProcess}.df ${outDir}/${DFPREFIX}_${nProcess}.df
  echo "ifdh cp ${thisOutputCreationDir}/log_${nProcess}.log ${outDir}/log_${nProcess}.log"
  ifdh cp ${thisOutputCreationDir}/log_${nProcess}.log ${outDir}/log_${nProcess}.log
  echo "@@ Done!"
else
  echo "File not exist"
fi

