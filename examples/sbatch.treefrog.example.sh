#!/bin/bash -l
# This script runs TreeFrog. It can be paired with the load balancing output produced by .loadbalance script.
# You must ensure that the numpermpi, nsteps and mpi match those in the loadbalance run (with nmpi here the nmpithreads in the
# loadbalance script.

# alter this script to point to the appropriate directories and also allocate the desired
# resources by changing stuff like JOBNAME, PROJECT, NCPUS, MEM, QUEUE, TIME
# note that this example script is designed to run treefrog with mpi threads but
# assumes one OpenMP thread per MPI process.

#set resources
#SBATCH -J=JOBNAME
#SBATCH -n=NCPUS
#SBATCH -t=TIME
#SBATCH --mem=MEM
#SBATCH -p=PARTITION
#SBATCH --account=ACCOUNT


nomp=NCPUS
nmpi=1
export OMP_NUM_THREADS=$nomp

# make sure to load modules need to run treefrog, like mpi, gsl, hdf5
module list

# treefrog exe
treefrog=/dir/to/TreeFrog/builddir/treefrog

# input dir where halos are located
inputdir=/dir/to/halos/

# diretory where trees will be written
outputdir=/dir/to/tree/
# name of tree
treename=DescendantTree
outfile=${outputdir}/${treename}

# base configuration file
treefrogbaseconfig=/dir/to/baseconfig

#alter the initial and final snaphots as desired
isnap=0
fsnap=100

# produce input file which is an ascii file listing halo catalogs
# from snapshots isnap to fsnap
inputfilelist=${outfile}.snaplist.txt
rm ${inputfilelist}
for ((i=${isnap};i<=${fsnap};i++))
do
    #this writes a list of file names assuming some convention. Alter if necessary.
	ii=`printf %03d $i`
	echo ${inputdir}/snapshot_${ii}.VELOCIraptor >> ${inputfilelist}
done


# alter config as desired. Here we assume that we wish to alter
# some specific things in base config file related to particle and halo IDs
treefrogconfig=${outfile}.treefrog.config
cp ${treefrogbaseconfig} ${treefrogconfig}
# we set the 18609625000 id and don't construct and id to index map
maxid=18609625000
map=0
temporalid=1000000000000
sed -i 's/MAXID/'"${maxid}"'/g' ${treefrogconfig}
sed -i 's/IDTOINDEXMAP/'"${map}"'/g' ${treefrogconfig}
sed -i 's/TEMPORALHALOIDVAL/'"${temporalid}"'/g' ${treefrogconfig}
sed -i 's/HALOOFFSET/0/g' ${treefrogconfig}
sed -i 's/SNAPSHOTOFFSET/'"${isnap}"'/g' ${treefrogconfig}

#run descendant tree
inputargs=" -i ${inputfilelist} -s $((${fsnap}-${isnap}+1)) "
mpirun -np ${nmpi} ${treefrog} ${inputargs} -o ${outfile} -C ${treefrogconfig} > ${outfile}.log
