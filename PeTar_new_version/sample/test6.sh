
 

mcluster -M 14000 -m 0.7 -P 0 -R 0.457 -f 0 -u 1 -C 5 >mc.log


# use petar.init to create initial data for petar.
# the mcluster option '-u 1' generate data in astronomical unit (Msun, pc, km/s), but petar requires a self-consistent unit of velocity: pc/Myr, '-v kms2pcmyr' will do this.
petar.init -s bse -v kms2pcmyr -f input test.dat.10
# Use '-t 100.0' to run the simulation for 100 Myr.
# Use '-o 5' to generate output snapshots every 5 Myr.
# Use '-u 1' to set the units to astronomical units (Msun, pc, pc/Myr).
# Use '-b 500' to specify the number of primordial binaries as 500.
# By default, OpenMP utilizes all CPU threads. If you wish to use a specific number of threads, please add 'OMP_NUM_THREADS=[number of threads]'."
#OMP_STACKSIZE=128M petar -u 1 -b 500 -t 5.0 -o 0.05 input &>output
#OMP_STACKSIZE=128M OMP_NUM_THREADS=12 mpiexec -n 4 petar -u 1 --bse-t 50 -o 0.1 input_bh &> output
# after mode finished, gether the output data and do post-data process to detect binaries, obtain Lagrangian and core radii and corresponding properties.
# To maintain consistent units during post-processing, use '-G 0.00449830997959438' to set the gravitational constant to astronomical units.
#petar.data.gether data
#petar.data.process -G 0.00449830997959438 data.snap.lst
nohup bash -c "

# MPI also creates pthreads; OMP_STACKSIZE only covers OpenMP workers.
# BSE contributes about 27 MiB of static thread-local storage.
ulimit -S -s 131072 || exit 1
export UCX_VFS_ENABLE=n
export OMP_STACKSIZE=128M
export OMP_NUM_THREADS=4
mpirun --mca btl ^openib -np 4 petar --debug_lessmerger -a 0 -u 1 --stellar-evolution 2 \
--energy-err-hard 1e-3 \
--frozen-radius-factor 1.5 \
-t 500 -o 1 input &> output.log &&
petar.data.gether data &&
petar.data.process -G 0.00449830997959438 data.snap.lst
" &
