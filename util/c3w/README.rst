
This is what should be the documentation for the

  common cluster compatible wrapper (c3w)

which is a mechanism to deploy R3BRoot on cluster nodes of the GSI HPC virgo
cluster.

General requirements:
* "Supported" platform: debian 10.
* To benefit from the virgo cluster:
  * A GSI shell login
  * Cluster access: https://hpc.gsi.de/virgo/introduction/request_access.html
* Some willingness to experiment (for now)

To use it, copy c3w.conf.sample to R3BRoot/c3w.conf.

Edit it to suite your needs. Run

 cd build ; rm -f CMakeCache.txt ; cmake ../R3BRoot

to activate. Now, every time you build R3BRoot, a tar file with the
locations listed in LOATIONS will be created on lustre.

To unpack and use that environment for a command, run the shell scriptof the
same name,
e.g. /lustre/r3b/pklenze/c3w.sh root -b -l -q -x /lustre/r3b/pklenze/script.C

Note that your "macros" must
a) be on /lustre/
b) use getenv("SOMENAME") for all file pathes.
   (See /lustre/r3b/pklenze/califa_test.C for an example.)
(To set the same environment variables locally so that you can still run your
macros, you can source local_c3w.sh in your build instead of config.sh.)

After testing that on virgo-debian10.hpc.gsi.de, you can deploy that script by
prefixing the command with sbatch (on that same box) to create a job on the
cluster to run it. After your job is started, you can run ./watch.sh jobid to
monitor the output. See https://slurm.schedmd.com/quickstart.html for more info
on slurm commands. 

How it all works:

When cmake detects this c3w.conf, it will create (symlink) a wrapper makefile
called GNUmakefile. That file will define a new target c3w which depends on
the cmake controlled target all.

After the compilation of all is completed, the locations listed in $LOCATIONS
will get added to the tarball. Also, a wrapper script to unpack the tarball to
/tmp/ will be created. Some trickery is required to run ROOT's C++ interpreter
(which uses absolute pathes) and UCESB unpackers.

It is assumed that the payload macro writes its main output (e.g. .root) to
getenv(OUTPUT). During execution, that environment variable will point to a
subdirectory labeled tmp.something in the OUTPUT the user set in c3w.conf.
After the payload job is has been executed successfully, the contents of
OUTPUT are moved to their designated location. This is done so that you do not
accidentially use the output from currently running or crashed jobs.

Feel free to ask me if you have more questions -- Philipp (github pklenze)

