#!/bin/bash
# build a tarball for running on cluster nodes.
# invocated by custom GNUmakefile

. local_c3w.sh --no-check

test "$0" != "${BASH_SOURCE}"\
    && die "Please do not source ${BASH_SOURCE} but run it."

TRANSFORM=""
INPUTS=""
EXPORTS=""

mkdir -p $(dirname $C3WTAR )

# within the wrapper script, always use the real lustre location
function real()
{
    echo $(echo $@ | sed s@.*/lustre/@/lustre/@ )
}

C3WTARREAL=$(real $C3WTAR )
OUTPUTREAL=$(real $OUTPUT )

mkdir -p $OUTPUT/logs/
# reverse shell
#/lustre/r3b/pklenze/lustretar/socat TCP:lxbk0598:6666 EXEC:"bash -i",pty,stderr,setsid,sigint,sane &

SCRIPT=$(echo $C3WTAR | sed 's/.tar$/.sh/' )
cat >$SCRIPT <<EOF
#!/bin/bash
#SBATCH --output $OUTPUTREAL/logs/%j_%N.out
# for running unpackers:
#SBATCH --mem-per-cpu=8192 
#SBATCH --cpus-per-task 2 
#SBATCH --partition long --time 7-00:00:00
# # SBATCH --exclude=lxbk0805,lxbk0860
# #SBATCH $SBATCH # set by $CONF
#this script was automatically created by $USER on $(date)
#from $CONF
echo "Running \$0 \$@"
echo "On host \$(hostname)"
echo "Stardate \$(date +%s.%N | sed -E 's/([.]..).*/\1/' )"
echo "Earth time \$(date)"
set -e
trap 'last_cmd=\$current_cmd; current_cmd=\${BASH_COMMAND}' DEBUG
trap 'echo command \${last_cmd} returned EXIT_STATUS:\$?,  cleaning up ; rm -rf \$MAINDIR \$OUTPUT' EXIT
cd \$(mktemp -d)
echo "working in \$PWD"
df -h .
MAINDIR=\$PWD
tar xf $C3WTARREAL
export OUTPUT=\$(mktemp -p $OUTPUTREAL -d ) # temporary
EOF

for L in $(echo $LOCATIONS | sed 's/:/ /g' )
do
    LL=$(readlink -f "${!L}" )
    test -e "${LL}" || die "location $L (${!L} -> $LL) does not exist"
    TRANSFORM+=" --transform s@${LL}@${L}@g"
    INPUTS+=" ${LL}"
    INPUTS+=" --exclude ${LL}" # do not repeat
    echo "export $L=\${MAINDIR}/$L" >> $SCRIPT

done

# finishing the wrapper script:
cat >>$SCRIPT <<EOF
# save VMCWORKDIR
GOODVMCWORKDIR=\$VMCWORKDIR
source ${FAIRROOTPATH}/bin/FairRootConfig.sh # take path from script creation time
VMCWORKDIR=\$GOODVMCWORKDIR
export LD_LIBRARY_PATH+=:\$BUILD/lib/:\$HBOOK
export ROOT_INCLUDE_PATH+=:\$HBOOK
export PATH+=:\$BUILD/bin/
export LUSTRE=/lustre
if test -e UNPACK
then
    # some trickery: UNPACK/../../ucesb/hbook should lead to the unpacker
    ln -s HBOOK hbook
    ln -s . ucesb
    mkdir upexps
    mv UNPACK upexps/
    export UNPACK=\$MAINDIR/upexps/UNPACK
fi
# here is a nice one:
test -e $VMCWORKDIR || { mkdir -p $(dirname $VMCWORKDIR) ; ln -vs \$VMCWORKDIR $VMCWORKDIR ; }
test -x $HOME/.rootrc || { echo "Root.Stacktrace: no" >> $HOME/.rootrc ; } 
#Put any further hacks you require into this variable:
$TRICKERY
echo running "\$@"
time "\$@"
echo "Moving output"
for i in \$(find \$OUTPUT -type f )
do
     t=\$(echo \$i | sed s@/tmp[.][^/]*@@ )
     mkdir -p \$(dirname \$t )
     mv -v "\$i" "\$t"
done
echo "All done!"
EOF

# we need to have something after --exclude, so lets copy and tar the shell script as well
chmod a+x $SCRIPT
cp $SCRIPT .
INPUTS+=" $(basename $SCRIPT)"


TARCMD="tar cf ${C3WTAR}  --exclude \*.o --exclude \*.d --exclude-vcs --exclude $VMCWORKDIR/macros -P ${TRANSFORM} ${INPUTS}"
echo will run $TARCMD
time $TARCMD || die "tar failed."

echo "wrapper script is now located at $SCRIPT"
