#!/bin/bash
# run mutations of preprotein using scwrl
# this is specific for Tom AYEG with the snip subs thing

####### Changes
### Adaped for full subs on 22/02/2017
### Adapted for three timepoints on 27/02/2017


### Just put standard stuff here. Will sort it out below

DESCIN=RQH_full		
GROIN=GGG	 ### Must match the name of input GRO. e.g. GGG.110ns.gro >> GGG
TOPIN=GGG_full.top	


# Don't change owt
###############################################

# Ok. Let's set the loop up.

TIME=(70 80 110)   # always
CD=`pwd`

cd $CD

for l in `seq 0 2`; do
	# Tidy all variables
	DESC=${DESCIN}_${TIME[l]}
	GRO2=${GROIN}.${TIME[l]}
	if [[ ! -f ${GROIN}.${TIME[l]}.gro ]]; then
		echo ${GROIN}.${TIME[l]}.gro not found
		exit 0
	fi
	TOP2=${TOPIN//.top}
	editconf -f ${GRO2}.gro -o ${GRO2}.pdb >& edc
	PDB=${GRO2}.pdb
# 1. make input files

if [[ ! -d ${DESC} ]]; then
        mkdir $DESC
fi

cd $DESC
cp ../*itp .

if [[ -f mutateto ]]; then
	rm mutateto
fi

for i in `seq 0 $len2`; do
	echo -n ${NEW[i]} >> mutateto
done

# make one big seq input

grep CA ../${PDB} | awk '{print $4}' | tr '\n' ' ' | sed 's/ALA/A/g;s/CYS/C/g;s/ASP/D/g;s/GLU/E/g;s/PHE/F/g;s/GLY/G/g;s/HIS/H/g;s/ILE/I/g;s/LYS/K/g;s/LEU/L/g;s/MET/M/g;s/ASN/N/g;s/PRO/P/g;s/GLN/Q/g;s/ARG/R/g;s/SER/S/g;s/THR/T/g;s/VAL/V/g;s/TRP/W/g;s/TYR/Y/g' | sed 's/ //g' > ${PDB}.seq

grep CA ../tip.pdb | awk '{print $4}' | tr '\n' ' ' | sed 's/ALA/A/g;s/CYS/C/g;s/ASP/D/g;s/GLU/E/g;s/PHE/F/g;s/GLY/G/g;s/HIS/H/g;s/ILE/I/g;s/LYS/K/g;s/LEU/L/g;s/MET/M/g;s/ASN/N/g;s/PRO/P/g;s/GLN/Q/g;s/ARG/R/g;s/SER/S/g;s/THR/T/g;s/VAL/V/g;s/TRP/W/g;s/TYR/Y/g' | sed 's/ //g' > tip.seq

newtrim=`cat mutateto`
sed "s/.\{$len3\}$/-replace-/" ${PDB}.seq -i
sed "s/-replace-/${newtrim}/" ${PDB}.seq -i
cat tip.seq >> ${PDB}.seq 

grep -v "POP\|NA\|CL\|MG\|SOL\|ATP\|ADP" ../${PDB} > input.pdb

sed 's/O1/O /g' input.pdb -i
#sed 's/O2/OT2/g' input.pdb -i

# 2. do scwrl and make new gro

if [[ -f mutated.pdb ]]; then
	rm mutated.pdb
fi

Scwrl4 -i input.pdb -t -s ${PDB}.seq -o mutated.pdb >& scrwloutput.txt ## do in conext of ayeg to provent horros! 

sed '/HD1 HIS/d' mutated.pdb -i

if [[ ! -f mutated.pdb ]]; then
	echo scrwl failed
	exit 0
else
	echo Scwrl worked!
fi

if [[ -f xx01.pdb ]]; then
	rm xx01.pdb xx00.pdb
fi

#Trim into two pdbs
#split=`grep "N   LYS  1231" mutated.pdb | awk '{print $2}'`
#split2=`echo "$split - 6" | bc`
#
#csplit mutated.pdb "$split2" # makes xx00 and xx01
#

csplit -k mutated.pdb '/N   LYS  1231/'

# 3. make new topology       

mv xx01 xx01.pdb
mv xx00 xx00.pdb

pdb2gmx -f xx01 -o subs_mut.gro -p subs_mut.top -i subs_mut_posre.itp -water spc -ff oplsaa >& pdb2gmxout

cp ../${TOP2}.top newtop.top

mv subs_mut.top subs_mut.itp
sed '/forcefield/d' subs_mut.itp -i
sed '/water topology/,$d' subs_mut.itp -i
sed 's/Protein/Protein5/g' subs_mut.itp -i 
#itpfile=`grep -A1 SS.itp newtop.top |tail -n 1`
#sed "s/${itpfile}/#include \"subs_mut.itp\"/g" newtop.top -i
sed "s/model_adj_mv_Protein5.itp/subs_mut.itp/g" newtop.top -i

#so fucking fiddly!!

sed '/N   ASP   757/i \TER \' xx00.pdb -i
sed '/N   LYS  1175/i \TER \' xx00.pdb -i
sed '/N   LYS  1231/i \TER \' xx00.pdb -i 
sed '/12141  N/i \TER \' xx00.pdb -i

if [[ -f prot.gro ]]; then
	rm prot.gro
fi

pdb2gmx -f xx00.pdb -o prot.gro -p del.top -i del -ignh -water spc -ff oplsaa >& pdb2gmxout2
his=`grep -c "CA  HIS" xx00.pdb`
hl=`echo ${#his[@]} -1 | bc`

#rm protonation
#
#for l in `seq 0 $his`; do
#	echo 1 >> protonation
#done

pdb2gmx -f xx00.pdb -o prot.gro -p del.top -i del -water spc -ff oplsaa -his -ignh  < ../HISPROT >& pdb2gmxout2
# 4. combine

#276HIS    HD2
tail -n +3 prot.gro | head -n -1 > prot2.gro
egrep -e MG -e ATP -e ADP ../${GRO2}.gro > nuc.gro
egrep -e POP -e SOL -e NA -e CL ../${GRO2}.gro > sol.gro
tail -n +3 subs_mut.gro | head -n -1 > subs_mut2.gro
cat nuc.gro prot2.gro subs_mut2.gro sol.gro > all.gro
count=`grep -c [0-9] all.gro`
sed "1 i\ ${count}" all.gro -i
sed "1 i\ ${DESC}" all.gro -i
tail -n 1 ../${GRO2}.gro >> all.gro

#5 prepare
echo -e "2 | 3 | 7" '\n' "2 | 3 | 7 | 4" '\n' q | make_ndx -f all.gro -o groups.ndx >& md1

newpop=`grep -c P8 all.gro`
newsol=`grep -c OW all.gro`
newna=`grep -c NA all.gro`
newcl=`grep -c CL all.gro`

sed '/POP/d' newtop.top -i
sed '/SOL/d' newtop.top -i
sed '/NA/d' newtop.top -i
sed '/CL/d' newtop.top -i

sed '/^Protein5/a \POP \' newtop.top -i
sed '/POP/a \SOL \' newtop.top -i
sed '/SOL/a \NA \' newtop.top -i
sed '/NA/a \CL  \' newtop.top -i

sed "s/POP/POP   ${newpop}/" newtop.top -i
sed "s/SOL/SOL   ${newsol}/" newtop.top -i
sed "s/NA/NA   ${newna}/" newtop.top -i
sed "s/CL/CL   ${newcl}/" newtop.top -i

grep NA all.gro > newna.gro
grep CL all.gro > newcl.gro
cat newna.gro newcl.gro > newions.gro
tail -n 1 all.gro >> newions.gro
sed '/NA/d' all.gro -i
sed '/CL/d' all.gro -i
head -n -1 all.gro > pre2.gro
cat pre2.gro newions.gro > all2.gro

${GMX51} grompp -f ../pr_GPU_5.1.mdp -c all2.gro -p newtop.top -o ${DESC} -n groups.ndx >& grompp1

echo SOL | ${GMX51} genion -s ${DESC} -o all_ions.gro  -pname NA -nname CL -neutral >& genion1

if [[ ! -f all_ions.gro ]]; then
	echo genion failed
	exit 0
fi

echo -e "2 | 3 | 7" '\n' "2 | 3 | 7 | 4" '\n' "4 | 27 " '\n' q | make_ndx -f all_ions.gro -o groups.ndx >& mnd1
newna=`grep -c NA all_ions.gro`
newcl=`grep -c CL all_ions.gro`
newsol=`grep -c OW all_ions.gro`

sed "s/SOL /SOL ${newsol} ;/g" newtop.top -i
sed "s/NA /NA ${newna} ;/g" newtop.top -i
sed "s/CL /CL ${newcl} ;/g" newtop.top -i

# 6 prevent overlapping atoms

grep NA all_ions.gro > newions_na.gro
grep CL all_ions.gro > newions_cl.gro
tail -n 1 all_ions.gro >> newions_cl.gro
sed '/NA/d' all_ions.gro -i
sed '/CL/d' all_ions.gro -i
head -n -1 all_ions.gro > pre_a.gro
cat pre_a.gro newions_na.gro newions_cl.gro > all_ions_a.gro


sed 's/Nuc/MG_ATP_Protein/g' ~/Desktop/MD_KIT/gmembed_nuc.mdp > gmembed_preprotein.mdp
sed 's/Mem/POP/g' gmembed_preprotein.mdp -i 
sed 's/SOL/Water_and_ions/g' gmembed_preprotein.mdp -i

if [[ -f gmembed.tpr ]]; then
	rm gmembed.tpr
fi

#head -n -5 newtop.top > temp.top    ## gromacs 4 don't like no dist res
grompp -v -f ../gmembed_preprotein.mdp -o gmembed -c all_ions_a.gro -p newtop.top -n groups -maxwarn 2 >& grompp2


# This removes overlapping atoms, removing a few lipids (but not too many) and some solvent
echo MG_ATP_Protein Water_and_ions | mdrun -s gmembed.tpr -membed ~/Desktop/MD_KIT/membed.dat -o traj.trr -c membedded.gro -e ener.edr -cpt -1 -v -stepout 100 -mn groups 
#need to update top again
newpop=`grep -c P8 membedded.gro`
newsol=`grep -c OW membedded.gro`
newna=`grep -c NA membedded.gro`
newcl=`grep -c CL membedded.gro`

if [[ ! -f membedded.gro ]]; then
	echo membed fail
	exit 0
fi

sed '/POP/d' newtop.top -i
sed '/SOL/d' newtop.top -i
sed '/NA/d' newtop.top -i
sed '/CL/d' newtop.top -i

sed '/Protein5/a \ POP \' newtop.top -i 
sed '/POP/a \ SOL \' newtop.top -i
sed '/SOL/a \ NA \' newtop.top -i
sed '/NA/a \ CL  \' newtop.top -i

sed "s/POP/POP   ${newpop}/" newtop.top -i
sed "s/SOL/SOL   ${newsol}/" newtop.top -i
sed "s/NA/NA   ${newna}/" newtop.top -i
sed "s/CL/CL   ${newcl}/" newtop.top -i

egrep -e NA -e CL membedded.gro > newions.gro
tail -n 1 membedded.gro >> newions.gro
sed '/NA/d' membedded.gro -i
sed '/CL/d' membedded.gro -i
head -n -1 membedded.gro > pre2.gro

cat pre2.gro newions.gro > BUILT.gro

sed '/CA/d' newtop.top > ${DESC}.top

#cp ../*itp .

#dist1=`grep "^ 1268" BUILT.gro | grep CA | awk '{print $2}' | tr -d CA`
#dist2=`grep "^  822" BUILT.gro | grep CA | awk '{print $2}' | tr -d CA`
#dist3=`grep "^ 1269" BUILT.gro | grep CA | awk '{print $2}' | tr -d CA`
#dist4=`grep "^ 1149" BUILT.gro | grep CA | awk '{print $2}' | tr -d CA`
#
#echo "
#$dist2 $dist1   1    0    1    0.0   0.3   0.6  1.0     ; ILE822 CA to X1268 CA
#$dist4 $dist3   1    1    1    0.0   0.8   1.0  1.0     ; ILE1149 CA to X1269 CA
#
#" >> ${DESC}.top

echo -e "2 | 3 | 7" '\n' "2 | 3 | 7 | 4" '\n' "4 | 27 " '\n' q | make_ndx -f BUILT.gro  -o groups.ndx
${GMX51} grompp -f ../em.mdp -c BUILT.gro -p ${DESC}.top -o ${DESC}_em -n groups.ndx
${GMX51} mdrun -v -deffnm ${DESC}_em

${GMX51} grompp -f ../pr_GPU_5.1.mdp -c ${DESC}_em -p ${DESC}.top -o ${DESC}_pr -n groups.ndx
${GMX51} mdrun -v -deffnm ${DESC}_pr

done
