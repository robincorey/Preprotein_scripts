#!/bin/bash

# To resize a simulation box 
# Currently only good for POPC, SPC etc. Easy to make more robust, if anyone cared.

SYSTEM=          # Name of tpr file
COORD=ATP.gro                   # Name of gro file
PNAME=P8			# ID of Phos in gro file (e.g. P8)
TOP=ATP.top 			# Name of toplogy file
IONS=y				# y,n

# first, make selection.dat

x1=0    #
x2=12   # Define atoms to keep
y1=0    # innit
y2=12   #
z1=0
z2=6

# Get protein molecule

echo Protein | editconf -f ${COORD}.gro -n ${COORD}.ndx -o resizeprotein.gro -bt triclinic

## make index for building individual ndxs

echo q '\n' | make_ndx -f ${COORD}.gro -o ${COORD}.ndx
cp ${COORD}.ndx ${COORD}_count.ndx
groups=`grep -c "\[" ${COORD}_count.ndx`

## prepare temporary toplogy file

cp ${TOP}.top del.top
sed 's/Protein        /;Protein        /g' del.top -i 

## Prepare for main building of system
## For each vector (x1,x2,y1,y2,z1,z2) follow this process:
## 1. Create selection.dat
## 2. Run g_select and append to original index file
## 3. Create box according to reduced coordinates
## 4. Update toplogy
## 5. Create simple tpr as input for next vector

rm del*.gro
rm resizesolemem.gro

#x1
echo -e 'resname "POP|SOL|NA|CL" and x>'${x1}'' > selection_one.dat
g_select -f ${COORD}.gro -s ${SYSTEM} -sf selection_one.dat -on ${SYSTEM}_one -selrpos whole_res_com -rmpbc
cat ${COORD}.ndx ${SYSTEM}_one.ndx > one.ndx
echo ${groups} | editconf -f ${COORD}.gro -n one.ndx -o del1.gro -bt triclinic
sed '/POP/d' del.top -i; sed '/SOL/d' del.top -i; sed '/NA/d' del.top -i; sed '/CL/d' del.top -i;
echo "POPC `grep -c ${PNAME} del1.gro`" >> del.top
echo "SOL `grep -c OW del1.gro`" >> del.top
if [[ ${IONS} == "y" ]]; then
	echo "NA `grep -c NA del1.gro`" >> del.top
	echo "CL `grep -c CL del1.gro`" >> del.top
fi
grompp -f ~/Desktop/MD_KIT/em1.mdp -p del.top -c del1.gro -o del1.tpr -maxwarn 5

### checkpoint

if [[ ! -f del1.tpr ]]; then 
	echo error innit del1
	exit 0
fi

#x2
echo -e 'resname "POP|SOL|NA|CL" and x<'${x2}'' > selection_two.dat
g_select -f del1.gro -s del1.tpr -sf selection_two.dat -on ${SYSTEM}_two -selrpos whole_res_com -rmpbc
cat ${COORD}.ndx ${SYSTEM}_two.ndx > two.ndx
echo ${groups} | editconf -f del1.gro -n two.ndx -o del2.gro -bt triclinic
sed '/POP/d' del.top -i; sed '/SOL/d' del.top -i; sed '/NA/d' del.top -i; sed '/CL/d' del.top -i;
echo "POPC `grep -c ${PNAME} del2.gro`" >> del.top
echo "SOL `grep -c OW del2.gro`" >> del.top
if [[ ${IONS} == "y" ]]; then
        echo "NA `grep -c NA del1.gro`" >> del.top
        echo "CL `grep -c CL del1.gro`" >> del.top
fi
grompp -f ~/Desktop/MD_KIT/em1.mdp -p del.top -c del2.gro -o del2.tpr -maxwarn 5

### checkpoint

if [[ ! -f del2.tpr ]]; then
        echo error innit del2
        exit 0
fi

##y1
echo -e 'resname "POP|SOL|NA|CL" and y>'${y1}'' > selection_three.dat
g_select -f del2.gro -s del2.tpr -sf selection_three.dat -on ${SYSTEM}_three -selrpos whole_res_com -rmpbc
cat ${COORD}.ndx ${SYSTEM}_three.ndx > three.ndx
echo ${groups} | editconf -f del2.gro -n three.ndx -o del3.gro -bt triclinic
sed '/POP/d' del.top -i; sed '/SOL/d' del.top -i; sed '/NA/d' del.top -i; sed '/CL/d' del.top -i;
echo "POPC `grep -c ${PNAME} del3.gro`" >> del.top
echo "SOL `grep -c OW del3.gro`" >> del.top
if [[ ${IONS} == "y" ]]; then
        echo "NA `grep -c NA del1.gro`" >> del.top
        echo "CL `grep -c CL del1.gro`" >> del.top
fi
grompp -f ~/Desktop/MD_KIT/em1.mdp -p del.top -c del3.gro -o del3.tpr -maxwarn 5

### checkpoint

if [[ ! -f del3.tpr ]]; then
        echo error innit del3
        exit 0
fi

#y2
echo -e 'resname "POP|SOL|NA|CL" and y<'${y2}'' > selection_four.dat
g_select -f del3.gro -s del3.tpr -sf selection_four.dat -on ${SYSTEM}_four -selrpos whole_res_com -rmpbc
cat ${COORD}.ndx ${SYSTEM}_four.ndx > four.ndx
echo ${groups} | editconf -f del3.gro -n four.ndx -o del4.gro -bt triclinic
sed '/POP/d' del.top -i; sed '/SOL/d' del.top -i; sed '/NA/d' del.top -i; sed '/CL/d' del.top -i;
echo "POPC `grep -c ${PNAME} del4.gro`" >> del.top
echo "SOL `grep -c OW del4.gro`" >> del.top
if [[ ${IONS} == "y" ]]; then
        echo "NA `grep -c NA del1.gro`" >> del.top
        echo "CL `grep -c CL del1.gro`" >> del.top
fi
grompp -f ~/Desktop/MD_KIT/em1.mdp -p del.top -c del4.gro -o del4.tpr -maxwarn 5

### checkpoint

if [[ ! -f del4.tpr ]]; then
        echo error innit del4
        exit 0
fi

#z1
echo -e 'resname "POP|SOL|NA|CL" and z>'${z1}'' > selection_five.dat
g_select -f del4.gro -s del4.tpr -sf selection_five.dat -on ${SYSTEM}_five -selrpos whole_res_com -rmpbc
cat ${COORD}.ndx ${SYSTEM}_five.ndx > five.ndx
echo ${groups} | editconf -f del4.gro -n five.ndx -o del5.gro -bt triclinic
sed '/POP/d' del.top -i; sed '/SOL/d' del.top -i; sed '/NA/d' del.top -i; sed '/CL/d' del.top -i;
echo "POPC `grep -c ${PNAME} del5.gro`" >> del.top
echo "SOL `grep -c OW del5.gro`" >> del.top
if [[ ${IONS} == "y" ]]; then
        echo "NA `grep -c NA del1.gro`" >> del.top
        echo "CL `grep -c CL del1.gro`" >> del.top
fi
grompp -f ~/Desktop/MD_KIT/em1.mdp -p del.top -c del5.gro -o del5.tpr -maxwarn 5

### checkpoint

if [[ ! -f del5.tpr ]]; then
        echo error innit del5
        exit 0
fi

#z2
echo -e 'resname "POP|SOL|NA|CL" and z<'${z2}'' > selection_six.dat
g_select -f del5.gro -s del5.tpr -sf selection_six.dat -on ${SYSTEM}_six -selrpos whole_res_com -rmpbc
cat ${COORD}.ndx ${SYSTEM}_six.ndx > six.ndx
echo ${groups} | editconf -f del5.gro -n six.ndx -o del6.gro -bt triclinic
sed '/POP/d' del.top -i; sed '/SOL/d' del.top -i; sed '/NA/d' del.top -i; sed '/CL/d' del.top -i;
echo "POPC `grep -c ${PNAME} del6.gro`" >> del.top
echo "SOL `grep -c OW del6.gro`" >> del.top
if [[ ${IONS} == "y" ]]; then
        echo "NA `grep -c NA del1.gro`" >> del.top
        echo "CL `grep -c CL del1.gro`" >> del.top
fi

## Resizing done!

## Check

if [[ ! -f del6.gro ]]; then echo Resizing has not worked. Check output files! ; exit 0; fi

## Combine into one gro file (with pdb out)

rm resizeall.gro
tail -n+3 del6.gro > del6_crop.gro
tail -n+3 resizeprotein.gro > delprot1.gro
sed '$ d' delprot1.gro > delprot2.gro
cat delprot2.gro del6_crop.gro >> resizeall.gro
atoms=`grep -c [0-9] resizeall.gro`; atoms2=`echo ${atoms} - 1 | bc` ### A bit ropey for some reason...
sed '1i '${atoms2}'' resizeall.gro -i
sed '1i Resized box' resizeall.gro -i
editconf -f resizeall.gro -o resizeall.pdb
