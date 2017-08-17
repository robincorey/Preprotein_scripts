#!/bin/bash
# not necessarily working :q

pdb2gmx -f ${NAME} -o ${NAME}.gro -p ${NAME}.top -i ${NAME}_posre.itp -water spc -ff oplsaa >& pdb2gmxout
genbox -cp ${NAME}.gro -o ${NAME}_sol.gro -p ${NAME}.top -cs
${GMX51} grompp -f ${SETUP}/for_genion.mdp -c ${NAME}_sol.gro -p ${NAME}.top -o forgenion >& grompp1
echo SOL | ${GMX51} genion -s forgenion -p ${NAME}.top -o ions.gro  -pname NA -nname CL -neutral >& genion1

if [[ ! -f ions.gro ]]; then
	echo genion failed
	exit 0
fi

echo -e "Waters | IONS" '\n' q | make_ndx -f ions.gro -o groups.ndx >& mnd1

newna=`grep -c NA all_ions.gro`
newcl=`grep -c CL all_ions.gro`
newsol=`grep -c OW all_ions.gro`

#sed "s/SOL /SOL ${newsol} ;/g" newtop.top -i
#sed "s/NA /NA ${newna} ;/g" newtop.top -i
#sed "s/CL /CL ${newcl} ;/g" newtop.top -i

${GMX51} grompp -f ${SETUP}/em.mdp -c ions.gro -p ${NAME}.top -o ${NAME}_em -n groups.ndx
${GMX51} mdrun -v -deffnm ${NAME}_em
${GMX51} grompp -f ${SETUP}/pr_GPU_5.1_fold.mdp -c ${NAME}_em -p ${NAME}.top -o ${NAME}_pr -n groups.ndx
${GMX51} mdrun -v -deffnm ${NAME}_pr

${GMX51} grompp -f ../md_GPU_5.1.mdp -c ${NAME}_pr -p ${NAME}.top -o ${NAME}_md -n groups.ndx
