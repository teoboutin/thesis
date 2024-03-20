#!/bin/bash
#~ pwd


THESIS_DIR=$(dirname $(realpath -s $0))
echo "Installing thesis template from path $THESIS_DIR"



LYXLAYOUTDIR=$HOME/.lyx/layouts

echo "Installing lyx files in ${LYXLAYOUTDIR}"
rm $LYXLAYOUTDIR/mythesis.layout
ln -s ${THESIS_DIR}/template/lyx/mythesis.layout $LYXLAYOUTDIR/mythesis.layout


TEXDIR=$HOME/texmf
echo "Installing latex files in ${TEXDIR}"

LINK_NAME_LATEX=$TEXDIR/tex/latex/thesis
rm $LINK_NAME_LATEX
ln -s ${THESIS_DIR}/template/latex $LINK_NAME_LATEX
texhash $TEXDIR



echo "Preparing symlinks to macros and bib"
DIRLIST="chapter1 chapter2 chapter3 chapter4 annex_asymptotics annex_lbm annex_solutions_sharp_interface"

for DIR in ${DIRLIST}
do
    ln -sf "${THESIS_DIR}/macros.lyx" "${THESIS_DIR}/${DIR}/macros.lyx"
    ln -sf "${THESIS_DIR}/thesis.bib" "${THESIS_DIR}/${DIR}/thesis.bib"
done

### prepare textidote run
echo "Preparing files for running textidote"
DIRLIST="chapter1 chapter2 chapter3 chapter4 ccl"

script=$(cat ./.textidote)

for DIR in ${DIRLIST}
do
    local_script=${script//dummy/${DIR}}
    echo $local_script > "${THESIS_DIR}/${DIR}/.textidote"
    ln -sf "${THESIS_DIR}/dico.txt" "${THESIS_DIR}/${DIR}/dico.txt"
done

### create symlinks to thesis chapter figures in defense directory
echo "Preparing figures for defense"

fig_path="defense/Figures/main"
mkdir -p ${fig_path}
DIRLIST="chapter1 chapter2 chapter3 chapter4"
for DIR in ${DIRLIST}
do
    rm -f "${THESIS_DIR}/${fig_path}/${DIR}"
    ln -s "${THESIS_DIR}/${DIR}/Figures" "${THESIS_DIR}/${fig_path}/${DIR}"
done



### build all tikz figures
list=$(find -wholename *Tikz*.tex)
for file in ${list}
do
	echo "${file}"
    r=$(cd $(dirname $file) && ${THESIS_DIR}/make_tikz_figure.py $(basename $file))
done
