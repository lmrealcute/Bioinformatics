
cd /Users/chayahboyd/Documents/GitHub/Bioinformatics/RProjects/HW7Phyl


/Applications/raxml-ng_v1.2.1_macos_M1/raxml-ng \
--all \
--msa primate.phy \
--model GTR+G \
--prefix T1 \
--tree pars{10} \
--bs-trees 200 \
--data-type DNA \
--msa-format PHYLIP

