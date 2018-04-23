set -e

mkdir -p dataset-levure
(
    cd dataset-levure &&
    wget -O Snf2_countdata.tar.gz "https://github.com/bartongroup/profDGE48/blob/master/Preprocessed_data/Snf2_countdata.tar.gz?raw=true" &&
    wget -O WT_countdata.tar.gz "https://github.com/bartongroup/profDGE48/blob/master/Preprocessed_data/WT_countdata.tar.gz?raw=true" &&
    tar xvfz Snf2_countdata.tar.gz &&
    tar xvfz WT_countdata.tar.gz &&
    Rscript ../dataset_levure.R &&
    rm *.gbgout *.tar.gz
)
